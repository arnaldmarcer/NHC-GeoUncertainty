source("scripts/R/init.R")

#### FUNCTIONS -------------------------------------------------------------------------------- ####
# Get Australia mainland mask (only main Australian island)
getAustraliaMask <- function(model.taxon){
    f.out <- paste0(data.processed.dir, model.taxon$region.name, "_mask.tif")
    if(file.exists(f.out))
        return(terra::rast(f.out))
    aus.mpol <- gadm_sf_loadCountries("AUS", basefile = tmp.dir)$sf
    sf::st_crs(aus) <- 4326
    aus.pol <- st_cast(aus.mpol, "POLYGON")
    centroid <- st_centroid(aus)
    aus.boundary <- aus.pol[st_intersects(aus.pol, centroid, sparse = FALSE),][1] # We only want the mainland perimeter of Australia
    r <- terra::rast(paste0(data.raw.dir, "worldclim/wc2.1_30s_bio_1.tif"))
    r <- crop(r, model.taxon$region.ext)
    r.mask <- terra::rasterize(vect(aus.boundary), r)
    writeRaster(r.mask, filename = f.out, overwrite = T)
    cat("    ==>", f.out, "\n")
    return(r.mask)
}

# Get northern Canada and Greenland mask (oceans and lakes as NAs)
# Data source: https://www.worldwildlife.org/pages/global-lakes-and-wetlands-database
# Used Level 3 (GLWD-3), which "comprises lakes, reservoirs, rivers and different wetland types
# in the form of a global raster map at 30-second resolution."
getNorthernCanadaGreenlandMask <- function(model.taxon){
    f.out <- paste0(data.processed.dir, model.taxon$region.name, "_mask.tif")
    if(file.exists(f.out))
        return(raster(f.out))
    r.mask <- raster(paste0(data.processed.dir, "predictors/northern_canada_and_greenland/30s/bio1.tif"))
    f.glwd.tmp <- "tmp/glwd_lv3.zip"
    download.file(url="https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/9slil0ww7t_GLWD_level3.zip",
                  destfile=f.glwd.tmp)
    d.glwd.tmp <- paste0(tmp.dir, "glwd")
    dir.create(tmp.dir, recursive = T, showWarnings = F)
    unzip(zipfile = f.glwd.tmp, exdir=d.glwd.tmp)
    water.bodies <- raster(paste0(d.glwd.tmp, "/glwd_3"))
    raster::crs(water.bodies) <- "+init=epsg:4326"
    water.bodies <- raster::crop(water.bodies, model.taxon$region.ext)
    water.bodies[water.bodies <= 3] <- 1
    water.bodies[water.bodies > 3] <- NA
    water.bodies <- raster::projectRaster(water.bodies, r.mask, method = 'ngb')
    r.mask[water.bodies==1] <- NA
    r.mask[!is.na(r.mask)] <- 1
    raster::writeRaster(r.mask, filename = f.out, overwrite = T, options = c("COMPRESS=LZW"))
    return(r.mask)
    cat("    ==>", f.out, "\n")
}

# Prepares bioclim predictors
# Need to have worldclim-bioclim predictors at 30s resolution downloaded in data/raw/worldclim
# https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip
prepareBioclimPredictors <- function(model.taxon){
    dir.out <- paste0(data.processed.dir, "predictors/", model.taxon$region.name, "/", model.taxon$resolution , "/")
    dir.create(dir.out, recursive=T, showWarnings = F)
    for(i in c(1:19)){
        cat(paste0("  preparing predictor variable bio", i), "\n")
        f.in <- paste0(data.raw.dir, "worldclim/wc2.1_", model.taxon$resolution, "_bio_", i, ".tif")
        f.out <- paste0(dir.out, "bio", i, ".tif")
        cmd <- paste("gdal_translate -projwin",
                     model.taxon$region.ext@xmin, model.taxon$region.ext@ymax,
                     model.taxon$region.ext@xmax, model.taxon$region.ext@ymin,
                     f.in, f.out)
        system(cmd)
        unlink(paste0(dir.out, "bio", i, ".prj"))
        unlink(paste0(dir.out, "bio", i, ".asc.aux.xml"))
    }
    predictor.stack <- stack(list.files(dir.out, full.names = T, pattern = "*.tif"))
    return(predictor.stack)
}

# Apply mask to predictors
applyMask <- function(model.taxon){
    if(model.taxon$region.name == "australia"){
        r.mask <- getAustraliaMask(model.taxon)
    } else if(model.taxon$region.name == "northern_canada_and_greenland"){
        r.mask <- getNorthernCanadaGreenlandMask(model.taxon)
    }
    cat("  ")
    predictors.dir <- paste0(data.processed.dir, "predictors/", model.taxon$region.name, "/", model.taxon$resolution)
    tifs <- list.files(predictors.dir, full.names = T, pattern = "tif")
    for(t in tifs){
        cat("  -->", t, "\n")
        r <- raster(t)
        r <- r * r.mask
        raster::writeRaster(r, filename = t, overwrite = T)
    }
    cat("\n")
}

# Gets georeferenced GBIF data previously saved as a csv file (see geo-dataset-preparation::generateUncertaintyCSVFile)
getGBIFData <- function(){
    f.in <- paste0(data.processed.dir, "gbif/", date.string, "_coordinates_data.csv")
    gbif.df <- vroom(f.in, col_types = cols())
    return(gbif.df)
}

# Gets a sample of n.background points of environmental background data out of the predictor stack
getBackgroundData <- function(predictor.stack, n.background){
    set.seed(1234)
    bg <- sampleRandom(predictor.stack[[1]], n.background, xy=T) %>%
        as.data.frame()
    bg <- bg %>% dplyr::select(x, y) %>%
        bind_cols(data.frame(raster::extract(predictor.stack, bg[, c("x", "y")]))) %>%
        mutate(pa = 0)
    names(bg) <- c("lon", "lat", names(predictor.stack), "pa")
    return(bg)
}

# Runs a maxent model and saves model results
runMaxentModel <- function(model, force.run=T){
    output.dir <- paste0(model$output.dir, "maxunc_", model$max.uncertainty, "/", model$code)
    f.model.out <- paste0(output.dir, ".rds")
    if(file.exists(output.dir) & file.exists(f.model.out) & !force.run){
        cat(paste0("Model '", model$code, "' already exists. Skipping model generation ...\n"))
        return()
    }
    data.all <- vroom::vroom(model$data.file, col_types = cols())
    data.pa <- data.all %>% filter(pa == 1) %>% filter(unc <= model$max.uncertainty)
    data.bg <- data.all %>% filter(pa == 0)
    data <- bind_rows(data.pa, data.bg) %>% dplyr::select(-id, -unc)

    # cat("\n--------------------------------------------------------------------------------")
    cat("\nGENERATING MAXENT MODEL", model$code)
    cat("\n  PREDICTORS:", model$predictors)
    cat("\n  MAXENT ARGS:", model$maxent.args)
    cat(paste0("\n  DATA: ", table(data$pa)["1"], " presences, " , table(data$pa)["0"], " background points"), "\n")

    df.env <- data %>% dplyr::select(one_of(model$predictors)) %>% data.frame()

    set.seed(1234)
    mx <- dismo::maxent(df.env, p=data$pa, path=output.dir, args=model$maxent.args)
    saveRDS(mx, f.model.out)
    cat("  ==>", output.dir, "\n  ==>", paste0(output.dir, '.rds\n'))
    return(mx)
}

# Generate all n models for each modelling dataset prepared
fitMaxentModels <- function(model.taxon, force.run = F){
    dir.create(model.output.dir, recursive = T, showWarnings = F)
    predictor.stack <- getAllPredictorsStack(model.taxon)
    for(n in 1:model.taxon$models.nr){
        modelling.data.file <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_modelling_data_sample_", n, ".csv")
        model.code <- paste0("m_", n)
        model <- list("code" = model.code,
                      "name" = model.code,
                      "data.file" = modelling.data.file,
                      "predictors" = model.taxon$predictors,
                      "maxent.args" = model.taxon$maxent.args,
                      "output.dir" = model.output.dir,
                      "max.uncertainty" = model.taxon$max.uncertainty,
                      "force.run" = FALSE)
        f.out <- paste0(model.output.dir, model.code, ".rds")
        set.seed(1234)
        m <- runMaxentModel(model, force.run)
        if(!is.null(model.taxon$maxent.args.final)){
            model$code <- paste0(model$code, "_final")
            model$maxent.args <- model.taxon$maxent.args.final
            f.out <- paste0(model.output.dir, model.code, ".rds")
            m <- runMaxentModel(model, force.run)
        }
    }
}

# Prepares a background data frame ready for modelling
getBackgroundDF <- function(model.taxon){
    predictor.stack <- getAllPredictorsStack(model.taxon)
    f.out <- paste0(modelling.dir, "/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_background.csv")
    if(!file.exists(f.out)){
        set.seed(1234)
        bg <- getBackgroundData(predictor.stack, model.taxon$n.background)
        bg <- bg %>% mutate(species=model.taxon$taxon.name, unc=NA) %>%
            dplyr::select(pa, species, lon, lat, unc, all_of(names(predictor.stack)))
        vroom::vroom_write(bg, f.out)
        cat("    ==>", f.out)
    } else {
        bg <- vroom::vroom(f.out, col_types = cols())
    }
    return(bg)
}

# Adds predictor data to occurrences data frame
getOccurrencesEnvDF <- function(taxon.occ, model.taxon){
    predictor.stack <- getAllPredictorsStack(model.taxon)
    f.out <- paste0(modelling.dir, "/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_occ_env.csv")
    if(!file.exists(f.out)){
        v <- raster::extract(predictor.stack, taxon.occ[,c("lon", "lat")])
        taxon.env <- cbind(taxon.occ %>% mutate(pa=1) %>% dplyr::select(pa, species, lon, lat, unc), v)
        vroom::vroom_write(taxon.env, f.out)
        cat("    ==>", f.out, "\n")
    } else {
        taxon.env <- vroom::vroom(f.out, col_types = cols())
    }
    return(taxon.env)

}

# Prepares a data frame of occurrence coordinates and uncertainty
getOccurrences <- function(model.taxon, only.with.unc = F){
    f.out <- paste0(modelling.dir, "/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_gbif_data.csv")
    if(!file.exists(f.out)){
        if(!exists("gbif.df")){
            cat("\n  reading GBIF data file ...\n")
            gbif.df <- getGBIFData()
        }
        cat("  filtering taxon occurrences ...\n")
        taxon.df <- gbif.df %>% filter(str_starts(tolower(taxon), tolower(!!model.taxon$taxon.name)))
        taxon.region.df <- taxon.df %>% filter(lon >= model.taxon$region.ext@xmin & lon <= model.taxon$region.ext@xmax &
                                                   lat >= model.taxon$region.ext@ymin & lat <= model.taxon$region.ext@ymax)
        vroom::vroom_write(taxon.region.df, f.out)
        cat("    ==>", f.out)
    } else {
        cat("  reading file", f.out, "\n")
        taxon.region.df <- vroom::vroom(f.out, col_types = cols(), progress = F)
    }
    taxon.region.df <- taxon.region.df %>% dplyr::select(species=taxon, lon, lat, unc)
    if(only.with.unc){
        taxon.region.df <- taxon.region.df %>% filter(!is.na(unc))
    }
    return(taxon.region.df)
}

# Filter occurrences down to only one per grid cell. If more than one, select the one with least uncertainty
filterOccurrences <- function(taxon.occ, force.run = F){
    f.out <- paste0(modelling.dir, "/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_gbif_data_filtered_to_1_per_square_km.csv")
    if(file.exists(f.out) & !force.run){
        return(vroom(f.out, col_types = cols()))
    }
    coords <- taxon.occ[, c("lon", "lat")] %>% as.data.frame()
    predictor.stack <- getAllPredictorsStack(model.taxon)
    taxon.occ <- bind_cols(taxon.occ, cell = raster::cellFromXY(raster(predictor.stack), coords)) %>%
        distinct() %>%
        group_by(cell) %>%
        arrange(unc) %>%
        filter(row_number()==1) %>%
        ungroup()
    vroom::vroom_write(taxon.occ, f.out)
    cat("    ==>", f.out, "\n")
    return(taxon.occ)
}

# Gets the predictor set for the given taxon and region
getAllPredictorsStack <- function(model.taxon){
    d.predictors <- paste0(data.processed.dir, "predictors/", model.taxon$region.name, "/", model.taxon$resolution)
    f.predictors <- list.files(d.predictors, pattern="*.tif", full.names = T)
    predictor.stack <- stack(f.predictors)
    return(predictor.stack)
}
# Determines which variables to use according to a VIF less than 10 (Dormann et al. (2013)
# see https://onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2012.07348.x)
getPredictorsSelectedByVIF <- function(env.df, model.taxon, force.run = F){
    predictor.stack <- getAllPredictorsStack(model.taxon)
    f.out <- paste0(output.dir, "maxent/", model.taxon$taxon.abbrv, "-uncertainty/", model.taxon$taxon.abbrv, "_vif_analysis.rds")
    if(file.exists(f.out) & !force.run){
        vif.analysis <- readRDS(f.out)
    } else {
        set.seed(1234)
        vif.analysis <- vifstep(env.df %>% dplyr::select(names(predictor.stack)) %>% na.omit() %>% as.data.frame(), th=10)
        saveRDS(vif.analysis, f.out)
        cat("    ==>", f.out, "\n")
    }
    final.predictor.set <- vif.analysis@results$Variables
    env.df <- env.df %>% dplyr::select(pa, species, lon, lat, unc, all_of(final.predictor.set))
    return(env.df)
}

# Selects and lays out the columns properly
filterBackgroundDF <- function(bg, final.predictor.set){
    bg <- bg %>% dplyr::select(pa, species, lon, lat, all_of(final.predictor.set)) # Filter bg accordingly
    return(bg)
}

# Prepares a data set for all values of the grid cells contained withint the uncertainty boundary around each point
getPredictorValuesAroundUncertaintyBoundaries <- function(model.taxon, taxon.occ){
    f.out <- paste0(modelling.dir, "/data/",
                    model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_presences_env_data.csv")
    if(!file.exists(f.out)){
        nr.predictors <- dim(predictor.stack)[3]
        df.pa <- tibble()
        for(i in 1:nrow(taxon.occ)){
            occ <- taxon.occ[i,]
            coord <- cbind(occ$lon, occ$lat)
            values <- raster::extract(predictor.stack, coord, buffer=occ$unc)
            nr.values <- length(values[[1]])
            cat("  processing occurrence ", i, " at lon: ", occ$lon, " lat: ", occ$lat, ", uncertainty: ", occ$unc, ", number of values: ", nr.values/nr.predictors, "\n", sep = "")
            if(nr.values == 0)
                next()
            # When only values from one cell are returned they differ in str than when more than one is returned
            if(nr.values == nr.predictors){
                values = rbind(unlist(values))
            }
            row <- bind_cols(id=i, pa=1, occ, values)
            df.pa <- bind_rows(df.pa, row)
        }
        vroom::vroom_write(df.pa, f.out)
        cat("  ==>", f.out, "\n")
    } else {
        df.pa <- vroom::vroom(f.out, col_types=cols())
    }
    return(df.pa)
}

# Prepares modelling datasets, each one with only one sample value from each occurrence uncertainty area
# NOTE: This function accelerates the process but I did have problems for the first models corresponding to
#       the number of cores to use. With these, only the coordinates appeared in the resulting file. Needs
# checking
prepareModellingDatasets <- function(model.taxon, df.pa){
    # Sample files already done
    existing.ids <- list.files(paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/data/"), pattern = "*sample*") %>%
        str_split("\\_") %>% unlist() %>%
        str_replace("\\.csv", "")
    existing.ids <- sort(as.numeric(existing.ids[grepl("^-?[0-9.]+$", existing.ids)]))

    todo.ids <- 1:model.taxon$models.nr
    todo.ids <- todo.ids[!todo.ids %in% existing.ids]
    if(length(todo.ids) == 0)
        return()

    f.bg <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/data/",
                   model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_background.csv")
    bg <- read_delim(f.bg, delim = "\t", col_types = cols()) %>%
        mutate(id = row_number() + 10000) %>%
        relocate(id, .before = pa)

    cl <- parallel::makeCluster(nr.cores.to.use)
    doParallel::registerDoParallel(cl)
    datasets.results <- foreach(n=todo.ids, .combine = 'c', .packages='tidyverse') %dopar% {
        f.out <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_modelling_data_sample_", n, ".csv")
        df.mod <- df.pa %>%
            group_by(id) %>%
            do(sample_n(., 1)) %>%
            ungroup() %>%
            bind_rows(bg) %>%
            dplyr::select(id, pa, lon, lat, unc, all_of(model.taxon$predictors)) %>%
            filter_at(vars(model.taxon$predictors), all_vars(!is.na(.)))
        vroom::vroom_write(df.mod, f.out, progress = F)
    }
    parallel::stopCluster(cl)
}

# Generates predictions of each model for the study area
makeModelPredictionsParallel <- function(model.taxon, model.ids = 1:250){
    predictor.stack <- getAllPredictorsStack(model.taxon)[[model.taxon$predictors]]
    prediction.output.dir <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_", model.taxon$max.uncertainty, "/")
    model.output.dir <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/models/maxunc_", model.taxon$max.uncertainty, "/")
    dir.create(prediction.output.dir, showWarnings = F, recursive = T)

    cl <- parallel::makeCluster(nr.cores.to.use)
    doParallel::registerDoParallel(cl)
    prediction.results <- foreach(n=model.ids, .combine = 'c', .packages=c("dismo", "terra")) %dopar% {
        f.out <- paste0(prediction.output.dir, "m_", n, "_final_prediction.tif")
        if(!file.exists(f.out)){
            m <- readRDS(paste0(model.output.dir, "m_", n, "_final.rds"))
            prediction <- dismo::predict(m, predictor.stack)
            raster::writeRaster(prediction, f.out, overwrite = T, options = c("COMPRESS = LZW"))
        }
    }
    parallel::stopCluster(cl)
}

makeModelPredictions <- function(model.taxon){
    predictor.stack <- getAllPredictorsStack(model.taxon)[[model.taxon$predictors]]
    model.output.dir <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/predictions/")
    dir.create(model.output.dir, showWarnings = F, recursive = T)
    f.out <- paste0(model.output.dir, "m_", n, "_prediction.tif")
    cat("Predicting for model", n, "...\n")
    for(i in 1:model.taxon$models.nr){
        if(!file.exists(f.out)){
            m <- readRDS(paste0(model.output.dir, "../models/m_", n, ".rds"))
            prediction <- dismo::predict(m, predictor.stack)
            terra::writeRaster(prediction, f.out, overwrite = T)
        }
    }
}

# Collects maxent results from maxentResults.csv from all fitted models.
collectMaxentModelResults <- function(model.taxon, force.run = T){
    f.out <- paste0(modelling.dir, "/", model.taxon$taxon.abbrv, "_maxunc_", model.taxon$max.uncertainty, "_maxent_results.csv")
    if(!file.exists(f.out) | force.run){
        df.maxent.results <- tibble()
        for(n in 1:model.taxon$models.nr){
            f.in <- paste0(model.output.dir, "maxunc_", model.taxon$max.uncertainty, "/m_", n, "/maxentResults.csv")
            if(file.exists(f.in)){
                df <- vroom::vroom(f.in, col_types = cols())
                df.maxent.results <- rbind(df.maxent.results, cbind(model.id=n, df))
            }
        }
        vroom::vroom_write(df.maxent.results, f.out)
        cat("    ==>", f.out, "\n")
    } else {
        cat("    Reading", f.out, "...\n")
        df.maxent.results <- vroom::vroom(f.out, col_types = cols())
    }
    return(df.maxent.results)
}

calculateRanges <- function(model.taxon){
    cl <- parallel::makeCluster(nr.cores.to.use)
    doParallel::registerDoParallel(cl)
    f.out <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/", model.taxon$taxon.abbrv,
                    "_maxunc_", model.taxon$max.uncertainty, "_suitability_range_results.csv")
    model.output.dir <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/models/")
    prediction.output.dir <- paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/predictions/")
    df.range.results <- foreach(n=1:model.taxon$models.nr, .combine=rbind, .packages = c("vroom", "terra", "magrittr")) %dopar% {
        prediction <- terra::rast(paste0(prediction.output.dir, "maxunc_", model.taxon$max.uncertainty, "/m_", n, "_final_prediction.tif"))
        suitability <- sum(na.omit(terra::values(prediction)))
        df <- vroom::vroom(paste0(model.output.dir, "maxunc_", model.taxon$max.uncertainty, "/m_", n, "_final/maxentResults.csv"), col_types = vroom::cols())
        threshold <- df$`Maximum training sensitivity plus specificity Cloglog threshold`
        range.area <- sum(na.omit(terra::values(prediction >= threshold)))
        cbind(n, maxSSS=threshold, suitability=suitability, range.area=range.area)
    }
    parallel::stopCluster(cl)
    vroom::vroom_write(df.range.results %>% as_tibble(), f.out)
    cat("    ==>", f.out, "\n")
    cat(paste0(rep("-", 80), collapse = ""), "\n")
}

# Prepare predictors
prepareTaxonPredictors <- function(){
    cat("Preparing predictors ...\n")
    cat("  bioclim predictors ...\n")
    for(i in 1:length(taxons.to.model)){
        model.taxon <- taxons.to.model[[i]]
        predictor.stack <- prepareBioclimPredictors(model.taxon)
        # Apply mask if needed
        if(model.taxon$region.name %in% c("australia", "northern_canada_and_greenland")){
            cat("  applying mask ...\n")
            applyMask(model.taxon)
        }
    }
    cat("  tree cover predictor ...\n")
    source("scripts/R/tree_cover.R")
}

#### MODELLING -------------------------------------------------------------------------------- ####
# Setting up modelling environment
pnr.total.cores <- detectCores()
nr.cores.to.use <- pnr.total.cores - 4
cat(paste0(rep("=", 80), collapse = ""), "\n")

for(i in c(12)){#length(taxons.to.model)){
    model.taxon <- taxons.to.model[[i]]
    cat("Modelling ", model.taxon$taxon.name, "\n")
    cat("  setting up script parameters and folders ..\n")
    modelling.dir <- paste0(output.dir, "maxent/", model.taxon$taxon.abbrv, "-uncertainty")
    model.output.dir <- paste0(modelling.dir, "/models/")
    dir.create(model.output.dir, recursive=T, showWarnings = FALSE)
    model.data.dir <- paste0(modelling.dir, "/data")
    dir.create(model.data.dir, recursive=T, showWarnings = FALSE)
    predictor.stack <- getAllPredictorsStack(model.taxon)

    cat("  getting background data ...\n")
    bg <- getBackgroundDF(model.taxon)

    cat("  getting occurrence from GBIF data ...\n")
    taxon.occ <- getOccurrences(model.taxon, T)

    cat("  filtering occurrences to one per grid cell ...\n")
    taxon.occ <- filterOccurrences(taxon.occ, F)

    cat("  getting environmental data for the occurrences ...\n")
    taxon.env <- getOccurrencesEnvDF(taxon.occ, model.taxon)

    cat("  building occurrence + background dataset ...\n")
    occ.bg.df <- rbind(taxon.env, bg)

    cat("  subsetting data to final set of predictors as determined by VIF analysis ...\n")
    occ.bg.df <- getPredictorsSelectedByVIF(occ.bg.df, model.taxon, F)

    cat("  adapting predictor stack to selected predictors ...\n")
    final.predictor.set <- names(occ.bg.df)[!names(occ.bg.df) %in% c("pa", "species", "lon", "lat", "unc")]
    predictor.stack <- getAllPredictorsStack(model.taxon)
    predictor.stack <- predictor.stack[[final.predictor.set]]
    model.taxon$predictors <- final.predictor.set

    cat("  gathering predictor values at occurrences ...\n")
    df.pa <- getPredictorValuesAroundUncertaintyBoundaries(model.taxon, taxon.occ)

    cat("  generating modelling datasets ...\n")
    prepareModellingDatasets(model.taxon, df.pa)

    cat("  generating models ...\n")
    fitMaxentModels(model.taxon, T)

    cat("  predicting suitability ...\n")
    makeModelPredictionsParallel(model.taxon, 1:250)
    makeModelPredictionsParallel(model.taxon, 251:500)

    cat("  gathering maxent results ...\n")
    cmd <- paste0("mv ../NHC-GeoUncertainty_storage/outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/models/maxunc_", model.taxon$max.uncertainty,
                  " outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/models")
    system(cmd)
    cmd <- paste0("mv ../NHC-GeoUncertainty_storage/outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_", model.taxon$max.uncertainty,
                  " outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/predictions/")
    system(cmd)
    df.maxent.results <- collectMaxentModelResults(model.taxon, T)

    cat("  calculating suitability results ...\n")
    calculateRanges(model.taxon)
    cmd <- paste0("mv outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/models/maxunc_", model.taxon$max.uncertainty,
                  " ../NHC-GeoUncertainty_storage/outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/models/")
    system(cmd)
    cmd <- paste0("mv outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_", model.taxon$max.uncertainty,
                  " ../NHC-GeoUncertainty_storage/outputs/maxent/", model.taxon$taxon.abbrv, "-uncertainty/predictions/")
    system(cmd)
}
