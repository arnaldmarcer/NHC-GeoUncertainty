source("scripts/R/init.R")

# Generate list of extents that define each tile to process
getTileExtents <- function(size, unc){
    prediction.dir <- paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_", unc, "/")
    r.template <- rast(paste0(prediction.dir, "m_1_final_prediction.tif"))
    sa.ext <- ext(r.template)
    xs <- seq(sa.ext[1], sa.ext[2] + size, size)
    ys <- seq(sa.ext[3], sa.ext[4] + size, size)
    tiles <- list()
    i <- 1
    while(i < length(xs)){
        j <- 1
        while(j < length(ys)){
            ext.terra = ext(xs[i], xs[i+1], ys[j], ys[j+1])
            tiles <- c(tiles, ext.terra)
            j <- j + 1
        }
        i <- i + 1
    }
    return(tiles)
}

# Mosaics all calculated sd tiles into an sd raster for the whole area
mosaicTilesGDAL <- function(taxon, aggregate="sd", unc = Inf){
    # Generate virtual mosaicked layer
    vrt.file <- paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/", aggregate, "/mosaic.vrt")
    cat("  mosaicking tiles ...\n")
    cmd <- paste("gdalbuildvrt", vrt.file, paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/", aggregate, "/*.tif"))
    system(cmd)
    # Translate into final tif layer
    mosaicked.file <- paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_", unc, "/", taxon$taxon.abbrv, "_prediction_", aggregate, ".tif")
    cmd <- paste("gdal_translate -of GTiff -co 'COMPRESS=LZW'", vrt.file, mosaicked.file)
    system(cmd)
    cat("  ==>", mosaicked.file, "\n")
    unlink(vrt.file)
    return(rast(mosaicked.file))
}

# Calculates sd for a given tile, the function knows where to find files
# Warning: End of function deletes tif files for the tile from which the sd have been calculated
generateTileSD <- function(taxon, t){
    # Calculate standard deviation
    cat("Calculting sd for tile", t, "\n")
    f.std.out <- paste0(paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/sd/tile_", t, ".tif"))
    if(file.exists(f.std.out)){
        cat("  sd tiles", t, "already exist, skipping ...\n")
        return()
    }
    tifs <- list.files(paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/sd/"), pattern = paste0("*_", t, ".tif"), full.names = T)
    if(length(tifs) == 1){
        cat("  tile", t, "is NA, skipping from calculating sd, writing tile as is ...\n")
        file.copy(from=tifs, to=f.std.out)
        cat("    ==>", f.std.out, "...\n")
        cat(" removing files ...\n")
        do.call(file.remove, list(tifs))
        return()
    }
    st.tifs <- rast(tifs)
    if(!file.exists(f.std.out)){
        cat(" calculating tile standard deviation ...\n")
        tile.std.r <- app(st.tifs, fun=sd)
        writeRaster(tile.std.r, filename = f.std.out, gdal=c("COMPRESS=LZW"), overwrite = T)
        cat("    ==>", f.std.out, "...\n")
    }
    cat(" removing files ...\n")
    do.call(file.remove, list(tifs))
}

# Break tifs into tiles, calculate std for each tile, and write raster to disk
# There are 500 tifs, for each tile extent we have 500 tiles corresponding to the 500 whole area tifs
# Then, with these we can calculate the standard deviation for each tile, corresponding to its 500 models
# Finally, we mosaic all standard deviation tiles into one single piece again, the whole area
calculatePredictionsSD <- function(taxon, tiles, unc){
    dir.create(paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/sd"), recursive = T, showWarnings = F)
    prediction.dir <- paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_", unc)
    prediction.tifs <- list.files(prediction.dir, pattern=".tif", full.names = T)
    for(t in 1:length(tiles)){
        # Check if tile already exists
        f.std.out <- paste0(paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/sd/tile_", t, ".tif"))
        if(file.exists(f.std.out)){
            cat("  sd tile", t, "already exists, skipping ...\n")
            next
        }
        tile <- tiles[[t]]
        cat("  processing tile", t, paste(tile), "\n")
        for(i in 1:length(prediction.tifs)){
            tif <- prediction.tifs[i]
            cat("    processing tif file", tif, "\n")
            model <- rev(unlist(str_split(tif, "\\_")))[3]
            f.out <- paste0(tmp.dir, "tiles/", taxon$taxon.abbrv, "/sd/model_", model, "_tile_", t, ".tif")
            # Check if tile is all NA, if it is we can skip this tif and just write to disk an NA piece
            cat("  checking if tile is all NA ...\n")
            v.tile <- values(terra::crop(rast(tif), tile))
            if(length(v.tile[!is.na(v.tile)]) == 0){
                cat("    writing NA tile, need only to write once, can skip rest of tifs ...\n")
                cmd <- paste("gdalwarp -t_srs EPSG:4326 -te", tile[1], tile[3], tile[2], tile[4], "-overwrite", tif, f.out)
                system(cmd)
                break
            }
            if(file.exists(f.out)){
                break
            }
            cat("  writing non-NA tile ...\n")
            cmd <- paste("gdalwarp -t_srs EPSG:4326 -te", tile[1], tile[3], tile[2], tile[4], "-overwrite", tif, f.out)
            system(cmd)
        }
        # Calculate standard deviation
        generateTileSD(taxon, t)
    }
}

# Iterate over each taxon: calculate tiles, process tiles, generate sd layer
size <- 10 # This size works well in my computer
for(i in 2:length(taxons.to.model)){
    taxon <- taxons.to.model[[i]]
    cat(paste0(rep("=", 80)), "\n", sep="")
    cat("Calculating standard deviation for predictions of", taxon$taxon.name, "...\n")
    tiles <- getTileExtents(size, Inf)
    calculatePredictionsSD(taxon, tiles, Inf)
    r <- mosaicTilesGDAL(taxon, "sd", Inf)
}
# Delete all temp files
unlink(paste0(tmp.dir, "tiles/"), recursive = T)

