library(RSQLite)
library(vroom)
library(rnaturalearth)
library(ggridges)
library(sf)
library(sp)
library(GADMTools)
library(tidyverse)
library(viridis)
library(ggplot2)
library(cowplot)
library(dismo)
library(terra)
library(luna)
library(doParallel)
library(rvest)
library(stringr)
library(httr)
library(gtools)
library(usdm)
library(rJava)
library(tictoc)
library(scales)
library(rlist)
library(vroom)

#### DIRECTORIES AND FILES ==================================================================== ####
date.string <- "20210311T1049"
output.dir <- "outputs/"
data.raw.dir <- "data/raw/"
data.processed.dir <- "data/processed/"
manuscript.dir <- "manuscript/"
tmp.dir <- "tmp/"
dataset <- "whole"

# Databases
f.db.whole <- paste0(data.raw.dir, "gbif/", date.string, "_preserved_specimen_sqlite.db")
f.db.whole.geo.processed <- paste0(data.processed.dir, "gbif/", date.string, "_processed_whole_dataset_georeferenced.sqlite")
f.db.whole.non.geo.processed <- paste0(data.processed.dir, "gbif/", date.string, "_processed_whole_dataset_non_georeferenced.sqlite")

# Results files
f.out.non.geo.whole.results <- paste0(output.dir, date.string, "_non_geo_whole_results.rds")
f.out.geo.whole.results <- paste0(output.dir, date.string, "_geo_whole_results.rds")
f.results <- paste0(manuscript.dir, "results.dat")

# Variables
gbif.cols <- data.frame(rbind(
    cbind(var="gbif.id", gbif.col="gbifID"),
    cbind("occ.id", "occurrenceID"),
    cbind("ds.name", "datasetName"),
    cbind("inst.code", "institutionCode"),
    cbind("pub.country", "publishingcountry"),
    cbind("pub.key", "publishingorgkey"),
    cbind("lon", "decimalLongitude"),
    cbind("lat", "decimalLatitude"),
    cbind("unc", "coordinateUncertaintyInMeters"),
    cbind("vb.coordsystem", "verbatimCoordinateSystem"),
    cbind("vb.srs", "verbatimSRS"),
    cbind("georef.date", "georeferencedDate"),
    cbind("georef.protocol", "georeferenceProtocol"),
    cbind("geo.issue", "hasGeospatialIssues"),
    cbind("issue", "issue"),
    cbind("event.date", "eventdate_readable"),
    cbind("continent", "continent"),
    cbind("country.code", "countryCode"),
    cbind("kingdom", "kingdom"),
    cbind("phylum", "phylum"),
    cbind("class", "class"),
    cbind("order", "order_"),
    cbind("family", "family"),
    cbind("genus", "genus"),
    cbind("sci.name", "acceptedScientificName"),
    cbind("hascoordinate", "hascoordinate"),
    cbind("event.year", "eventdate_readable")))

#### TAXONS TO MODEL ========================================================================== ####
# Model with all occurrences
taxons.to.model <- list(
    # RHODODENDRON GROENLANDICUM ----------------------------------------------------------------- #
    # 1- Model with all occurrences
    rhododendron_groenlandicum_Inf = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = Inf,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 2 - Model with all occurrences below 3536m (5km resolution)
    rhododendron_groenlandicum_3536 = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 3536,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 3 - Model with all occurrences below 7071m (10km resolution)
    rhododendron_groenlandicum_7071 = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 7071,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 4 - Model with all occurrences below 14142m (20km resolution)
    rhododendron_groenlandicum_14142 = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 14142,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 5 - Model with all occurrences below 28284 (40km resolution)
    rhododendron_groenlandicum_28284 = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 28284,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 6 - Model with all occurrences below 70711m (100km resolution)
    rhododendron_groenlandicum_70711 = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 70711,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 7 - Model with all occurrences below 141421m (200km resolution)
    rhododendron_groenlandicum_141421 = list(
        taxon.abbrv = "rhododendron_groenlandicum",
        taxon.name = "Rhododendron Groenlandicum", # Rhododendron Groenlandicum (Oeder) Kron & Judd
        region.name = "northern_canada_and_greenland",
        region.ext = extent(-168, -10, 43, 84),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 141421,
        predictors = paste0("bio", 1:19),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # GUAZUMA ULMIFOLIA -------------------------------------------------------------------------- #
    # 8 - Model with all occurrences
    guazuma_ulmifolia_Inf = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = Inf,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 9 - Model with all occurrences below 3536m (5km resolution)
    guazuma_ulmifolia_maxunc_3536 = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 3536,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 10 - Model with all occurrences below 7071m (10km resolution)
    guazuma_ulmifolia_maxunc_7071 = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 7071,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 11 - Model with all occurrences below 14142m (20km resolution)
    guazuma_ulmifolia_maxunc_14142 = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 14142,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 12 - Model with all occurrences below 28284m (40km resolution)
    guazuma_ulmifolia_maxunc_28284 = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 28284,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 13 - Model with all occurrences below 70711m (10km resolution)
    guazuma_ulmifolia_maxunc_70711 = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 70711,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 14 - Model with all occurrences below 141421m (100km resolution)
    guazuma_ulmifolia_maxunc_141421 = list(
        taxon.abbrv = "guazuma_ulmifolia",
        taxon.name = "Guazuma ulmifolia Lam.",
        region.name = "mexico_to_argentina",
        region.ext = extent(-118, -34, -38, 32),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 141421,
        predictors = c(paste0("bio", 1:19), "treecover"),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # EUCALYPTUS GONGYLOCARPA -------------------------------------------------------------------- #
    # 15 - Model with all occurrences
    eucalyptus_gongylocarpa_Inf = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = Inf,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 16 - Model with all occurrences below 3536m (5km resolution)
    eucalyptus_gongylocarpa_3536 = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 3536,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 17 - Model with all occurrences below 7071m (10km resolution)
    eucalyptus_gongylocarpa_7071 = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 7071,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 18 - Model with all occurrences below 14142m (20km resolution)
    eucalyptus_gongylocarpa_14142 = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 14142,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 19 - Model with all occurrences below 28284m (40km resolution)
    eucalyptus_gongylocarpa_28284 = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 28284,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 20 - Model with all occurrences below 70711m (100km resolution)
    eucalyptus_gongylocarpa_70711 = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 70711,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/")),
    # 21 - Model with all occurrences below 141421m (200km resolution)
    eucalyptus_gongylocarpa_141421 = list(
        taxon.abbrv = "eucalyptus_gongylocarpa",
        taxon.name = "Eucalyptus gongylocarpa Blakely",
        region.name = "australia",
        region.ext = extent(113, 154, -39.5, -10),
        models.nr = 500,
        n.background = 10000,
        max.uncertainty = 141421,
        predictors = c(paste0("bio", 1:19)),
        resolution = "30s",
        maxent.args = c("autofeature=TRUE", "randomtestpoints=20"),
        maxent.args.final = c("autofeature=TRUE"),
        modelling.dir = paste0(output.dir, "maxent/"))
    )


#### WORLD MAP ================================================================================ ####
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world %>%
    # Correct code for Kosovo
    mutate(iso_a2 = ifelse(name == 'Kosovo', 'XK', iso_a2)) %>%
    # Mark Somaliland as SO to match GBIF data
    mutate(iso_a2 = ifelse(name == 'Somaliland', 'SO', iso_a2))
world_coastline <- ne_coastline(scale = "medium", returnclass = "sf")

# This is to be used to repair codes ('colony.code') in GBIF which are not present in rnaturalearth countries.
# Reference: https://en.wikipedia.org/wiki/List_of_ISO_3166_country_codes
colonies <- tibble(
    colony.code = c("BQ", "BV", "CC", "CX", "GF", "GI", "GP", "MQ", "RE", "YT", "SJ", "TK", "TV", "UM", "XK"),
    country.code =     c("NL", "NW", "AU", "AU", "FR", "GB", "FR", "FR", "FR", "FR", "NW", "NZ", "TV", "US", "XK"),
    country.name = c("Bonaire, Sint Eustatius and Saba", "Bouvet Island", "The Territory of Cocos (Keeling) Islands",
                     "The Territory of Christmas Island", "Guyane", "Gibraltar", "Guadeloupe", "Martinique", "Réunion",
                     "The Department of Mayotte", "Svalbard and Jan Mayen", "Tokelau", "Tuvalu",
                     "Baker Island, Howland Island, Jarvis Island, Johnston Atoll, Kingman Reef, Midway Atoll, Navassa Island, Palmyra Atoll, and Wake Island",
                     "Kosovo"))

#### FUNCTIONS ================================================================================ ####
# Writes results in the form of key:value
writeResult <- function(key.string, new.value){
    key.string <- str_trim(key.string)
    if(!file.exists(f.results)){
        cat(key.string, ":", new.value, "\n", file = f.results, sep = "")
    } else {
        df <- read.csv2(f.results, sep=":", header=F, stringsAsFactors = F) %>%
            setNames(., c("key", "value"))

        # new.value <- paste(" ", str_trim(new.value))
        if(key.string %in% df$key){
            df <- df %>% mutate(value=ifelse(key==key.string, new.value, value))
        } else {
            df <- rbind(df, data.frame(cbind(key=key.string, value=new.value)))
        }
        write_delim(df, file=f.results, delim = ":", col_names = F)
    }
    # cat("'", key.string, ": ", new.value, "' --> ", f.results, "\n", sep="")
}

# Reproject coordinates
coordTrans <- function(coords, from.epsg, to.epsg){
    coordinates(coords) <- c(1,2)
    proj4string(coords) <- CRS(paste("+init=epsg:", from.epsg, sep=""))
    df.trans <- data.frame(spTransform(coords, CRS(paste("+init=epsg:", to.epsg, sep=""))))
    names(df.trans) <- c("x", "y")
    return(df.trans)
}

