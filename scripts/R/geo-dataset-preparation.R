# This script prepares a dataframe from the geo sqlite database and stores it as an sqlite database.
# The process is done by querying the sqlite database by chunks of records in order
# not to run into memory problems (1 million have worked for us but modify accordingly to your needs)
# Each chunk is processed and written to a new processed sqlite database.

source("scripts/R/init.R")

processGeoDF <- function(df){
    # 'null' --> ''
    cols <- c("gbif.id", "occ.id", "lon", "lat", "unc", "vb.coordsystem", "vb.srs",
              "georef.protocol", "geo.issue", "issue", "event.date")
    for(col in cols){
        df <- df %>% mutate(!!col := ifelse(!!as.name(col) == 'null', '', !!as.name(col)))
    }
    # 'null' --> 'UNKNOWN'
    cols <- c("ds.name", "inst.code", "pub.country", "pub.key", "continent", "country.code",
              "kingdom", "phylum", "class", "order", "family", "sci.name")
    for(col in cols){
        df <- df %>% mutate(!!col := ifelse(!!as.name(col) == 'null', 'UNKNOWN', !!as.name(col)))
    }
    # other
    df <- df %>% mutate(country.code = ifelse(country.code == 'ZZ' | country.code == 'null', 'UNKNOWN', country.code)) %>%
        mutate(country.code = ifelse(country.code == 'XZ', 'INTL. WATERS', country.code)) %>%
        mutate(event.year = ifelse(event.year == '', '', str_sub(event.year, 1, 4)))
    return(df)
}

processGeoUncDF <- function(df){
    df <- df %>% mutate(unc=as.numeric(unc)) %>%
        mutate(unc.cat = cut(unc, breaks = c(0, 1, 10, 100, 250, 1000, 5000, 10000, 50000, 100000, 25000000),
                             include.lowest = T)) %>%
        mutate(unc.cat = fct_recode(unc.cat,
                                    "<= 1m"="[0,1]", "10m"="(1,10]", "100m"="(10,100]",
                                    "250m"="(100,250]", "1km"="(250,1e+03]", "5km"="(1e+03,5e+03]",
                                    "10km"="(5e+03,1e+04]", "50km"="(1e+04,5e+04]", "100km"="(5e+04,1e+05]",
                                    "> 100km"="(1e+05,2.5e+07]"))
    return(df)
}

processGeoColTypesDF <- function(df){
    df <- df %>% mutate(lon = as.double(lon), lat = as.double(lat), unc = as.double(unc),
                        event.date = as.Date(event.date, "%Y-%m-%d"), event.year = as.integer(event.year),
                        unc.cat = as.character(unc.cat))
    return(df)
}

generateGeoreferencedDatabase <- function(){
    cat(rep("=", 80), "\n", "Generating geo records database ...\n", sep="")
    geo.out.cols <- c("gbif.id", "occ.id", "ds.name", "inst.code", "pub.country", "pub.key",
                      "lon", "lat", "unc", "vb.coordsystem", "vb.srs",
                      "georef.protocol", "geo.issue", "issue",
                      "event.year",
                      "continent", "country.code",
                      "kingdom", "phylum", "class", "order", "family", "genus", "sci.name",
                      "unc.cat")

    sql <- paste("select", paste0(gbif.cols$gbif.col, collapse=", "),
                 "from ps where hascoordinate = 'true'")
    con <- dbConnect(RSQLite::SQLite(), f.db.whole)
    con.w <- dbConnect(RSQLite::SQLite(), f.db.whole.geo.processed)
    nrows <- 1000000
    offset <- 0
    n <- 0
    repeat{
        sql.rows <- paste0(sql, " limit ", nrows, " offset ", offset)
        df <- dbGetQuery(con, sql.rows) %>% setNames(., gbif.cols$var)
        if(nrow(df) == 0){
            break
        }
        df <- processGeoDF(df)
        df <- processGeoUncDF(df)
        df <- processGeoColTypesDF(df)
        # We consider records at [lon=0, lat=0] as errors
        df <- df %>% filter(lon != 0 & lat !=0)
        dbWriteTable(con.w, "ps", df %>% dplyr::select(all_of(geo.out.cols)) %>%
                         setNames(., gsub("\\.", "\\_", geo.out.cols)), append=T)
        offset <- offset + nrows
        n <- n + nrow(df)
        cat(" ==>", n, "records inserted so far ...\n")
    }
    dbDisconnect(con)
    cat("Creating index for kingdom column\n")
    dbSendQuery(con.w, "create index idx_kingdom on ps(kingdom);")
    dbDisconnect(con.w)
    cat(rep("-", 80), "\n", sep="")
}

generateUncertaintyCSVFile <- function(only.with.unc = T){
    cat(rep("=", 80), "\n", "Generating uncertainty csv file ...\n", sep="")
    if(only.with.unc){
        f.out <- paste0(data.processed.dir, "gbif/", date.string, "_uncertainty_data.csv")
    } else {
        f.out <- paste0(data.processed.dir, "gbif/", date.string, "_coordinates_data.csv")
    }
    unlink(f.out)
    tibble(occ.id = character(), pub.country = character(), pub.key = character(), dataset = character(),
           lon = double(), lat = double(),
           unc = double(), unc.cat = character(),
           kingdom = character(), phylum = character(), class = character(), order = character(),
           family = character(), genus = character(), taxon = character(), year = integer()) %>%
        vroom_write(f.out, delim=",", col_names = T)

    con <- dbConnect(RSQLite::SQLite(), f.db.whole.geo.processed)
    nrows <- 1000000
    offset <- 0
    n.processed <- 0
    n.inserted <- 0
    repeat{
        sql.rows <- paste0("select * from ps limit ", nrows, " offset ", offset)
        df.geo <- dbGetQuery(con, sql.rows)
        if(nrow(df.geo) == 0){
            break
        }
        df <- df.geo

        if(only.with.unc)
            df <- df %>% filter(!is.na(unc))

        df <- df %>%
            dplyr::select(occ.id=occ_id, pub.country=pub_country, pub.key=pub_key, dataset=ds_name,
                          lon, lat, unc, unc.cat=unc_cat,
                          kingdom, phylum, class, order, family, genus, taxon=sci_name, year=event_year) %>%
            vroom_write(f.out, delim=",", append = T)
        offset <- offset + nrows
        n.processed <- n.processed + nrow(df.geo)
        n.inserted <- n.inserted + nrow(df)
        cat(" ==>", n.processed, "records processed,", nrow(df), "inserted -->", n.inserted, "\n")
    }
    dbDisconnect(con)
    cat(rep("-", 80), "\n", sep="")
}

generateUncertaintyByDecimalDegreeCSVFile <- function(){
    cat(rep("=", 80), "\n", "Generating uncertainty by decimal degree csv file ...\n", sep="")
    # The following peaks have been previously identified by analysing the distribution of unc values.
    peaks <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 71, 100, 150, 200, 250, 300, 402, 500, 707, 1000,
           1415, 2000, 2500, 3036, 3535, 5000, 7071, 10000, 14143, 20000, 25000, 30000, 50000, 100000, 999999)

    f <- paste0(data.processed.dir, "gbif/", date.string, "_uncertainty_data.csv")
    df <- vroom(f, col_types = cols())
    dfa <- df %>% mutate(unc = factor(round(unc), levels = as.character(peaks)))
    df.res <- tibble()
    for(peak in peaks){
        cat("   processing uncertainty peak", peak, '...\n')
        dfb  <- dfa %>%
            filter(unc == !!peak) %>%
            mutate(lon = floor(lon) + 0.5, lat = floor(lat) + 0.5) %>%
            group_by(lon, lat) %>%
            count() %>%
            mutate(peak = !!peak)
        df.res <- bind_rows(df.res, dfb)
    }
    f.out <- paste0(data.processed.dir, "gbif/", date.string, "_uncertainty_by_decimal_degree.csv")
    vroom_write(df.res, f.out, delim = ',')
    cat("==>", f.out, "\n")
    cat(rep("-", 80), "\n", sep="")
}

generateGeoreferencedDatabase()
generateUncertaintyCSVFile()
generateUncertaintyByDecimalDegreeCSVFile()
