# This script prepares a dataframe from the non-geo sqlite database and stores it as an sqlite database.
# The process is done by querying the sqlite database by chunks of records in order not to run into
# memory problems (1 million have worked for us but modify accordingly to your needs).
# Each chunk is processed and written to a new processed sqlite database.

source("scripts/R/init.R")
dir.create(paste0(data.processed.dir, "gbif"), recursive = T, showWarnings = F)
processNonGeoDF <- function(df){
    # 'null' --> ''
    cols <- c("gbif.id", "occ.id", "issue", "event.date", "event.year")
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

geo.out.cols <- c("gbif.id", "occ.id", "ds.name", "inst.code", "pub.country", "pub.key", "issue",
                  "continent", "country.code",
                  "kingdom", "phylum", "class", "order", "family", "genus", "sci.name",
                  "event.date", "event.year")


sql <- paste0("select gbifid, occurrenceid, datasetname, institutioncode, publishingcountry, publishingorgkey, ",
              "issue, continent, countrycode, ",
              "kingdom, phylum, class, order_, family, genus, acceptedscientificname, eventdate_readable, ",
              "eventdate_readable from ps where hascoordinate = 'false'")

cat(rep("=", 80), "\n", "Generating non-geo records database ...\n", sep="")
con <- dbConnect(RSQLite::SQLite(), f.db.whole)
con.w <- dbConnect(RSQLite::SQLite(), f.db.whole.non.geo.processed)
nrows <- 1000000
offset <- 0
n <- 0
repeat{
    sql.rows <- paste0(sql, " limit ", nrows, " offset ", offset)
    df <- dbGetQuery(con, sql.rows) %>% setNames(., all_of(geo.out.cols))
    if(nrow(df) == 0){
        break
    }
    df <- processNonGeoDF(df)
    dbWriteTable(con.w, "ps", df %>% dplyr::select(geo.out.cols) %>%
                     setNames(., gsub("\\.", "\\_", all_of(geo.out.cols))), append=T)
    offset <- offset + nrows
    n <- n + nrow(df)
    cat(" ==>", n, "records inserted ...\n")
}
dbDisconnect(con)
dbDisconnect(con.w)
cat(rep("-", 80), "\n", sep="")
