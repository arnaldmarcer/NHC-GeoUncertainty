# This script reads the geo database and computes partial results so that scripts run faster down the line

source("scripts/R/init.R")
results <- list()

# Initialize number of records to 0
results[['world.numbers']]$georeferenced.all <- 0 # With coordinates
results[['world.numbers']]$georeferenced.unc <- 0 # With coordinates and uncertainty

# We compute partial results in record chunks of 'nrows' to make it more palatable for the server
nrows <- 1000000
cat(rep("=", 80), "\n", sep = "")
cat("Generating results for geo dataset ...\n")
cat("  getting partial results in", formatC(nrows, big.mark = ' ', digits = 0, format = 'f'), "record chunks\n")
con <- dbConnect(RSQLite::SQLite(), f.db.whole.geo.processed)

offset <- 0
n <- 0
repeat{
    sql.rows <- paste0("select * from ps", " limit ", nrows, " offset ", offset)
    df.geo <- dbGetQuery(con, sql.rows)
    names(df.geo) <- gsub("\\_", ".", names(df.geo))

    if(nrow(df.geo) == 0) # run out of records to process
        break
    # SUMMARY NUMBERS ---------------------------------------------------------------------------- #
    results[['summary.numbers']]$summary.taxonomy.country.dataset.geo <- bind_rows(
        results[['summary.numbers']]$summary.taxonomy.country.dataset.geo,
        df.geo %>%
            group_by(kingdom, phylum, class, country.code, ds.name, pub.country, pub.key) %>%
            summarise(n = n(), .groups = 'drop'))

    # EVENT DATE NUMBERS ------------------------------------------------------------------------- #
    results[['event.numbers']]$event.year.georeferenced.all <- bind_rows(
        results[['event.numbers']]$event.year.georeferenced.all,
        df.geo %>%
            mutate(event.year = as.numeric(event.year)) %>%
            mutate(has.unc = ifelse(is.na(unc), 'N', 'Y')) %>%
            group_by(country.code, kingdom, phylum, class, event.year, has.unc) %>%
            count() %>%
            pivot_wider(names_from = has.unc, values_from = n) %>%
            rename("n.geo"="N", "n.geo.unc"="Y"))

    results[['event.numbers']]$event.year.georeferenced.unc <- bind_rows(
        results[['event.numbers']]$event.year.georeferenced.unc,
        df.geo %>%
            filter(unc != '') %>%
            mutate(unc = as.double(unc)) %>%
            group_by(country.code, kingdom, phylum, class, event.year) %>%
            summarise(mean.unc=mean(unc), sd.unc=sd(unc), n=n(), .groups = 'drop'))

    # WORLD-WIDE NUMBERS OF RECORDS WITH COORDINATES AND RECORDS WITH COORDINATES AND GEOREFERENCED
    results[['world.numbers']]$georeferenced.all <-
        results[['world.numbers']]$georeferenced.all + nrow(df.geo)
    results[['world.numbers']]$georeferenced.unc <-
        results[['world.numbers']]$georeferenced.unc + df.geo %>% filter(unc != '') %>% nrow()

    results[['world.numbers']]$georeferenced.unc.cat <- bind_rows(
        results[['world.numbers']]$georeferenced.unc.cat,
        df.geo %>% filter(unc != '') %>%
            group_by(unc.cat) %>%
            summarise(n=n(), .groups = 'drop'))

    results[['world.numbers']]$georeferenced.by.degree <- bind_rows(
        results[['world.numbers']]$georeferenced.by.degree,
        df.geo %>%
            dplyr::select(lon, lat) %>%
            mutate(lon=floor(as.double(lon)) + 0.5, lat=floor(as.double(lat)) + 0.5) %>%
            group_by(lon, lat) %>%
            summarise(n = n(), .groups = 'drop'))

    results[['world.numbers']]$georeferenced.with.unc.by.degree <- bind_rows(
        results[['world.numbers']]$georeferenced.with.unc.by.degree,
        df.geo %>%
            filter(unc != '') %>%
            dplyr::select(lon, lat) %>%
            mutate(lon=floor(as.double(lon)) + 0.5, lat=floor(as.double(lat)) + 0.5) %>%
            group_by(lon, lat) %>%
            summarise(n = n(), .groups = 'drop'))

    # COUNTRY-LEVEL NUMBERS ---------------------------------------------------------------------- #
    results[['country.numbers']]$country.georeferenced.all <- bind_rows(
        results[['country.numbers']]$country.georeferenced.all,
        df.geo %>%
            group_by(country.code) %>%
            summarise(n=n(), .groups = 'drop'))

    results[['country.numbers']]$country.georeferenced.unc.percentages <- bind_rows(
        results[['country.numbers']]$country.georeferenced.unc.percentages,
        df.geo %>%
            dplyr::select(country.code, unc) %>%
            mutate(unc = as.double(unc)) %>%
            mutate(has.unc = ifelse(is.na(unc), 'N', 'Y')) %>%
            group_by(country.code, has.unc) %>%
            summarise(n=n(), .groups = 'drop'))

    results[['country.numbers']]$country.georeferenced.unc <- bind_rows(
        results[['country.numbers']]$country.georeferenced.unc,
        df.geo %>% filter(unc != '') %>%
            group_by(country.code) %>%
            summarise(n=n(), .groups = 'drop'))

    results[['country.numbers']]$country.georeferenced.unc.cat <- bind_rows(
        results[['country.numbers']]$country.georeferenced.unc.cat,
        df.geo %>% filter(unc != '') %>%
            group_by(country.code, unc.cat) %>%
            summarise(n=n(), .groups = 'drop'))

    # TAXONOMY-LEVEL NUMBERS --------------------------------------------------------------------- #
    results[['taxonomy.numbers']]$family.georeferenced.all <- bind_rows(
        results[['taxonomy.numbers']]$family.georeferenced.all,
        df.geo %>%
            mutate(unc = as.double(unc),
                   has.unc = ifelse(is.na(unc), 'N', 'Y')) %>%
            group_by(kingdom, phylum, class, order, family, has.unc) %>%
            summarise(n=n(), .groups = 'drop'))

    results[['taxonomy.numbers']]$sci.name.georeferenced.n <- bind_rows(
        results[['taxonomy.numbers']]$sci.name.georeferenced.n,
        df.geo %>%
            dplyr::select(kingdom, phylum, class, order, family, sci.name, unc) %>%
            mutate(unc=ifelse(unc == '', 'N', "Y")) %>%
            rename(has.unc=unc) %>%
            group_by(kingdom, phylum, class, order, family, sci.name, has.unc) %>%
            summarise(n=n(), .groups = 'drop'))

    results[['taxonomy.numbers']]$kingdom.georeferenced.unc.cat <- bind_rows(
        results[['taxonomy.numbers']]$kingdom.georeferenced.unc.cat,
        df.geo %>% filter(unc != '') %>%
            group_by(kingdom, unc.cat) %>%
            summarise(n=n(), .groups = 'drop'))

    # PUBLISHER AND DATASET NUMBERS -------------------------------------------------------------- #
    results[['institution.numbers']]$inst.country.taxonomy.geo <- bind_rows(
        results[['institution.numbers']]$inst.country.taxonomy.geo,
        df.geo %>%
            mutate(lon=floor(as.double(lon)) + 0.5, lat=floor(as.double(lat)) + 0.5) %>%
            mutate(has.unc = ifelse(is.na(unc), 'N', 'Y')) %>%
            group_by(pub.country, pub.key, ds.name, kingdom, phylum, class, country.code, lon, lat, has.unc) %>%
            summarise(n=n(), .groups = 'drop'))

    offset <- offset + nrows
    n <- n + nrow(df.geo)
    cat("    ==>", formatC(n, big.mark = ' ', digits = 0, format = 'f'), "records processed so far ...\n")
}
dbDisconnect(con)

# Finally compute overall final results from partial results above
cat(rep("-", 80), "\n\n", sep = "")
cat("  summarising final results from partial results ...\n")
results[['summary.numbers']]$summary.taxonomy.country.dataset.geo <-
    results[['summary.numbers']]$summary.taxonomy.country.dataset.geo %>%
    group_by(kingdom, phylum, class, country.code, ds.name, pub.country, pub.key) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['world.numbers']]$georeferenced.unc.cat <-
    results[['world.numbers']]$georeferenced.unc.cat %>%
    group_by(unc.cat) %>%
    summarise(n=sum(n)) %>%
    mutate(unc.cat = factor(unc.cat,
                            levels = c("<= 1m", "10m", "100m", "250m", "1km", "5km", "10km", "50km", "100km", "> 100km"))) %>%
    arrange(unc.cat) %>%
    mutate(cs = cumsum(n), cum.perc=round(cs/sum(n)*100,1), .groups = 'drop')

results[['world.numbers']]$georeferenced.by.degree <-
    results[['world.numbers']]$georeferenced.by.degree %>%
    group_by(lon, lat) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['world.numbers']]$georeferenced.with.unc.by.degree <-
    results[['world.numbers']]$georeferenced.with.unc.by.degree %>%
    group_by(lon, lat) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['event.numbers']]$event.year.georeferenced.all <-
    results[['event.numbers']]$event.year.georeferenced.all %>%
    group_by(country.code, kingdom, phylum, class, event.year) %>%
    summarise(n.geo = sum(n.geo), n.geo.unc = sum(n.geo.unc), .groups = 'drop')

# Weighted mean returns the correct mean of means, as if done with all data at once and
# the following function returns the SD out of a group of samples SDs
grand.sd = function(N, M, S) {
    Mbar = weighted.mean(M, N)
    SSw = (N-1)*S^2
    SSb = N*(M-Mbar)^2
    SSbtwn = sum(SSb)
    SSwithin = sum(SSw)
    SStot = SSbtwn + SSwithin
    SD =sqrt(SStot/(sum(N)-1))
    return(SD)
}
results[['event.numbers']]$event.year.georeferenced.unc <-
    results[['event.numbers']]$event.year.georeferenced.unc %>%
    group_by(country.code, kingdom, phylum, class, event.year) %>%
    summarise(mean.unc=weighted.mean(mean.unc, n), .groups = 'drop') %>%
    inner_join(results[['event.numbers']]$event.year.georeferenced.unc %>%
                   group_by(country.code, kingdom, phylum, class, event.year) %>%
                   summarise(sd.unc=grand.sd(n, mean.unc, sd.unc), .groups = 'drop'),
               by=c("country.code", "kingdom", "phylum", "class", "event.year"))

results[['country.numbers']]$country.georeferenced.all <-
    results[['country.numbers']]$country.georeferenced.all %>%
    group_by(country.code) %>%
    summarise(n=sum(n)) %>%
    mutate(n.cat = cut(n, breaks = c(-Inf, 10^3, 10^4, 10^5, 10^6, 10^7, Inf),
                       labels = c("<= 1", "(1, 10]", "(10, 100]", "(100, 1 000]", "(1 000, 10 000]", "> 10 000"),
                       include.lowest = T)) %>%
    ungroup()

results[['country.numbers']]$country.georeferenced.unc.percentages <-
    results[['country.numbers']]$country.georeferenced.unc.percentages %>%
    group_by(country.code, has.unc) %>%
    summarise(n=sum(n), .groups = 'drop') %>%
    dplyr::select(country.code, has.unc, n) %>%
    tidyr::spread(has.unc, n) %>%
    setNames(c("country.code", "no.unc", "unc")) %>%
    mutate(georef = no.unc + unc,
           perc.coords = round(georef / results[['world.numbers']]$georeferenced.all * 100, 1),
           perc.unc = round(unc / georef * 100, 1),
           perc.unc.cat = cut(perc.unc, breaks = c(-Inf, 25, 50, 75, 100),
                              labels = c("< 25", "(25 -50]", "(50, 75]", "(75, 100]")))

results[['country.numbers']]$country.georeferenced.unc <-
    results[['country.numbers']]$country.georeferenced.unc %>%
    group_by(country.code) %>%
    summarise(n=sum(n)) %>%
    mutate(n.cat = cut(n, breaks = c(-Inf, 10^3, 10^4, 10^5, 10^6, 10^7, Inf),
                       labels = c("<= 1", "(1, 10]", "(10, 100]", "(100, 1 000]", "(1 000, 10 000]", "> 10 000"),
                       include.lowest = T)) %>%
    ungroup()

results[['country.numbers']]$country.georeferenced.unc.cat <-
    results[['country.numbers']]$country.georeferenced.unc.cat %>%
    group_by(country.code, unc.cat) %>%
    summarise(n=sum(n), .groups = 'drop') %>%
    arrange(country.code, desc(n)) %>%
    mutate(cs = cumsum(n), cum.perc=round(cs/sum(n)*100,1))

results[['taxonomy.numbers']]$family.georeferenced.all <-
    results[['taxonomy.numbers']]$family.georeferenced.all %>%
    group_by(kingdom, phylum, class, order, family, has.unc) %>%
    summarise(n=sum(n), .groups = 'drop') %>%
    tidyr::spread(has.unc, n) %>%
    # Following mutate is because there may be records without georeference or without uncertainty
    mutate(Y = ifelse(is.na(Y), 0, Y),
           N = ifelse(is.na(N), 0, N)) %>%
    rename("non.unc"='N', 'unc'='Y') %>%
    mutate(georef = non.unc + unc,
           perc.coords = round(georef / results[['world.numbers']]$georeferenced.all * 100, 1),
           perc.unc = round(unc / georef * 100, 1)) %>%
    arrange(desc(georef))

results[['taxonomy.numbers']]$sci.name.georeferenced.n <-
    results[['taxonomy.numbers']]$sci.name.georeferenced.n %>%
    group_by(kingdom, phylum, class, order, family, sci.name, has.unc) %>%
    summarise(n=sum(n), .groups = 'drop')

unc.cat.labels <- c("<= 1m", "10m", "100m", "250m", "1km", "5km", "10km", "50km", "100km", "> 100km")
results[['taxonomy.numbers']]$kingdom.georeferenced.unc.cat <-
    results[['taxonomy.numbers']]$kingdom.georeferenced.unc.cat %>%
    mutate(unc.cat = factor(unc.cat, levels = unc.cat.labels)) %>%
    group_by(kingdom, unc.cat) %>%
    summarise(n=sum(n), .groups = 'drop') %>%
    arrange(kingdom, unc.cat) %>%
    mutate(cs = cumsum(n), cum.perc=round(cs/sum(n)*100,1))

# PUBLISHER AND DATASET NUMBERS ------------------------------------------------------------------ #
results[['institution.numbers']]$inst.country.taxonomy.geo <-
    results[['institution.numbers']]$inst.country.taxonomy.geo %>%
        group_by(pub.country, pub.key, ds.name, kingdom, phylum, class, country.code, lon, lat, has.unc) %>%
        summarise(n=sum(n), .groups = 'drop')

cat("\n  saving results ...\n")
saveRDS(results, f.out.geo.whole.results)
cat("    ==>", f.out.geo.whole.results, "\n")
cat(rep("-", 80), "\n", sep = "")
