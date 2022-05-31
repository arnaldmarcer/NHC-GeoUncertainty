# This script reads the non-geo database and computes partial results so that
# scripts run faster down the line

source("scripts/R/init.R")
dir.create(output.dir, recursive = T, showWarnings = F)

results <- list()
results[['world.numbers']]$non.georeferenced <- 0

nrows <- 1000000 # Process in chunks of this size
cat(rep("=", 80), "\n", sep = "")
cat("Generating results for non-geo dataset ...\n")
cat("  getting partial results in", formatC(nrows, big.mark = ' ', digits = 0, format = 'f'), "record chunks\n")
con <- dbConnect(RSQLite::SQLite(), f.db.whole.non.geo.processed)
offset <- 0
n <- 0
repeat{
    sql.rows <- paste0("select * from ps", " limit ", nrows, " offset ", offset)
    df.non.geo <- dbGetQuery(con, sql.rows)
    names(df.non.geo) <- gsub("\\_", ".", names(df.non.geo))

    if(nrow(df.non.geo) == 0) # run out of records to process
        break

    results[['world.numbers']]$non.georeferenced <-
        results[['world.numbers']]$non.georeferenced + nrow(df.non.geo)

    # SUMMARY DATA ------------------------------------------------------------------------------- #
    # Number of records by kingdom, phylum, class, country.code, ds.name, inst.code
    results[['summary.numbers']]$summary.taxonomy.country.dataset.non.geo <- bind_rows(
        results[['summary.numbers']]$summary.taxonomy.country.dataset.non.geo,
        df.non.geo %>%
            group_by(kingdom, phylum, class, country.code, ds.name, inst.code) %>%
            summarise(n = n(), .groups = 'drop'))

    # EVENT DATE --------------------------------------------------------------------------------- #
    results[['event.numbers']]$event.year.non.georeferenced <- bind_rows(
        results[['event.numbers']]$event.year.non.georeferenced,
        df.non.geo %>%
            group_by(country.code, kingdom, phylum, class, event.year) %>%
            count())

    # BY COUNTRY --------------------------------------------------------------------------------- #
    # Number of non.georeferenced records by country
    # Non-geo: country, n
    results[['country.numbers']]$country.non.georeferenced <- bind_rows(
        results[['country.numbers']]$country.non.georeferenced,
        df.non.geo %>%
            group_by(country.code) %>%
            summarise(n=n(), .groups = 'drop'))

    # BY TAXONOMY -------------------------------------------------------------------------------- #
    # Number of non.georeferenced records by taxonomy, down to scientific name
    # has.unc and subset to be able to bind rows with georeferenced counterpart below
    # Non-geo: kingdom, phylum, class, order, family, has.unc, n, subset
    results[['taxonomy.numbers']]$sci.name.non.geo <- bind_rows(
        results[['taxonomy.numbers']]$sci.name.non.geo,
        df.non.geo %>%
            group_by(kingdom, phylum, class, order, family, sci.name) %>%
            summarise(n=n(), .groups = 'drop'))

    # BY PUBLISHING INSTITUTION ------------------------------------------------------------------ #
    results[['institution.numbers']]$inst.non.geo <- bind_rows(
        results[['institution.numbers']]$inst.non.geo,
        df.non.geo %>%
            group_by(pub.country, pub.key, ds.name, country.code, kingdom, phylum, class, event.year) %>%
            summarise(n=n(), .groups = 'drop'))

    offset <- offset + nrows
    n <- n + nrow(df.non.geo)
    cat("    ==>", formatC(n, big.mark = ' ', digits = 0, format = 'f'), "records processed so far ...\n")
}
dbDisconnect(con)
cat(rep("-", 80), "\n\n", sep = "")

cat(rep("=", 80), "\n", sep = "")
cat("  summarising final results from partial results ...\n")
results[['summary.numbers']]$summary.taxonomy.country.dataset.non.geo <-
    results[['summary.numbers']]$summary.taxonomy.country.dataset.non.geo %>%
    group_by(kingdom, phylum, class, country.code, ds.name, inst.code) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['event.numbers']]$event.year.non.georeferenced <-
    results[['event.numbers']]$event.year.non.georeferenced %>%
    group_by(country.code, kingdom, phylum, class, event.year) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['country.numbers']]$country.non.georeferenced <-
    results[['country.numbers']]$country.non.georeferenced %>%
    group_by(country.code) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['taxonomy.numbers']]$sci.name.non.geo <-
    results[['taxonomy.numbers']]$sci.name.non.geo %>%
    group_by(kingdom, phylum, class, order, family, sci.name) %>%
    summarise(n = sum(n), .groups = 'drop')

results[['institution.numbers']]$inst.non.geo <-
    results[['institution.numbers']]$inst.non.geo %>%
        group_by(pub.country, pub.key, ds.name, country.code, kingdom, phylum, class, event.year) %>%
        summarise(n = sum(n), .groups = 'drop')

# SAVE RESULTS ----------------------------------------------------------------------------------- #
cat("\n  saving results ...\n")
saveRDS(results, f.out.non.geo.whole.results)
cat("    ==>", f.out.non.geo.whole.results, "\n")
cat(rep("-", 80), "\n\n", sep = "")

