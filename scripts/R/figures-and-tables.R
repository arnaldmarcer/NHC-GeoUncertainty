# This scripts creates the table and figures of the manuscript and writes results in file results.dat
source("scripts/R/init.R")

gf.col <- c("NC"="#BBBBBB", "CO"="#4477AA", "CU"="#AA3477")
gt.col <- c("NC"="#8A8A8A", "CO"="#4477AA", "CU"="#AA3477")

plot_label_maker = function(breaks,unit_MK,unit_scale) {
    labels = scales::unit_format(unit = unit_MK, scale = unit_scale,accuracy = 1,sep="")(breaks) %>%
        stringr::str_replace_all(stringr::regex("^0M$"),"0") %>%
        stringr::str_replace_all(stringr::regex("^0K$"),"0") %>%
        stringr::str_replace_all(stringr::regex("^1 000K$"),"1M") %>%
        stringr::str_replace_all(stringr::regex("^2 000K$"),"2M") %>%
        stringr::str_replace_all(stringr::regex("^3 000K$"),"3M") %>%
        stringr::str_replace_all(stringr::regex("^4 000K$"),"4M") %>%
        stringr::str_replace_all(stringr::regex("^5 000K$"),"5M") %>%
        stringr::str_replace_all(stringr::regex("^1 000M$"),"1B") %>%
        stringr::str_replace_all(stringr::regex("^2 000M$"),"2B")
    return(labels)
}

# Read results files
cat(rep("=", 80), "\n", sep = "")
cat("Reading intermediate results files ...\n")
results.non.geo <- readRDS(paste0(output.dir, date.string, "_non_geo_whole_results.rds"))
results.geo <- readRDS(paste0(output.dir, date.string, "_geo_whole_results.rds"))
results <- mapply(c, results.geo, results.non.geo, SIMPLIFY = FALSE)

#### MAIN TEXT ================================================================================ ####
#### TABLE 1 ---------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Generating Table 1 data ...\n")
ps.total <- results.non.geo$world.numbers$non.georeferenced + results.geo$world.numbers$georeferenced.all
ps.nc <- results.non.geo$world.numbers$non.georeferenced
ps.co <- results.geo$world.numbers$georeferenced.all - results.geo$world.numbers$georeferenced.unc
ps.cu <- results.geo$world.numbers$georeferenced.unc
writeResult("table-1-number-preserved-specimens-total", ps.total)
writeResult("table-1-number-preserved-specimens-without-coordinates", ps.nc)
writeResult("table-1-number-preserved-specimens-coordinates-only", ps.co)
writeResult("table-1-number-preserved-specimens-coordinates-with-uncertainty", ps.cu)
writeResult("table-1-percentage-NC-over-total", round(ps.nc/ps.total * 100, 2))
writeResult("table-1-percentage-CO-over-total", round(ps.co/ps.total * 100, 2))
writeResult("table-1-percentage-CU-over-total", round(ps.cu/ps.total * 100, 2))

#### RESULTS' TEXT ---------------------------------------------------------------------------- ####
f.data <- paste0(data.processed.dir, "gbif/", date.string, "_coordinate_uncertainty_snapshots.tsv")
d <- readr::read_tsv(f.data, col_types = cols()) %>%
    dplyr::select(date, n.g=count_total, n.ngu=count_null, nps=total_specimens) %>%
    mutate(n.ng = nps - n.g, n.gu = n.g - n.ngu) %>%
    mutate(ratio.g.gu = round(n.gu / n.g, 4),
           ratio.gu.ngu = round(n.gu / n.ngu, 4)) %>%
    tidyr::pivot_longer(-date, names_to = "subset", values_to = "value") %>%
    mutate(subset = factor(subset, levels = c("n.ng", "n.g", "n.ngu", "n.gu", "nps", "ratio.g.gu", "ratio.gu.ngu"),
                           labels = c("NC", "C",
                                      "CO", "CU",
                                      "Nr of specimens",
                                      "Ratio CU/C", "Ratio CU/CO")))

breaks <- seq(0, 1.8*10^8, 2*10^7)
labels = plot_label_maker(breaks,unit_MK="M",unit_scale=1e-6)

p <- ggplot() +
    geom_area(data=d %>% filter(subset %in% c("NC", "CO", "CU")),
              aes(x=date, y=value, fill=subset), position = "stack", alpha = 0.75) +
    geom_line(data=d %>% filter(subset == 'Ratio CU/C'),
              aes(x=date, y=value * 10^8), colour="blue", linetype = 'dashed') +
    geom_line(data=d %>% filter(subset == 'Ratio CU/CO'),
              aes(x=date, y=value * 10^8), colour = "red", linetype = 'dashed') +
    scale_y_continuous("Number of specimens", expand = c(0,0), breaks=breaks, labels=labels, limits=c(0, max(breaks)),
                       sec.axis = sec_axis(~ ./10^8, name="", breaks=seq(0, 1, 0.1))) +
    scale_x_date("Year", date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    scale_fill_manual(values = gf.col[c("NC", "CO", "CU")],
                      labels = c("Records without coordinates","Records with coordinates (CO)","Records with uncertainty (CU)")) +
    # scale_colour_manual("Ratio CU/(CO + CU)", values = c("Ratio CU/C"='darkblue'), guide = F) +
    theme(panel.border = element_blank(), legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.position=c(0.25,0.90),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 5),
          axis.ticks = element_line(size = 0.2),
          axis.title.y.right = element_text(hjust = 0.85, vjust = 0.8),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 4))
p
f.out <- paste0(manuscript.dir, "figures/figure_1.tif")
ggsave(plot=p, device="tiff", filename = f.out, width = 9, height = 7, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Data on records with coordinates in January 2015
nr.coords.2015 <- d %>% arrange(date) %>%
    filter(date >= "2015-01-01" & date <= "2015-01-31" & (subset == 'CO' | subset == 'CU')) %>%
    pull(value) %>% sum()
nr.coords.2015
nr.all.2015 <- d %>% arrange(date) %>%
    filter(date >= "2015-01-01" & date <= "2015-01-31" & subset == 'Nr of specimens') %>%
    pull(value)
nr.all.2015
perc.cooords.2015 <- round(nr.coords.2015/nr.all.2015 * 100, 2)
perc.cooords.2015
writeResult("results-text-nr-coords-2015", nr.coords.2015)
writeResult("results-text-nr-all-2015", nr.all.2015)
writeResult("results-text-percentage-coords-2015", perc.cooords.2015)

# Data on records with coordinates in March 2021
nr.coords.2021 <- d %>% arrange(date) %>%
    filter(date >= "2021-03-01" & date <= "2021-03-31" & (subset == 'CO' | subset == 'CU')) %>%
    pull(value) %>% sum()
nr.coords.2021
nr.all.2021 <- d %>% arrange(date) %>%
    filter(date >= "2021-03-01" & date <= "2021-03-31" & subset == 'Nr of specimens') %>%
    pull(value)
nr.all.2021
perc.cooords.2021 <- round(nr.coords.2021/nr.all.2021 * 100, 2)
perc.cooords.2021
writeResult("results-text-nr-coords-2021", nr.coords.2021)
writeResult("results-text-nr-all-2021", nr.all.2021)
writeResult("results-text-percentage-coords-2021", perc.cooords.2021)

# Data on records with coordinates and uncertainty in January 2015
nr.unc.2015 <- d %>% arrange(date) %>%
    filter(date >= "2015-01-01" & date <= "2015-01-31" & (subset == 'CU')) %>%
    pull(value)
nr.unc.2015
perc.unc.2015 <- round(nr.unc.2015/nr.coords.2015 * 100, 2)
perc.unc.2015
writeResult("results-text-nr-unc-2015", nr.unc.2015)
writeResult("results-text-percentage-unc-2015", perc.unc.2015)

# Data on records with coordinates and uncertainty in January 2021
nr.unc.2021 <- d %>% arrange(date) %>%
    filter(date >= "2021-03-01" & date <= "2021-03-31" & (subset == 'CU')) %>%
    pull(value)
nr.unc.2021
perc.unc.2021 <- round(nr.unc.2021/nr.coords.2021 * 100, 2)
perc.unc.2021
writeResult("results-text-nr-unc-2021", nr.unc.2021)
writeResult("results-text-percentage-unc-2021", perc.unc.2021)


#### FIGURE 1 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting Figure 1 ...\n")

breaks <- seq(0, 2*10^8, 2*10^7)
labels = plot_label_maker(breaks,unit_MK="M",unit_scale=1e-6)

p <- ggplot() +
    geom_area(data=d %>% filter(subset %in% c("NC", "CO", "CU")),
              aes(x=date, y=value, fill=subset), position = "stack", alpha = 0.75) +
    geom_line(data=d %>% filter(subset == 'Ratio CU/C'),
              aes(x=date, y=value * 10^8, colour=subset), linetype = 'dashed') +
    scale_y_continuous("Number of specimens", expand = c(0,0), breaks=breaks, labels=labels, limits=c(0, max(breaks)),
                       sec.axis = sec_axis(~ ./10^8, name="Ratio CU/(CO + CU)", breaks=seq(0, 0.75, 0.1))) +
    scale_x_date("Year", date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    scale_fill_manual(values = gf.col[c("NC", "CO", "CU")],
                      labels = c("Records without coordinates","Records with coordinates (CO)","Records with uncertainty (CU)")) +
    scale_colour_manual("Ratio CU/(CO + CU)", values = c("Ratio CU/C"='darkblue'), guide = F) +
    theme(panel.border = element_blank(), legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.position=c(0.30,0.85),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 5),
          axis.ticks = element_line(size = 0.2),
          axis.title.y.right = element_text(hjust = 0.85, vjust = 0.8),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 4))
p
f.out <- paste0(manuscript.dir, "figures/figure_1.png")
ggsave(plot=p, device="png", filename = f.out, width = 9, height = 7, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

#### FIGURE 2 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting Figure 2 ...\n")

perc.interval <- 20
df.coords <- results$world.numbers$georeferenced.by.degree %>%
    rename("n.geo" = "n") %>%
    left_join(results$world.numbers$georeferenced.with.unc.by.degree %>%
                  rename("n.unc" = "n"), by=c("lon", "lat")) %>%
    group_by(lon, lat) %>%
    mutate(n.geo = ifelse(is.na(n.geo), 0, n.geo),
           n.unc = ifelse(is.na(n.unc), 0, n.unc)) %>%
    summarise(n.geo = sum(n.geo),
              n.unc = sum(n.unc), .groups = 'drop') %>%
    mutate(unc.perc = round(n.unc / n.geo * 100, 1)) %>%
    mutate(unc.perc.c = as.character(ifelse(unc.perc == 100, 100/perc.interval - 1, floor(unc.perc/perc.interval))))

# Fig 2a
breaks <- c(1, 1 %o% 10^(1:floor(log10(max(df.coords$n.geo)))))
p.geo <- ggplot() +
    geom_tile(data = df.coords, aes(x = lon, y = lat, fill = n.geo)) +
    geom_sf(data = world_coastline, fill = NA, colour = 'lightyellow2', size = 0.4) +
    scale_fill_viridis(trans = "log",
                       option = 'viridis',
                       breaks = breaks,
                       labels = formatC(breaks, digits=0, format="e")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_colourbar(title = "Number of records with coordinates",
                                  title.position = "top",
                                  label.theme = element_text(hjust=0.5, size=8),
                                  title.theme = element_text(size=9, hjust = 0.5),
                                  barwidth = unit(7, "cm"),
                                  barheight = unit(0.25, "cm")))
# p.geo
f.out <- paste0(manuscript.dir, "figures/figure_2a.tif")
ggsave(plot=p.geo, device="tiff", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig 2b
unc.cat.colors <- c("#404040", "#ca0020", "yellow", "green", "blue")
p.unc <- ggplot() +
    geom_tile(data = df.coords, aes(x = lon, y = lat, fill = unc.perc.c, alpha = scales::rescale(n.geo))) +
    geom_sf(data = world_coastline, fill = NA, colour = 'grey40', size = 0.4) +
    scale_fill_manual(values = unc.cat.colors, labels = c("[0-20)", "[20-40)", "[40-60)", "[60-80)", "[80-100]")) +
    scale_alpha_continuous(guide = 'none', trans = 'log', breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank(), legend.title = element_text(size = 7),
        legend.text = element_text(size = 6)) +
    guides(
        fill = guide_legend(title.position = "top", order = 1, title = "% records with uncertainty",
                            label.theme = element_text(hjust=0.5, size=8),
                            title.theme = element_text(size=9, hjust = 0.5))
    )
# p.unc
f.out <- paste0(manuscript.dir, "figures/figure_2b.png")
ggsave(plot=p.unc, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

#### FIGURE 3 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure 3 ...\n")

country.cum.perc.threshold <- 80
df.bar.country <- results[['country.numbers']]$country.non.georeferenced %>%
    rename("non.georef"="n") %>%
    inner_join(results[['country.numbers']]$country.georeferenced.unc.percentages %>%
                   filter(!is.na(georef)), by = c("country.code")) %>%
    mutate(n.total = non.georef + georef) %>%
    mutate(perc.geo = round(georef/n.total * 100, 1)) %>%
    arrange(desc(n.total)) %>%
    mutate(perc.total = round(n.total/sum(across(n.total)) * 100, 1)) %>%
    mutate(perc.unc.total = round(unc/n.total * 100, 1)) %>%
    mutate(perc.cum = cumsum(perc.total)) %>%
    slice(1:which(perc.cum >= country.cum.perc.threshold)[1]) %>%
    inner_join(world %>%
                   dplyr::select(iso_a2, name) %>%
                   st_drop_geometry() %>%
                   mutate(name=ifelse(name == 'Korea', 'South Korea', name)), by = c("country.code"="iso_a2")) %>%
    dplyr::select(country.code, name, n.total, non.georef, georef, perc.geo, no.unc, unc, perc.unc, perc.unc.total, perc.cum)

writeResult("figure-3a-cumulative-percentage-threshold-used-to-select-countries", country.cum.perc.threshold)
writeResult("figure-3a-number-of-countries-represented", df.bar.country %>% pull(country.code) %>% length())

# Fig 3a
leg.x <- 6
leg.y <- 18*10^6
p.bar.country <- ggplot() +
    # total specimens
    geom_bar(data = df.bar.country, aes(x = reorder(name, n.total), y = n.total), fill = gf.col["NC"], stat = 'identity') +
    # With coordinates
    geom_bar(data = df.bar.country, aes(x = reorder(name, georef), y = georef), fill = gf.col["CO"], stat = 'identity') +
    # With coordinates and uncertainty
    geom_bar(data = df.bar.country, aes(x = reorder(name, georef), y = unc), fill = gf.col["CU"], stat = 'identity') +
    annotate('text',
             x = 1:24,
             y = -6*10^6,
             label = paste0("[", formatC(df.bar.country %>% arrange(n.total) %>% pull(perc.unc.total),
                                         flag = '0', format = 'f', digits=1, width=4), "%, "),
             colour = gt.col["CU"], size = 4, fontface = "bold") +
    annotate('text',
             x = 1:24,
             y = -3.5*10^6,
             label = paste0(formatC(df.bar.country %>% arrange(n.total) %>% pull(perc.geo), format = 'f', digits=1, width=4), "%] "),
             colour = gt.col["CO"], size = 4, fontface = "bold") +
    annotate('text',
             x = 1:24,
             y = -1.5*10^6,
             label = paste0(formatC(df.bar.country %>% arrange(n.total) %>% pull(n.total) / 10^6, format = 'f', digits=1, width=4), "M"),
             colour = gt.col["NC"], size = 4, fontface = "bold") +
    scale_y_continuous("Number of specimens (in millions)", breaks = seq(0, 35000000, 5000000),
                       labels = c("0M", "5M", "10M", "15M", "20M", "25M", "30M", "35M"),
                       limits=c(-7.5*10^6,36000000), expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    # white background
    geom_rect(aes(xmin=leg.x - 4, ymin=leg.y - 10^6, xmax=leg.x + 2, ymax=leg.y + 15.5^6), fill="white") +
    # Total
    geom_rect(aes(xmin=leg.x, ymin=leg.y, xmax=leg.x + 1, ymax=leg.y + 10^6), fill=gf.col["NC"]) +
    annotate("text", x=leg.x + 0.5, y=leg.y + 1500000, label="Preserved specimens",
             hjust = 0, size=4) +
    # Georef
    geom_rect(aes(xmin=leg.x - 0.5, ymin=leg.y, xmax=leg.x - 1.5, ymax=leg.y + 10^6), fill=gf.col["CO"]) +
    annotate("text", x=leg.x - 1, y=leg.y + 1500000, label="Records with coordinates",
             hjust = 0, size=4) +
    # Uncertainty
    geom_rect(aes(xmin=leg.x - 2, ymin=leg.y, xmax=leg.x - 3, ymax=leg.y + 10^6), fill=gf.col["CU"]) +
    annotate("text", x=leg.x - 2.5, y=leg.y + 1500000, label="Records with uncertainty",
             hjust = 0, size=4) +
    theme(axis.text.x = element_text(angle = 0, size=11),
          panel.grid.major.x = element_line(colour="grey70", linetype = "dashed", size = 0.2),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          axis.text.y = element_text(hjust=1, size=11),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(0.8, 0.2))
p.bar.country
f.out <- paste0(manuscript.dir, "figures/figure_3a.png")
ggsave(plot=p.bar.country, device="png", filename = f.out, width = 28, height = 16, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig 3b
publishing.country.cum.perc.threshold <- 80
df.bar.publishing.country <- inner_join(results$institution.numbers$inst.non.geo %>%
                                            group_by(pub.country) %>%
                                            summarise(non.georef = sum(n), .groups = "drop"),
                                        results$institution.numbers$inst.country.taxonomy.geo %>%
                                            group_by(pub.country, has.unc) %>%
                                            summarise(n = sum(n), .groups = "drop") %>%
                                            pivot_wider(names_from = has.unc, values_from = n) %>%
                                            mutate(N = ifelse(is.na(N), 0, N),
                                                   Y = ifelse(is.na(Y), 0, Y)) %>%
                                            mutate(georef = N + Y) %>%
                                            dplyr::select(pub.country, georef, unc=Y), by = c("pub.country")) %>%
    mutate(n.total = non.georef + georef) %>%
    mutate(perc.geo = round(georef/n.total * 100, 1),
           no.unc = georef - unc,
           perc.unc = round(unc/georef * 100, 1),
           perc.unc.total = round(unc / (georef + non.georef) * 100, 1),
           perc.total = round(n.total/sum(across(n.total)) * 100, 1)) %>%
    arrange(desc(n.total)) %>%
    mutate(perc.cum = cumsum(perc.total)) %>%
    inner_join(world %>% st_drop_geometry %>% dplyr::select(name, pub.country = iso_a2), by=c("pub.country")) %>%
    slice(1:which(perc.cum >= publishing.country.cum.perc.threshold)[1])

writeResult("figure-3b-cumulative-percentage-threshold-used-to-select-countries", publishing.country.cum.perc.threshold)
writeResult("figure-3b-number-of-countries-represented", nrow(df.bar.publishing.country))

writeResult("publishing_country_maximum_number_of_records",
            df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(1) %>% pull(name))
writeResult("figure-3b-publishing_country_maximum_number_of_records_percentage",
            df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(1) %>% pull(perc.total))
writeResult("figure-3b-publishing_country_maximum_number_of_records_number",
            formatC(df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(1) %>% pull(n.total),
                    format = "f", digits = 0, big.mark = " "))

writeResult("figure-3b-publishing_country_second_maximum_number_of_records",
            df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(2) %>% pull(name))
writeResult("figure-3b-publishing_country_second_maximum_number_of_records_percentage",
            df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(2) %>% pull(perc.total))
writeResult("figure-3b-publishing_country_second_maximum_number_of_records_number",
            formatC(df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(2) %>% pull(n.total),
                    format = "f", digits = 0, big.mark = " "))

writeResult("figure-3b-publishing_country_third_maximum_number_of_records",
            df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(3) %>% pull(name))
writeResult("figure-3b-publishing_country_third_maximum_number_of_records_percentage",
            df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(3) %>% pull(perc.total))
writeResult("figure-3b-publishing_country_third_maximum_number_of_records_number",
            formatC(df.bar.publishing.country %>% arrange(desc(n.total)) %>% slice(3) %>% pull(n.total),
                    format = "f", digits = 0, big.mark = " "))

writeResult("figure-3b-publishing_country_maximum_number_of_records_georef",
            df.bar.publishing.country %>% arrange(desc(perc.geo)) %>% slice(1) %>% pull(name))
writeResult("figure-3b-publishing_country_second_maximum_number_of_records_georef",
            df.bar.publishing.country %>% arrange(desc(perc.geo)) %>% slice(2) %>% pull(name))
writeResult("figure-3b-publishing_country_third_maximum_number_of_records_georef",
            df.bar.publishing.country %>% arrange(desc(perc.geo)) %>% slice(3) %>% pull(name))

writeResult("figure-3b-publishing_country_maximum_number_of_records_unc",
            df.bar.publishing.country %>% arrange(desc(perc.unc)) %>% slice(1) %>% pull(name))
writeResult("figure-3b-publishing_country_second_maximum_number_of_records_unc",
            df.bar.publishing.country %>% arrange(desc(perc.unc)) %>% slice(2) %>% pull(name))
writeResult("figure-3b-publishing_country_third_maximum_number_of_records_unc",
            df.bar.publishing.country %>% arrange(desc(perc.unc)) %>% slice(3) %>% pull(name))

# Fig 3b
leg.x <- 5
leg.y <- 32.5*10^6
p.bar.publishing.country <- ggplot() +
    # total specimens
    geom_bar(data = df.bar.publishing.country, aes(x = reorder(name, n.total), y = n.total), fill = gf.col["NC"], stat = 'identity') +
    # Georeferenced
    geom_bar(data = df.bar.publishing.country, aes(x = reorder(name, georef), y = georef), fill = gf.col["CO"], stat = 'identity') +
    # Georeferenced with uncertainty
    geom_bar(data = df.bar.publishing.country, aes(x = reorder(name, georef), y = unc), fill = gf.col["CU"], stat = 'identity') +
    annotate('text',
             x = 1:13,
             y = -6*11^6,
             label = paste0("[", formatC(df.bar.publishing.country %>% arrange(n.total) %>% pull(perc.unc.total),
                                         flag = '0', format = 'f', digits=1, width=4), "%, "),
             colour = gt.col["CU"], size = 4, fontface = "bold") +
    annotate('text',
             x = 1:13,
             y = -3.5*11^6,
             label = paste0(formatC(df.bar.publishing.country %>% arrange(n.total) %>% pull(perc.geo), format = 'f', digits=1, width=4), "%] "),
             colour = gt.col["CO"], size = 4, fontface = "bold") +
    annotate('text',
             x = 1:13,
             y = -1.5*11^6,
             label = paste0(formatC(df.bar.publishing.country %>% arrange(n.total) %>% pull(n.total) / 10^6, format = 'f', digits=1, width=4), "M"),
             colour = gt.col["NC"], size = 4, fontface = "bold") +
    scale_y_continuous("Nr of specimens (in millions)", breaks = seq(0, 65000000, 5000000),
                       labels = c("0M", "5M", "10M", "15M", "20M", "25M", "30M", "35M", "40M", "45M", "50M", "55M", "60M", "65M"),
                       limits=c(-7.5*11^6,66000000), expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    theme(axis.text.x = element_text(size=11),
          panel.grid.major.x = element_line(colour="grey70", linetype = "dashed", size = 0.2),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          axis.text.y = element_text(hjust=1, size=11),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(0.8, 0.2))
p.bar.publishing.country
f.out <- paste0(manuscript.dir, "figures/figure_3b.png")
ggsave(plot=p.bar.publishing.country, device="png", filename = f.out, width = 30, height = 10, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

#### FIGURE 4 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure 4 ...\n")
df.kingdom <- results[['taxonomy.numbers']]$sci.name.non.geo %>%
    group_by(kingdom) %>%
    summarise(non.geo=sum(n), .groups = "drop") %>%
    inner_join(
        results[['taxonomy.numbers']]$family.georeferenced.all %>%
            group_by(kingdom) %>%
            summarise(geo=sum(georef), unc=sum(unc), .groups = "drop"), by=c("kingdom")) %>%
    mutate(total = non.geo + geo,
           perc.geo = round(geo / total * 100, 1),
           perc.unc = round(unc / total * 100, 1),
           kingdom = ifelse(kingdom == 'incertae sedis', 'UNDETERMINED', kingdom))

perc.an.pl <- round(df.kingdom %>% filter(kingdom %in% c('Animalia', 'Plantae')) %>% pull(total) %>% sum() / df.kingdom %>%
                        pull(total) %>% sum() * 100, 1)
writeResult("figure-4a-percentage_ps_animalia_and_plantae", formatC(perc.an.pl, flag = '0', format = 'f', digits=1))
perc.an <- round(df.kingdom %>% filter(kingdom == 'Animalia') %>% pull(total) / df.kingdom %>%
                     pull(total) %>% sum() * 100, 1)
writeResult("figure-4a-percentage_ps_animalia", formatC(perc.an, flag = '0', format = 'f', digits=1))
perc.pl <- round(df.kingdom %>% filter(kingdom == 'Plantae') %>% pull(total) / df.kingdom %>%
                     pull(total) %>% sum() * 100, 1)
writeResult("figure-4a-percentage_ps_plantae", formatC(perc.pl, flag = '0', format = 'f', digits=1))
perc.fu <- round(df.kingdom %>% filter(kingdom == 'Fungi') %>% pull(total) / df.kingdom %>%
                     pull(total) %>% sum() * 100, 1)
writeResult("figure-4a-percentage_ps_fungi", formatC(perc.fu, flag = '0', format = 'f', digits=1))
perc.an.geo <- df.kingdom %>% filter(kingdom == 'Animalia') %>% pull(perc.geo)
writeResult("figure-4a-percentage_ps_animalia_georeferenced", formatC(perc.an.geo, flag = '0', format = 'f', digits=1))

perc.an.unc <- df.kingdom %>% filter(kingdom == 'Animalia') %>% pull(perc.unc)
writeResult("figure-4a-percentage__ps_animalia_CU", formatC(perc.an.unc, flag = '0', format = 'f', digits=1))

perc.pl.geo <- df.kingdom %>% filter(kingdom == 'Plantae') %>% pull(perc.geo)
writeResult("figure-4a-percentage_ps_plantae_C", formatC(perc.pl.geo, flag = '0', format = 'f', digits=1))

perc.pl.unc <- df.kingdom %>% filter(kingdom == 'Plantae') %>% pull(perc.unc)
writeResult("figure-4a-percentage_ps_plantae_CU", formatC(perc.pl.unc, flag = '0', format = 'f', digits=1))

perc.fu.geo <- df.kingdom %>% filter(kingdom == 'Fungi') %>% pull(perc.geo)
writeResult("figure-4a-percentage_ps_fungi_C", formatC(perc.fu.geo, flag = '0', format = 'f', digits=1))

perc.fu.unc <- df.kingdom %>% filter(kingdom == 'Fungi') %>% pull(perc.unc)
writeResult("figure-4a-percentage_ps_fungi_CU", formatC(perc.fu.unc, flag = '0', format = 'f', digits=1))

# Fig 4a
leg.x <- 4
leg.y <- 40*10^6
p.bar.kingdom <- ggplot(df.kingdom) +
    # total
    geom_bar(aes(x = reorder(kingdom, total), y = total), fill = gf.col["NC"], stat = 'identity') +
    annotate('text',
             x = 1:9,
             y = -13*10^6,
             label = paste0("[", formatC(df.kingdom %>% arrange(total) %>% pull(perc.unc),
                                         flag = '0', format = 'f', digits=1, width=4), "%, "),
             colour = gt.col["CU"], size = 4, fontface = "bold") +
    # georef
    geom_bar(aes(x = reorder(kingdom, total), y = geo), fill = gf.col["CO"], stat = 'identity') +
    annotate('text',
             x = 1:9,
             y = -7.5*10^6,
             label = paste0(formatC(df.kingdom %>% arrange(total) %>% pull(perc.geo),
                                    flag = '0', format = 'f', digits=1, width=4), "%] "),
             colour = gt.col["CO"], size = 4, fontface = "bold") +
    # uncertainty
    geom_bar(aes(x = reorder(kingdom, total), y = unc), fill = gf.col["CU"], stat = 'identity') +
    annotate('text',
             x = 1:9,
             y = -2.6*10^6,
             label = paste0(formatC(df.kingdom %>% arrange(total) %>% pull(total) / 10^6,
                                    flag = '0', format = 'f', digits=1, width=4), "M"),
             colour = gt.col["NC"], size = 4, fontface = "bold") +
    # white background
    geom_rect(aes(xmin=leg.x -2.5, ymin=leg.y - 10^6, xmax = leg.x + 1, ymax = leg.y + 30*10^6), fill="white") +
    # # Total
    geom_rect(aes(xmin=leg.x - 0.5, ymin=leg.y, xmax=leg.x -0.5 + 0.5, ymax=leg.y + 3.5*10^6), fill=gf.col["NC"]) +
    annotate("text", x=leg.x - 0.5 + 0.3, y=leg.y + 4000000, label="Preserved specimens",
             vjust = 1, hjust = 0, size=4) +
    # # Georef
    geom_rect(aes(xmin=leg.x - 1.25, ymin=leg.y, xmax=leg.x - 1.25 + 0.5, ymax=leg.y + 3.5*10^6), fill=gf.col["CO"]) +
    annotate("text", x=leg.x - 1.25 + 0.3, y=leg.y + 4000000, label="Records with coordinates",
             vjust = 1, hjust = 0, size=4) +
    # # Uncertainty
    geom_rect(aes(xmin=leg.x - 2, ymin=leg.y, xmax=leg.x - 2 + 0.5, ymax=leg.y + 3.5*10^6), fill=gf.col["CU"]) +
    annotate("text", x=leg.x - 2 + 0.3, y=leg.y + 4000000, label="Records with uncertainty",
             vjust = 1, hjust = 0, size=4) +

    # scales
    scale_y_continuous("Nr of specimens", breaks=seq(0, 9*10^7, 2*10^7), labels=c("0M", "20M", "40M", "60M", "80M"),
                       limits = c(-16000000, 9*10^7), expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    theme(axis.text = element_text(size=11),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(colour="grey70", linetype = "dashed", size = 0.2),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("Kingdom")
p.bar.kingdom
f.out <- paste0(manuscript.dir, "figures/figure_4a.png")
ggsave(plot=p.bar.kingdom, device="png", filename = f.out, width = 30, height = 10, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig 4b
# Adapted from: https://www.r-graph-gallery.com/299-circular-stacked-barplot.html
df.kingdom.phylum <- results[['taxonomy.numbers']]$sci.name.non.geo %>%
    group_by(kingdom, phylum) %>%
    summarise(non.geo=sum(n), .groups = "drop") %>%
    inner_join(
        results[['taxonomy.numbers']]$family.georeferenced.all %>%
            group_by(kingdom, phylum) %>%
            summarise(geo=sum(georef), unc=sum(unc), .groups = "drop"), by=c("kingdom", "phylum")) %>%
    mutate(total = non.geo + geo,
           perc.geo = round(geo / total * 100, 1),
           perc.unc = round(unc / total * 100, 1)) %>%
    filter(kingdom != 'incertae sedis' & phylum != 'UNKNOWN') %>%
    arrange(phylum, total)

df.phylum <- df.kingdom.phylum %>%
    dplyr::select(-perc.geo, -perc.unc, -non.geo) %>%
    gather(key = "var", value="value", -c(1,2)) %>%
    ungroup() %>%
    filter(kingdom != 'incertae sedis' & phylum != 'UNKNOWN') %>%
    mutate(kingdom = as.factor(kingdom))

df.phylum.total <- df.phylum %>% filter(var == 'total')
df.phylum.geo <- df.phylum %>% filter(var == 'geo')
df.phylum.unc <- df.phylum %>% filter(var == 'unc')

empty_bar <- 4
nObsType <- nlevels(as.factor(df.phylum.total$var))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(df.phylum.total$kingdom)*nObsType, ncol(df.phylum.total)) )
colnames(to_add) <- colnames(df.phylum.total)
to_add$kingdom <- rep(levels(df.phylum.total$kingdom), each=empty_bar*nObsType )
df.phylum.total <- rbind(df.phylum.total, to_add)
df.phylum.total <- df.phylum.total %>% arrange(kingdom, phylum)
df.phylum.total$id <- rep( seq(1, nrow(df.phylum.total)/nObsType) , each=nObsType)

# Need to add id to geo and unc df in order to correctly place bars
df.phylum.geo <- df.phylum.geo %>%
    inner_join(df.phylum.total %>% dplyr::select(id, kingdom, phylum), by=c("kingdom", "phylum"))
df.phylum.unc <- df.phylum.unc %>%
    inner_join(df.phylum.total %>% dplyr::select(id, kingdom, phylum), by=c("kingdom", "phylum"))

# Get the name and the y position of each label
label_data <- df.phylum.total %>% rename("tot"="value")
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) / number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- df.phylum.total %>%
    group_by(kingdom) %>%
    summarize(start=min(id), end=max(id) - empty_bar, .groups = "drop") %>%
    rowwise() %>%
    mutate(title=mean(c(start, end))) %>%
    mutate(end = ifelse(end - start == 0, end + 0.5, end)) # Not perfect, but to avoid a 0

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

p.phylum.circular.bar <- ggplot() +

    # Concentric circles for Y axis shown labels
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0),
                 colour = "grey90", linetype = "solid", alpha=1, size=0.3, inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = 0, y = 0.5*10^7, xend = 105, yend = 0.5*10^7),
                 colour = "grey90", linetype = "solid", alpha=0.5, size=0.3, inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = 0, y = 2.5*10^7, xend = 105, yend = 2.5*10^7),
                 colour = "grey90", linetype = "solid", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = 0, y = 5.0*10^7, xend = 105, yend = 5.0*10^7),
                 colour = "grey90", linetype = "solid", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = 0, y = 7.5*10^7, xend = 105.5, yend = 7.5*10^7),
                 colour = "grey90", linetype = "solid", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = 0, y = 10^8, xend = 105.5, yend = 10^8),
                 colour = "grey90", linetype = "solid", alpha=1, size=0.3, inherit.aes = FALSE ) +

    # Y axis labels
    annotate("text", x = rep(107.5, 6), y = c(0, 0.5*10^7, 2.5*10^7, 5.0*10^7, 7.5*10^7, 10^8),
             label = c("0", "5M", "25M", "50M", "75M", "100M"),
             color="grey30", size=5 , angle=0, hjust=1)

# The actual data
p.phylum.circular.bar <- p.phylum.circular.bar +
    geom_bar(data=df.phylum.total, aes(x=id, y=value, fill=var),
             stat="identity", fill = gf.col["NC"], alpha=0.5) +
    geom_bar(data=df.phylum.geo, aes(x=id, y=value, fill=var),
             stat="identity", fill = gf.col["CO"], alpha=0.5) +
    geom_bar(data=df.phylum.unc, aes(x=id, y=value, fill=var),
             stat="identity", fill = gf.col["CU"], alpha=1)

# Radial dotted lines
p.phylum.circular.bar <- p.phylum.circular.bar +
    geom_segment(data=df.phylum.total %>% filter(!is.na(value)),
                 aes(x = id, y = 0.2*10^7, xend = id, yend = 4.8*10^7),
                 colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE, linetype='dotted')

# Radial phylum labels
p.phylum.circular.bar <- p.phylum.circular.bar +
    geom_text(data=label_data, aes(x=id, y= ifelse(tot + 52000000 > 52000000, 52000000, tot + 52000000),
                                   label=phylum, hjust=hjust),
              color="black", fontface="plain", alpha=0.6, size=6,
              angle= label_data$angle, inherit.aes = FALSE )

# Add base line information
p.phylum.circular.bar <- p.phylum.circular.bar +
    geom_segment(data=base_data, aes(x = 0, y = -1000000, xend = 105, yend = -1000000),
                 colour = "grey", alpha=0.8, size=0.4 , inherit.aes = FALSE )

kingdom.text.colour <- "darkblue"
p.phylum.circular.bar <- p.phylum.circular.bar +
    annotate("text", x = 15, y = 95000000,
             label = "Animalia",
             color=kingdom.text.colour, size=7 , angle=303, hjust=0, vjust=0) +
    annotate("text", x = 35, y = 95000000,
             label = "Archaea",
             color=kingdom.text.colour, size=7 , angle=239, hjust=0, vjust=0) +
    annotate("text", x = 45, y = 95000000,
             label = "Bacteria",
             color=kingdom.text.colour, size=7 , angle=207, hjust=0, vjust=0) +
    annotate("text", x = 64, y = 95000000,
             label = "Chromista",
             color=kingdom.text.colour, size=7 , angle=139, hjust=0, vjust=0) +
    annotate("text", x = 77, y = 95000000,
             label = "Fungi",
             color=kingdom.text.colour, size=7 , angle=101, hjust=0, vjust=0) +
    annotate("text", x = 87, y = 95000000,
             label = "Plantae",
             color=kingdom.text.colour, size=7 , angle=65, hjust=0, vjust=0) +
    annotate("text", x = 98, y = 95000000,
             label = "Protozoa",
             color=kingdom.text.colour, size=7 , angle=28, hjust=0, vjust=0)

p.phylum.circular.bar <- p.phylum.circular.bar +
    ylim(-20000000, 100000000) +
    coord_polar() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")
    )
p.phylum.circular.bar
f.out <- paste0(manuscript.dir, "figures/figure_4b.png")
ggsave(plot=p.phylum.circular.bar, device="png", filename = f.out, width = 35, height = 35, units = "cm", dpi = 300)
cat("  ==>", f.out, "\n")

#### FIGURE 5 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure 5 ...\n")

df.year.geo <- results$event.numbers$event.year.non.georeferenced %>%
    mutate(event.year = as.numeric(event.year)) %>%
    rename("n.non.geo"="n") %>%
    left_join(results$event.numbers$event.year.georeferenced.all,
              by = c("country.code", "kingdom", "phylum", "class", "event.year")) %>%
    mutate(n.non.geo = ifelse(is.na(n.non.geo), 0, n.non.geo),
           n.geo = ifelse(is.na(n.geo), 0, n.geo),
           n.geo.unc = ifelse(is.na(n.geo.unc), 0, n.geo.unc))

year.interval <- 10
df.year <- df.year.geo %>%
    group_by(event.year) %>%
    summarise(n.non.geo = sum(n.non.geo), n.geo = sum(n.geo), n.geo.unc = sum(n.geo.unc), .groups = "drop") %>%
    filter(!is.na(event.year)) %>%
    mutate(ey = floor(event.year/year.interval) * year.interval + year.interval / 2)

# Partial results: Prior to XIXth century (1800)
# NC records
year <- 1800
dfng <- results$event.numbers$event.year.non.georeferenced %>%
    mutate(event.year = as.integer(event.year)) %>%
    filter(!is.na(event.year) & event.year < year)

# CO and CU records
dfg <- results$event.numbers$event.year.georeferenced.all %>%
    filter(!is.na(event.year) & event.year < year) %>%
    mutate(n.geo = ifelse(is.na(n.geo), 0, n.geo), n.geo.unc = ifelse(is.na(n.geo.unc), 0, n.geo.unc)) %>%
    mutate(n = n.geo + n.geo.unc)

writeResult("figure-5a-number_ps_before_year_1800", sum(dfng$n) + sum(dfg$n))
writeResult("figure-5a-number_nc_before_year_1800", sum(dfng$n))
writeResult("figure-5a-number_c_before_year_1800", sum(dfg$n))
writeResult("figure-5a-number_co_before_year_1800", sum(dfg$n.geo))
writeResult("figure-5a-number_cu_before_year_1800", sum(dfg$n.geo.unc))

# Fig 5a
colors = c("total"="dodgerblue", "geo"="dodgerblue4", "unc"="brown3")
leg.x <- 1600
leg.y <- 4000000
p.year <- ggplot(df.year) +
    geom_bar(aes(x = ey, y=n.non.geo), stat = "identity", fill = gf.col["NC"]) +
    geom_bar(aes(x = ey, y=n.geo), stat = "identity", fill = gf.col["CO"]) +
    geom_bar(aes(x = ey, y=n.geo.unc), stat = "identity", fill = gf.col["CU"]) +
    scale_y_continuous("Nr of specimens", breaks=seq(0, 5*10^6, 1*10^6), labels=c("0", "1M", "2M", "3M", "4M", "5M"),
                       limits = c(0, 5*10^6)) +
    scale_x_continuous("Year of collection event", limits=c(1600, 2025)) +
    # white background
    geom_rect(aes(xmin=leg.x, ymin=leg.y - 300000, xmax = leg.x + 160, ymax = leg.y + 1*10^6), fill="white") +
    # Total
    geom_rect(aes(xmin=leg.x + 5, ymin=leg.y + 800000, xmax=leg.x + 20, ymax=leg.y + 500000), fill=gf.col["NC"]) +
    annotate("text", x=leg.x + 25, y=leg.y + 650000, label="Preserved specimens",
             hjust = 0, size=5) +
    # Georef
    geom_rect(aes(xmin=leg.x + 5, ymin=leg.y + 450000, xmax=leg.x + 20, ymax=leg.y + 150000), fill=gf.col["CO"]) +
    annotate("text", x=leg.x + 25, y=leg.y + 300000, label="Records with coordinates",
             hjust = 0, size=5) +
    # Uncertainty
    geom_rect(aes(xmin=leg.x + 5, ymin=leg.y + 100000, xmax=leg.x + 20, ymax=leg.y - 200000), fill=gf.col["CU"]) +
    annotate("text", x=leg.x + 25, y=leg.y - 50000, label="Records with uncertainty",
             hjust = 0, size=5) +

    theme(axis.text = element_text(size=11),
          panel.grid.major.x = element_line(colour="grey70", linetype = "dashed", size = 0.2),
          panel.grid.major.y = element_line(colour="grey70", linetype = "dashed", size = 0.2),
          axis.title = element_text(size = 12),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank())
p.year
f.out <- paste0(manuscript.dir, "figures/figure_5a.png")
ggsave(plot=p.year, device="png", filename = f.out, width = 20, height = 10, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig 5b
f.in <- paste0(data.processed.dir, "gbif/", date.string, "_uncertainty_data.csv")
df.unc <- vroom(f.in, delim=',', col_types = cols())
df.unc.2d <- df.unc %>%
    dplyr::select(kingdom, event.year=year, unc) %>%
    filter(!is.na(event.year) & !is.na(unc))

bins <- round((max(df.unc.2d$event.year) - min(df.unc.2d$event.year)) / 5,0)
p.stat.unc <- ggplot(df.unc.2d) +
    geom_bin2d(aes(unc, event.year), position="identity", bins=bins) +
    scale_x_log10("Coordinate uncertainty in meters",
                  breaks=c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
                  labels=c("0.001", "0.01", "0.1", "1", "10^1", "10^2", "10^3", "10^4", "10^5", "10^6", "10^7")) +
    scale_y_continuous("Year collected", breaks = c(seq(1600, 2000, 50)), limits=c(1600, 2020)) +
    scale_fill_continuous("N", trans = "log", type = "viridis",
                          breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000),
                          labels=c("1e+00", "1e+01", "1e+02", "1e+03", "1e+04", "1e+05", "1e+06")) +
    theme(panel.grid.major = element_line(colour="grey70", linetype = "dashed", size = 0.1),
          axis.text = element_text(hjust=1, size=11),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "N",
                                  title.position = "bottom",
                                  label.theme = element_text(hjust=0.5, size=6),
                                  title.theme = element_text(size=7, hjust = 0.5),
                                  barwidth = unit(8, "cm"),
                                  barheight = unit(0.2, "cm")))
p.stat.unc
f.out <- paste0(manuscript.dir, "figures/figure_5b.png")
ggsave(plot=p.stat.unc, device="png", filename = f.out, width = 20, height = 20, units = "cm", dpi = 300)
cat("  ==>", f.out, "\n")


#### FIGURE 6 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure 6 ...\n")
peaks <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 71, 100, 150, 200, 250, 300, 402, 500, 707,
           1000, 1415, 2000, 2500, 3036, 3535, 5000, 7071, 10000, 14143, 20000, 25000,
           30000, 50000, 100000, 999999)
# Fig 6a
p.unc.density <- ggplot(df.unc) + geom_density(aes(x=round(unc))) +
    scale_x_continuous('', trans = 'log', limits = c(0.5, 1200000), expand = c(0, 0)) +
    scale_y_continuous('', limits = c(-0.35, 0.9), expand = c(-0.02, 0), breaks = seq(0, 0.75, 0.25)) +
    geom_vline(xintercept=peaks, colour = gf.col["CU"], linetype = 'dashed', size = 0.2, alpha = 0.4) +
    geom_vline(xintercept=0, colour = "grey70", linetype = 'solid', size = 0.2) +
    annotate('rect', xmin = 0.5, xmax = 1150000, ymin = -0.35, ymax = -0.003, fill = 'white') +
    annotate('text', x = 700, y = -0.25, label = 'Uncertainty (m)', size = 5) +
    annotate(geom = "text",
             x = peaks, y = -0.01,
             label = formatC(peaks, format='f', digits = 0),
             color = gf.col["CU"], angle = 90, hjust = 1, size = 3.25) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.background = element_blank())
p.unc.density

f.out <- paste0(manuscript.dir, "figures/figure_6a.png")
# f.out <- paste0(tmp.dir, "figure_6a.png")
ggsave(plot=p.unc.density, device="png", filename = f.out, width = 29, height = 8, units = "cm", dpi = 300)
cat("  ==>", f.out, "\n")

writeResult("figure-6a-percentage-uncertainty-density-peaks", round((df.unc %>% filter(round(unc) %in% peaks) %>% nrow())/nrow(df.unc) * 100,1))

# Fig 6b
# Unique locations with same uncertainty per species
df1 <- df.unc %>%
    dplyr::select(taxon, lon, lat, unc, unc.cat) %>%
    group_by(taxon, lon, lat, unc, unc.cat) %>%
    summarise(n=n(), .groups = "drop")

# Number of records available for research at each resolution
df2 <- df1 %>%
    group_by(unc.cat) %>%
    summarise(n=n()) %>%
    mutate(unc.cat = factor(unc.cat,
                            levels = c("<= 1m", "10m", "100m", "250m", "1km", "5km", "10km", "50km", "100km", "> 100km"))) %>%
    arrange(unc.cat) %>%
    mutate(cs = cumsum(n), cum.perc=round(cs/sum(n)*100,1), .groups = 'drop')

df3 <- bind_rows(df2 %>% dplyr::select(-.groups) %>% mutate(subset = "Unique taxon locations"),
                 results$world.numbers$georeferenced.unc.cat %>% na.omit() %>%
                     dplyr::select(-.groups) %>% mutate(subset = "All records"))

p.bar.resolution.research <- ggplot() +
    geom_col(data=df3, aes(x = unc.cat, y = cs, fill = subset), width= 0.75,
             position = position_dodge2(width=0.9, reverse = T)) +
    geom_text(data=df3, aes(x = unc.cat, y = cs, colour = subset,
                            label = paste0(format(round(cs/1000000, 1), digits = 1,
                                                  scientific=F, big.mark = " "), paste0(", ", cum.perc, " %")),
                            angle = 90, size = 12),
              position = position_dodge2(width=0.7, reverse = T), hjust = -0.05, vjust = 0.5, size = 3.5) +
    scale_y_continuous("Nr of specimens (in millions of records)",
                       limits=c(0, 43000000),
                       breaks=seq(0, 43000000, 5000000),
                       labels = paste0(format(seq(0, 43000000, 5000000)/1000000,
                                              scientific = F, big.mark = " "), "M"),
                       expand = c(0,0)) +
    scale_x_discrete("Target research resolution") +
    scale_fill_manual(values = c("grey40", "firebrick")) +
    scale_colour_manual(values = c("grey40", "firebrick")) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 13),
          axis.title.x = element_text(margin = margin(t = -7, r = 0, b = 0, l = 0)),
          panel.grid.major = element_blank(),
          axis.line = element_line(size = 0.3, colour = "grey30"),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = c(0.15, 0.8),
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
p.bar.resolution.research
f.out <- paste0(manuscript.dir, "figures/figure_6b.png")
ggsave(plot=p.bar.resolution.research, device="png", filename = f.out,
       width = 29, height = 12, units = "cm", dpi = 300)
cat("  ==>", f.out, "\n")

# Fig 6c
df.unc.8857 <- df.unc %>% dplyr::select(lon, lat, kingdom, unc)
df.unc.8857 <- df.unc.8857 %>% bind_cols(coordTrans(df.unc.8857 %>% dplyr::select(lon, lat), 4326, 8857))
br <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
grid.sizes <- c(1000, 5000, 10000, 50000, 100000, 500000)
df.res.kgd <- tibble()
for(i in 1:length(grid.sizes)){
    gs <- grid.sizes[i]
    cat("  preparing data for grid size", gs, "...\n")
    df <- df.unc.8857 %>%
        filter(unc <= gs) %>%
        mutate(xc = floor(x/gs) * gs + gs/2,
               yc = floor(y/gs) * gs + gs/2) %>%
        group_by(kingdom, xc, yc) %>%
        count() %>%
        mutate(grid.size = gs)
    df.res.kgd <- bind_rows(df.res.kgd, df)
}
wc <- st_transform(world_coastline, 8857)
for(i in 1:length(grid.sizes)){
    gs <- grid.sizes[i]
    cat("  plotting grid size", gs, "...\n")
    dfg <- df.res.kgd %>% filter(grid.size == gs) %>% group_by(xc,yc) %>% summarise(n = sum(n), .groups = "drop")
    p <- ggplot(dfg) +
        geom_tile(aes(x=xc, y=yc, fill =n)) +
        geom_sf(data = wc, fill = NA, colour = 'white', size = 0.1) +
        scale_fill_gradient("N", trans = "log", low = "yellow", high = "red",
                            breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
                            limits = c(1, 10000000)) +
        theme(panel.background = element_rect(fill = "grey20"),
              axis.title = element_blank(),
              legend.title = element_blank(),
              legend.text = element_blank(),
              legend.position = "none",
              panel.grid = element_blank(),
              legend.background = element_rect(fill="transparent")) +
        guides(fill=guide_colourbar(title.vjust = 2, title.hjust = 0.7,
                                    barwidth = unit(5, "cm"), barheight = unit(0.1, "cm"),
                                    direction = "horizontal"))
    print(p)
    f.out <- paste0(manuscript.dir, "figures/figure_6c_", formatC(gs, digits=0, format="f"), ".png")
    ggsave(plot=p, device="png", filename = f.out, width = 9, height = 4.5, units = "cm", dpi = 600)
    cat("    ==>", f.out, "\n")
}

# Get legend, common to all, and save to compose figure externally with a graphics program
f.out <- paste0(manuscript.dir, "figures/figure_6_legend.png")
ggsave(plot=get_legend(p + theme(legend.position = c(0.5, 0.001))), device="png", filename = f.out, width = 6, height = 0.3, units = "cm", dpi = 600)
cat("    ==>", f.out, "\n")

#### FIGURE 7 --------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure 7 ...\n")

getFigureFilePath <- function(model.taxon, figure.number = "", figure.letter = "", subfig = "", suffix = ""){
    path <- paste0(manuscript.dir, "figures/")
    path <- paste0(path, "figure_", figure.number, "_", figure.letter, subfig, "_",
                   model.taxon$taxon.abbrv, "_", model.taxon$max.uncertainty, suffix, ".png")
    return(path)
}
getModellingDir <- function(model.taxon){
    return(paste0(model.taxon$modelling.dir, model.taxon$taxon.abbrv, "-uncertainty/"))
}
getQuantileModel <- function(range.results.df, quant){
    idx.q <- as.integer(nrow(range.results.df) * quant)
    q <- range.results.df %>% arrange(range.area) %>% slice(idx.q)
    return(q)
}
getPredictionRaster <- function(model.taxon, id){
    f <- paste0(getModellingDir(model.taxon), "predictions/maxunc_", model.taxon$max.uncertainty, "/m_", id, "_final_prediction.tif")
    r <- raster(f)
    return(r)
}
getPredictionDF <- function(model.taxon, id, threshold){
    r <- getPredictionRaster(model.taxon, id)
    r.df <- data.frame(as(r, "SpatialPixelsDataFrame")) %>%
        setNames(., c("v", "lon", "lat")) %>%
        mutate(v = ifelse(v <= threshold, NA, v))
    return(r.df)
}
plotSuitabilityRangeMap <- function(model.taxon, df, figure.number, subfig, figure.letter, only.legend = F){
    p <- ggplot(df) +
        geom_sf(data=world_coastline, size = 0.01, alpha = 1) +
        geom_raster(aes(x = lon, y = lat, fill = v)) +
        scale_fill_viridis(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
        scale_x_continuous("", limits = c(model.taxon$region.ext@xmin, model.taxon$region.ext@xmax)) +
        scale_y_continuous("", limits = c(model.taxon$region.ext@ymin, model.taxon$region.ext@ymax)) +
        coord_sf() +
        theme(panel.grid = element_blank(), panel.background = element_blank(),
              legend.position = "bottom",
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 3)) +
        guides(fill = guide_colourbar(title = "Suitability index",
                                      barheight = 0.30, barwidth=9, direction = "horizontal",
                                      title.position="top", title.hjust = 0.5,title.vjust = 0.5, nbin = 50,
                                      title.theme = element_text(size = 8),
                                      label.theme = element_text(size = 6)))
    if(only.legend){
        leg <- get_legend(p)
        f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig, "_legend")
        ggsave(plot=leg, device="png", filename = f.out, width = 6, height = 1.2, units = "cm", dpi = 600)
        cat("    ==>", f.out, "\n")
    } else {
        f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig)
        ggsave(plot=p + theme(legend.position = "none"), device="png", filename = f.out, width = 6, height = 5, units = "cm", dpi = 600)
        cat("    ==>", f.out, "\n")
    }
}
plotRangeDiffMap <- function(model.taxon, q.05, q.95, figure.number, subfig, figure.letter){
    r.q05.pr <- getPredictionRaster(model.taxon, q.05$model.id)
    r.q05.pr[r.q05.pr < q.05$maxSSS] <- 0
    r.q05.pr[r.q05.pr >= q.05$maxSSS] <- 1

    r.q95.pr <- getPredictionRaster(model.taxon, q.95$model.id)
    r.q95.pr[r.q95.pr < q.95$maxSSS] <- 0
    r.q95.pr[r.q95.pr >= q.95$maxSSS] <- 1

    r.diff.pr <- abs(r.q95.pr - r.q05.pr)
    r.q05.area <- sum(na.omit(getValues(r.q05.pr)))
    r.q95.area <- sum(na.omit(getValues(r.q95.pr)))
    times.q05 <- r.q95.area / r.q05.area
    writeResult(paste0("figure-", model.taxon$taxon.abbrv, "-times_q95_of_q05_ranges"), times.q05)

    range.diff <- as(r.diff.pr, "SpatialPixelsDataFrame") %>%
        as.data.frame() %>%
        setNames(., c("v", "lon", "lat")) %>%
        mutate(v = factor(v, levels = c("0", "1")))

    p <- ggplot(range.diff) +
        geom_sf(data=world_coastline, size = 0.01, colour = "grey70") +
        geom_raster(aes(x = lon, y = lat, fill = v)) +
        scale_fill_manual(values = c("0" = "grey70", "1" = "firebrick3")) +
        coord_sf() +
        scale_x_continuous("", limits = c(model.taxon$region.ext@xmin, model.taxon$region.ext@xmax)) +
        scale_y_continuous("", limits = c(model.taxon$region.ext@ymin, model.taxon$region.ext@ymax)) +
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              legend.position = "none",
              axis.text = element_blank(),
              axis.ticks = element_blank())
    # p
    f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig)
    ggsave(plot=p, device="png", filename = f.out, width = 6, height = 5, units = "cm", dpi = 600)
    cat("    ==>", f.out, "\n")
}
plotRangeAreaDensity <- function(model.taxon, range.results.df, q.05, q.95, figure.number, subfig, figure.letter){
    breaks <- c(min(range.results.df$range.area),q.05$range.area, q.95$range.area, max(range.results.df$range.area))
    breaks.labels <- round(breaks/1000)

    p <- ggplot(range.results.df, aes(x=range.area, y = ..scaled..)) +
        geom_density() +
        geom_area(
            aes(x = stage(range.area, after_stat = oob_censor(x, c(q.05$range.area, q.95$range.area)))),
            fill = "grey60",
            stat = "density"
        ) +
        geom_area(
            aes(x = stage(range.area, after_stat = oob_censor(x, c(0, q.05$range.area)))),
            fill = "grey80",
            stat = "density") +
        geom_area(
            aes(x = stage(range.area, after_stat = oob_censor(x, c(q.95$range.area, Inf)))),
            fill = "grey80",
            stat = "density") +
        scale_x_continuous(breaks = breaks, labels = breaks.labels, expand = c(0.015, 0.015)) +
        scale_y_continuous(expand = c(0, 0.02)) +
        annotate("segment", x = breaks[2], y = 0.6, xend = breaks[3], yend = 0.6,
                 colour = "red", size = 0.4, alpha = 0.5) +
        annotate("text", x = breaks[3], y = 0.66, colour = "red", size = 2.5, hjust = 1,
                 label = paste0("x", formatC(breaks[3]/breaks[2], format = "f", digits = 1))) +
        annotate('segment', breaks[1], y = 0.8, xend = breaks[4], yend = 0.8,
                 size = 0.4, colour = "blue",alpha = 0.5) +
        annotate("text", x = breaks[4], y = 0.86, colour = "blue", size = 2.5, hjust = 1,
                 label = paste0("x", formatC(breaks[4]/breaks[1], format = "f", digits = 1))) +
        theme(axis.text.x = element_text(size = 6, angle = 45), axis.text.y = element_blank(),
              axis.ticks.x = element_line(colour = "grey60", size = 0.3),
              axis.ticks.y = element_blank(),
              axis.title = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())
    # p
    f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig)
    ggsave(plot=p, device="png", filename = f.out, width = 6, height = 4, units = "cm", dpi = 300)
    cat("    ==>", f.out, "\n")
}
plotLegend <- function(p, figure.number, figure.letter){
    leg <- get_legend(p)
    f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig = "", suffix = "_legend")
    ggsave(plot=leg, device="png", filename = f.out, width = 6, height = 1.2, units = "cm", dpi = 600)
    cat("  ==>", f.out, "\n")
}
# Plot maps of minimum and maximum estimated distributions
plotOccurrencesUncertainty <- function(model.taxon, figure.number, figure.letter, only.legend = F){
    # Prepare occurrences data
    df.occ <- vroom::vroom(paste0(output.dir, "maxent/", model.taxon$taxon.abbrv, "-uncertainty",
                                  "/data/",
                                  model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_gbif_data.csv"),
                           col_types = cols()) %>%
        dplyr::select(lon, lat, unc) %>%
        mutate(unc.cat = cut(unc, breaks=c(0, 1000, 5000, 10000, 50000, 100000, 1000000),
                             labels = c("<=1k", "(1k-5k]", "(5k-10k]", "(10k-50k]", "(50k-100k]", ">100k")))

    # Output number of occurrences to results text file
    occ.with.unc.nr <- df.occ %>% filter(!is.na(unc)) %>% count() %>% pull(n)
    writeResult(paste0("figure-", model.taxon$taxon.abbrv, "-number_of_occurrences_total"), nrow(df.occ))
    writeResult(paste0("figure-", model.taxon$taxon.abbrv, "-number_of_occurrences_CU"), occ.with.unc.nr)

    # Plot occurrences map
    breaks <- c(1000, 5000, 10000, 50000, 100000, 500000, 1000000)
    labels <- paste0(breaks/1000, " km")

    df.occ.t <- coordTrans(df.occ %>% dplyr::select(lon, lat), 4326, 8857) %>% setNames(., c("lon", "lat")) %>%
        bind_cols(df.occ %>% dplyr::select(unc, unc.cat))

    p <- ggplot() +
        geom_sf(data=world_coastline, fill="red", size = 0.3) +
        geom_point(data = df.occ %>% filter(is.na(unc)), aes(x = lon, y = lat),
                   shape='+', size=1, colour = "grey40") +
        geom_point(data = df.occ %>% filter(!is.na(unc)) %>% arrange(desc(unc)),
                   aes(x = lon, y = lat, size = unc, fill=unc.cat), shape=21, colour="grey20") +
        scale_x_continuous("", limits = c(model.taxon$region.ext@xmin, model.taxon$region.ext@xmax)) +
        scale_y_continuous("", limits = c(model.taxon$region.ext@ymin, model.taxon$region.ext@ymax)) +
        scale_size_continuous("Uncertainty (m)", limits = c(0, 1000000), breaks = breaks,
                              labels = labels,
                              guide = guide_legend(title.position = "top", nrow = 1, order = 1,
                                                   title.vjust = 5.0,
                                                   keywidth = unit(0.4, "cm"), keyheight = unit(0.4, "cm"))) +
        scale_fill_discrete("",
                            guide = guide_legend(nrow = 1, order = 2,
                                                 keywidth = unit(0.4, "cm"), keyheight = unit(0.4, "cm"),
                                                 override.aes = list(size = 3))) +
        coord_sf() +
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              legend.position = "bottom",
              legend.box = "vertical",
              legend.direction = "horizontal",
              panel.grid.major = element_line(linetype = "dashed", colour = "grey70", size = 0.2),
              axis.text = element_text(size = 6),
              axis.ticks = element_blank(),
              legend.title = element_text(size = 7),
              legend.text = element_text(size = 6),
              legend.key = element_rect(fill = NA),
              legend.spacing.y = unit(-0.1, "cm"))
    # Save the legend, just once
    if(!only.legend){
        # print(p + theme(legend.position = "none"))
        f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig = "")
        ggsave(plot=p + theme(legend.position = "none"), device="png", filename = f.out, width = 9, height = 7.2, units = "cm", dpi = 600)
        cat("    ==>", f.out, "\n")
    } else {
        leg <- get_legend(p)
        f.out <- getFigureFilePath(model.taxon, figure.number, figure.letter, subfig = "", suffix = "_occurrences_legend")
        ggsave(plot=leg, device="png", filename = f.out, width = 12, height = 2, units = "cm", dpi = 600)
        cat("    ==>", f.out, "\n")
    }

}
saveFigureSummary <- function(model.taxon, range.results.df, q.05, q.95, figure.number){
    min.df <- range.results.df[which.min(range.results.df$range.area),]
    max.df <- range.results.df[which.max(range.results.df$range.area),]
    t <- tibble(
        timestamp = as.POSIXlt(Sys.time()),
        taxon = model.taxon$taxon.abbrv,
        region = model.taxon$region.name,
        models.nr = model.taxon$models.nr,
        background.nr = model.taxon$n.background,
        max.uncertainty = model.taxon$max.uncertainty,
        q = c("0.05", "0.95", "min", "max"),
        q.id = c(q.05$model.id, q.95$model.id, min.df$model.id, max.df$model.id),
        q.maxSSS = c(q.05$maxSSS, q.95$maxSSS, min.df$maxSSS, max.df$maxSSS),
        q.suitability = c(q.05$suitability, q.95$suitability, min.df$suitability, max.df$suitability),
        q.range.area = c(q.05$range.area, q.95$range.area, min.df$range.area, max.df$range.area)
    )
    f.out <- paste0(manuscript.dir, "figures/figure_", figure.number, "_summary.csv")
    if(!file.exists(f.out)){
        write_csv(t, f.out, append = F)
    } else {
        write_csv(t, f.out, append = T)
    }
    cat("    ==>", f.out, "\n")
}
generateTaxonFigures <- function(taxon, max.unc, figure.number){
    cat("Generating figure for", taxon, "and uncertainty threshold", max.unc, "...\n")
    i <- list.findi(taxons.to.model, taxon.abbrv == taxon & max.uncertainty == max.unc)
    model.taxon <- taxons.to.model[[i]]
    figure.letter <- switch (taxon,
                             "rhododendron_groenlandicum" = "a",
                             "guazuma_ulmifolia" = "b",
                             "eucalyptus_gongylocarpa" = "c"
    )
    f.range.results <- paste0(getModellingDir(model.taxon), model.taxon$taxon.abbrv,
                              "_maxunc_", model.taxon$max.uncertainty, "_suitability_range_results.csv")
    range.results.df <- vroom(f.range.results, col_types = cols())

    # Plot occurrences subfigure
    if(model.taxon$max.uncertainty == "Inf"){
        cat("  creating occurrences subfigure ...\n")
        plotOccurrencesUncertainty(model.taxon, figure.number, figure.letter, F)
        plotOccurrencesUncertainty(model.taxon, figure.number, figure.letter, T)
    }
    # Plot quantile 0.05 distribution map subfigure
    cat("  creating quantile 0.05 distribution map subfigure ...\n")
    q.05 <- getQuantileModel(range.results.df, 0.05)
    df <- getPredictionDF(model.taxon, q.05$model.id, q.05$maxSSS)
    plotSuitabilityRangeMap(model.taxon, df, figure.number, 1, figure.letter, F)

    # Plot quantile 0.95 distribution map subfigure
    cat("  creating quantile 0.95 distribution map subfigure ...\n")
    q.95 <- getQuantileModel(range.results.df, 0.95)
    df <- getPredictionDF(model.taxon, q.95$model.id, q.95$maxSSS)
    plotSuitabilityRangeMap(model.taxon, df, figure.number, 2, figure.letter, F)

    # Plot the distribution map legend only, just once is needed
    if(i == 1){
        plotSuitabilityRangeMap(model.taxon, df, figure.number, 1, figure.letter, T)
    }

    # Plot range difference map between quantiles 0.05 and 0.95
    cat("  creating range difference map subfigure ...\n")
    plotRangeDiffMap(model.taxon, q.05, q.95, figure.number, 3, figure.letter)

    # Plot range area density showing quantiles cut-offs.
    cat("  creating range area density subfigure ...\n")
    plotRangeAreaDensity(model.taxon, range.results.df, q.05, q.95, figure.number, 4, figure.letter)

    # Figure summary
    cat("  saving summary ...\n")
    saveFigureSummary(model.taxon, range.results.df, q.05, q.95, figure.number)
}

unc.thresholds <- c(Inf)
for(i in 1:length(unc.thresholds)){
    generateTaxonFigures("eucalyptus_gongylocarpa", unc.thresholds[i], "7")
    generateTaxonFigures("rhododendron_groenlandicum", unc.thresholds[i], "7")
    generateTaxonFigures("guazuma_ulmifolia", unc.thresholds[i], "7")
}

# Convenient output for data to be placed manually on publishable figure
df <- vroom(paste0(manuscript.dir, "figures/figure_7_summary.csv"), col_types = cols(), delim = ",")
for(taxon in c("rhododendron_groenlandicum", "guazuma_ulmifolia", "eucalyptus_gongylocarpa")){
    cat("=====", taxon, "=====\n")
    df.rg <- df %>% filter(taxon == !!taxon & max.uncertainty == Inf)
    cat("  Quantiles 0.05 and 0.95\n")
    q05.area <- df.rg[df.rg$q == "0.05", "q.range.area"] %>% as.numeric()
    q95.area <- df.rg[df.rg$q == "0.95", "q.range.area"] %>% as.numeric()
    q.dif <- q95.area - q05.area

    min.area <- df.rg[df.rg$q == "min", "q.range.area"] %>% as.numeric()
    max.area <- df.rg[df.rg$q == "max", "q.range.area"] %>% as.numeric()
    minmax.dif <- max.area - min.area

    cat("    q 0.05:", q05.area, "- q 0.95:", q95.area, "\n")
    cat("    difference", q.dif, "\n")
    cat("    % difference:", round(q.dif/q05.area * 100, 1), "%\n")

    cat("  Minimum and maximum\n")
    cat("    min:", min.area, "- max:", max.area, "\n")
    cat("    difference", round(minmax.dif/min.area * 100, 1), "%\n")
}

#### SUPPLEMENTARY MATERIALS ================================================================== ####
#### TABLE ST1 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Gathering VIF analysis results ...\n")
df.res <- data.frame()
for(i in 1:length(taxons.to.model)){
    model.taxon <- taxons.to.model[[i]]
    modelling.dir <- paste0(output.dir, "maxent/", model.taxon$taxon.abbrv, "-uncertainty/")
    f <- paste0(modelling.dir, model.taxon$taxon.abbrv, "_vif_analysis.rds")
    vif.df <- readRDS(f)@results %>% arrange(desc(VIF)) %>%
        mutate(taxon = model.taxon$taxon.abbrv)
    df.res <- bind_rows(df.res, vif.df)
}

f.out <- paste0(manuscript.dir, "tables/table_ST1.csv")
write_csv(df.res %>% dplyr::select(taxon, variable=Variables, vif=VIF), f.out)
cat("  ==>", f.out, "\n")

cat("\n", rep("-", 80), "\n", sep = "")
cat("Gathering model performance results ...\n")
df.res <- data.frame()
unc.thresholds <- c(3536, 7071, 14142, 28284, Inf)
taxon.models <- list.filter(taxons.to.model, max.uncertainty %in% unc.thresholds)
for(model.taxon in taxon.models){
    modelling.dir <- paste0(output.dir, "maxent/", model.taxon$taxon.abbrv, "-uncertainty/")
    f <- paste0(modelling.dir, model.taxon$taxon.abbrv, "_maxunc_", model.taxon$max.uncertainty, "_maxent_results.csv")
    df <- vroom::vroom(f, col_types = cols()) %>%
        mutate(max.uncertainty = model.taxon$max.uncertainty) %>%
        dplyr::select(Species, max.uncertainty, auc = `Test AUC`, auc.sd = `AUC Standard Deviation`) %>%
        mutate(Species = model.taxon$taxon.abbrv)
    df.res <- bind_rows(df.res, df)
}
f.out <- paste0(manuscript.dir, "tables/table_ST2.csv")
df.res %>% group_by(Species, max.uncertainty) %>%
    summarise(mean.auc = round(mean(auc), 3), mean.sd = round(mean(auc.sd), 3)) %>%
    write_csv(., f.out)
cat("  ==>", f.out, "\n")

#### FIGURE S1 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting Figure S1 ...\n")

perc.interval <- 20
sf::sf_use_s2(FALSE)
df.countries <- world %>% dplyr::select(iso_a2, name) %>%
    inner_join(results$country.numbers$country.non.georeferenced, by = c("iso_a2"="country.code")) %>%
    rename(ng=n) %>%
    inner_join(results$country.numbers$country.georeferenced.all %>%
                   dplyr::select(country.code, g=n), by = c("iso_a2"="country.code")) %>%
    inner_join(results$country.numbers$country.georeferenced.unc %>%
                   dplyr::select(country.code, gu=n), by = c("iso_a2"="country.code")) %>%
    bind_cols(area=st_area(.)) %>%
    mutate(nps = ng + g) %>%
    mutate(nps_alpha = (nps - min(nps)) / (max(nps) - min(nps))) %>%
    mutate(pg = round(g / nps *100, 1)) %>%
    mutate(pg.c = as.character(ifelse(pg == 100, 100/perc.interval - 1, floor(pg/perc.interval)))) %>%
    mutate(pgu = round(gu / g, 3) * 100) %>%
    mutate(pgu.c = as.character(ifelse(pgu == 100, 100/perc.interval - 1, floor(pgu/perc.interval)))) %>%
    mutate(area = round(as.double(area/10^6)),1) %>%
    rename(country.code = iso_a2)

# Fig S1a
breaks <- c(1 %o% 10^(1:floor(log10(max(df.countries$nps)))))
p.country.nps <- ggplot(df.countries) +
    geom_sf(aes(fill=nps), colour='grey30') +
    scale_fill_viridis(discrete = FALSE, trans = "log", breaks = breaks, labels = formatC(breaks, digits=0, format="E")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_colourbar(title = "Number of specimens",
                                  title.position = "top",
                                  label.theme = element_text(hjust=0.5, size=6),
                                  title.theme = element_text(face = 'bold', size=8, hjust = 0.5),
                                  barwidth = unit(7, "cm"),
                                  barheight = unit(0.25, "cm")))
p.country.nps
f.out <- paste0(manuscript.dir, "figures/figure_S1a.png")
ggsave(plot=p.country.nps, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig S1b
breaks <- c(0, 20, 40, 60, 80, 100)
p.country.perc.g <- ggplot(df.countries) +
    geom_sf(aes(fill=pg.c), colour = 'grey30') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(11,"Spectral")[c(6,8,4,2,11)], labels = c("[0-20)", "[20-40)", "[40-60)", "[60-80)", "[80-100]")) +
    scale_alpha_continuous(guide = 'none') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_legend(order = 1, title = "% records with coordinates"))
p.country.perc.g
f.out <- paste0(manuscript.dir, "figures/figure_S1b.png")
ggsave(plot=p.country.perc.g, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig S1c
p.country.perc.gu <- ggplot(df.countries) +
    geom_sf(aes(fill=pgu.c), colour = 'grey30') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(11,"Spectral")[c(6,8,4,2,11)], labels = c("[0-20)", "[20-40)", "[40-60)", "[60-80)", "[80-100]")) +
    scale_alpha_continuous(guide = 'none') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_legend(order = 1, title = "% records with uncertainty"))
p.country.perc.gu
f.out <- paste0(manuscript.dir, "figures/figure_S1c.png")
ggsave(plot=p.country.perc.gu, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

#### FIGURE S2 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S2 ...\n")

df.inst.geo <- world %>% dplyr::select(iso_a2, name) %>%
    left_join(
        pivot_wider(
            bind_rows(
                results$institution.numbers$inst.non.geo %>%
                    left_join(colonies, by = c("pub.country"="colony.code")) %>%
                    mutate(pub.country = ifelse(!is.na(country.code.y), country.code.y, pub.country)) %>%
                    group_by(pub.country) %>%
                    summarise(n = sum(n), .groups = 'drop') %>%
                    mutate(subset = "non.geo"),

                results$institution.numbers$inst.country.taxonomy.geo %>%
                    left_join(colonies, by = c("pub.country"="colony.code")) %>%
                    mutate(pub.country = ifelse(!is.na(country.code.y), country.code.y, pub.country)) %>%
                    group_by(pub.country) %>%
                    summarise(n = sum(n), .groups = "drop") %>%
                    mutate(subset = "geo")), names_from = subset, values_from = n) %>%
            # Need to do this, if not, next sum is NA if non.geo or geo are NA
            mutate(non.geo = ifelse(is.na(non.geo), 0, non.geo), geo = ifelse(is.na(geo), 0, geo)) %>%
            mutate(total = non.geo + geo) %>%
            mutate(pg = round(geo/total*100,1)) %>%
            mutate(pg.c = as.character(ifelse(pg == 100, 100/perc.interval - 1, floor(pg/perc.interval)))),
        by = c("iso_a2"="pub.country")) %>%
    # Need to set NAs to zero since there are countries coming from 'world' which do not have any preserved specimens
    mutate(non.geo = ifelse(is.na(non.geo), 0, non.geo), geo = ifelse(is.na(geo), 0, geo),
           total = ifelse(is.na(total), 0, total), pg = ifelse(is.na(pg), 0, pg),
           pg.c = ifelse(is.na(pg.c), 0, pg.c))

# Fig S2a
breaks <- c(1 %o% 10^(1:ceiling(log10(max(na.omit(df.inst.geo$total))))))
p.country.published <-
    ggplot() +
    geom_sf(data = df.inst.geo %>% filter(total != 0), aes(fill=total), colour='grey30') +
    geom_sf(data = df.inst.geo %>% filter(total == 0), fill='white', colour = 'grey30') +
    scale_fill_viridis(discrete = FALSE, trans = "log",
                       breaks = breaks, labels = formatC(breaks, digits=0, format="E"), na.value = 'grey92') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_colourbar(title = "Number of specimens",
                                  title.position = "top",
                                  label.theme = element_text(hjust=0.5, size=6),
                                  title.theme = element_text(face = 'bold', size=8, hjust = 0.5),
                                  barwidth = unit(7, "cm"),
                                  barheight = unit(0.25, "cm")))
p.country.published
f.out <- paste0(manuscript.dir, "figures/figure_S2a.png")
ggsave(plot=p.country.published, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig S2b
p.country.published.perc.g <-
    ggplot() +
    geom_sf(data = df.inst.geo %>% filter(geo != 0), aes(fill=pg.c), colour='grey30') +
    geom_sf(data = df.inst.geo %>% filter(geo == 0), fill='white', colour = 'grey30') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(11,"Spectral")[c(6,8,4,2,11)], labels = c("[0-20)", "[20-40)", "[40-60)", "[60-80)", "[80-100]")) +
    scale_alpha_continuous(guide = 'none') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_legend(order = 1, nrow = 1, title = "% records with coordinates"))
p.country.published.perc.g
f.out <- paste0(manuscript.dir, "figures/figure_S2b.png")
ggsave(plot=p.country.published.perc.g, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Fig S2c
df.inst.geo.unc <- world %>% dplyr::select(iso_a2, name) %>%
    left_join(
        results$institution.numbers$inst.country.taxonomy.geo %>%
            left_join(colonies, by = c("pub.country"="colony.code")) %>%
            mutate(pub.country = ifelse(!is.na(country.code.y), country.code.y, pub.country)) %>%
            group_by(pub.country, has.unc) %>%
            summarise(n = sum(n), .groups = "drop") %>%
            pivot_wider(names_from = has.unc, values_from = n) %>%
            rename("with.unc"="Y", "without.unc"="N") %>%
            mutate(with.unc = ifelse(is.na(with.unc), 0, with.unc)) %>%
            mutate(without.unc = ifelse(is.na(without.unc), 0, without.unc)) %>%
            mutate(total.geo = with.unc + without.unc) %>%
            mutate(perc.unc = round(with.unc / total.geo * 100, 1)) %>%
            mutate(perc.unc.c = as.character(ifelse(perc.unc == 100, 100/perc.interval - 1, floor(perc.unc/perc.interval)))),
        by = c("iso_a2"="pub.country")) %>%
    # Need to set NAs to zero since there are countries coming from 'world' which do not have any preserved specimens
    mutate(without.unc = ifelse(is.na(without.unc), 0, without.unc), with.unc = ifelse(is.na(with.unc), 0, with.unc),
           total.geo = ifelse(is.na(total.geo), 0, total.geo), perc.unc = ifelse(is.na(perc.unc), 0, perc.unc),
           perc.unc.c = ifelse(is.na(perc.unc.c), 0, perc.unc.c))

breaks <- c(0, 20, 40, 60, 80, 100)
p.country.published.perc.unc.g <-
    ggplot(df.inst.geo.unc) +
    geom_sf(data = df.inst.geo.unc %>% filter(with.unc != 0), aes(fill=perc.unc.c), colour='grey30') +
    geom_sf(data = df.inst.geo.unc %>% filter(with.unc == 0), fill='white', colour = 'grey30') +
    scale_fill_manual(values = RColorBrewer::brewer.pal(11,"Spectral")[c(6,8,4,2,11)], labels = c("[0-20)", "[20-40)", "[40-60)", "[60-80)", "[80-100]")) +
    scale_alpha_continuous(guide = 'none') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", panel.background = element_blank()) +
    guides(fill = guide_legend(order = 1, title = "% records with uncertainty"))
p.country.published.perc.unc.g
f.out <- paste0(manuscript.dir, "figures/figure_S2c.png")
ggsave(plot=p.country.published.perc.unc.g, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

#### FIGURE S3 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S3 ...\n")

f <- paste0(data.processed.dir, "gbif/", date.string, "_uncertainty_by_decimal_degree.csv")
df.unc.dd <- read_csv(f, col_types = cols())
breaks <- c(5, 50, 500, 5000, 50000, 500000, 5000000)
peaks.fig <- data.frame(cbind(peak=peaks,
                              subfig=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o",
                                       "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "aa", "ab", "ac", "ad",
                                       "ae", "af", "ag", "ah", "ai"))) %>% mutate(peak = as.numeric(peak))
for(i in 1:nrow(peaks.fig)){
    peak <- peaks.fig[i, "peak"]
    subfig <- peaks.fig[i, "subfig"]
    df.peak <- df.unc.dd %>% filter(peak == !!peak)
    br <- breaks[1:which(breaks > peaks[which(peaks > max(df.peak$n))[1]])[1]]


    title <- paste0('Uncertainty = ', formatC(peak, digits=0, format='f', big.mark = ' '), ' m. ',
                    '(Nr of specimens = ', formatC(sum(df.peak$n), digits=0, format='f', big.mark = ' '), ', ',
                    formatC(sum(df.peak$n)/sum(df.unc.dd$n) * 100, format = 'f', digits = 1), '%)')
    p.peak <- ggplot(df.peak) +
        geom_tile(aes(x = lon, y = lat, fill=n)) +
        geom_sf(data = world_coastline, fill = NA, colour = 'grey60', size = 0.2) +
        scale_x_continuous(expand = c(0.02, 0.02)) +
        scale_fill_viridis(trans = 'log', option = 'viridis', breaks = br) +
        annotate('text', x = -180, y = 88,
                 label = paste0('bold("', title, '")'),
                 parse = TRUE,
                 vjust = 0, hjust = 0, size = 2.5) +
        theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              legend.position = c(.7, .05), legend.background = element_rect(fill = NA),
              panel.grid.major = element_line(size = 0.2)) +
        guides(fill = guide_colourbar(title = "N",
                                      title.position = "left",
                                      direction = 'horizontal',
                                      title.theme = element_text(face = 'bold', size=7, hjust = 0, vjust = 0),
                                      label.theme = element_text(face = 'bold', size=6),
                                      label.position = 'top',
                                      barwidth = unit(5, "cm"),
                                      barheight = unit(0.15, "cm")))
    f.out <- paste0(manuscript.dir, "figures/figure_S3", subfig, ".png")
    ggsave(plot=p.peak, device="png", filename = f.out, width = 16, height = 16 * 0.6, units = "cm", dpi = 600)
    cat("  ==>", f.out, "\n")
}

#### FIGURE S4 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S4 ...\n")

breaks.unc <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
breaks.unc.labels <- formatC(breaks.unc, format = "e", digits = 0)

df.unc.country <- df.unc %>% filter(pub.country != 'UNKNOWN')
countries <- df.countries %>% st_drop_geometry() %>% dplyr::select(country.code, name) %>%
    mutate(name = ifelse(country.code == 'SO', 'Somalia', name)) %>%
    distinct()
df.unc.country <- df.unc.country %>% inner_join(countries, by = c("pub.country"="country.code"))

p.unc.country <- ggplot(df.unc.country,
                        aes(x = unc, y = name, fill = name)) +
    scale_x_continuous('Uncertainty (m)', trans="log10", breaks = breaks.unc, labels = breaks.unc.labels, expand = c(0,0)) +
    scale_y_discrete('Publishing country', expand = c(0,0), limits = rev) +
    geom_density_ridges(size = 0.2) +
    theme_ridges() +
    theme(legend.position = "none", axis.text = element_text(size = 7),
          axis.title = element_text(size = 8))

f.out <- paste0(manuscript.dir, "figures/figure_S4.png")
ggsave(plot=p.unc.country, device="png", filename = f.out, width = 16, height = 24, units = "cm", dpi = 300)
cat("  ==>", f.out, "\n")

#### FIGURE S5 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S5 ...\n")

kgd <- c('Animalia', 'Plantae', 'Fungi', 'Protozoa', 'Chromista', 'Bacteria', 'Archaea', 'Viruses', 'Unknown')
df.unc.kingdom <- df.unc %>%
    mutate(kingdom = ifelse(kingdom == "incertae sedis", "Unknown", kingdom)) %>%
    mutate(kingdom = factor(kingdom, levels = kgd))
p.unc.kingdom <- ggplot(df.unc.kingdom, aes(x = unc, y = kingdom, fill = kingdom)) +
    scale_x_continuous('Uncertainty (m)', trans="log10",
                       breaks = breaks.unc, labels = breaks.unc.labels, expand = c(0,0)) +
    scale_y_discrete('Kingdom', expand = c(0,0), limits = rev) +
    geom_density_ridges(size = 0.2) +
    theme_ridges() +
    theme(legend.position = "none", axis.text = element_text(size = 7),
          axis.title = element_text(size = 8))

f.out <- paste0(manuscript.dir, "figures/figure_S5.png")
ggsave(plot=p.unc.kingdom, device="png", filename = f.out, width = 16, height = 16, units = "cm", dpi = 300)
cat("  ==>", f.out, "\n")

#### FIGURE S6 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S6 ...\n")

# NUmber of records per kingdom and resolution
df1.res <- tibble()
for(kgd in c("Animalia", "Plantae", "Fungi")){
    df1 <- df.unc %>%
        filter(kingdom == kgd) %>%
        group_by(kingdom, unc.cat) %>%
        summarise(n=n(), .groups = "drop") %>%
        mutate(unc.cat = factor(unc.cat,
                                levels = c("<= 1m", "10m", "100m", "250m", "1km", "5km", "10km", "50km", "100km", "> 100km"))) %>%
        arrange(kingdom, unc.cat) %>%
        mutate(cs = cumsum(n), cum.perc=round(cs/sum(n)*100,1), .groups = 'drop')
    df1.res <- rbind(df1.res, df1)
}

# Number of unique taxon-location per kingdom and resolution
df2 <- df.unc %>%
    filter(kingdom %in% c("Animalia", "Plantae", "Fungi")) %>%
    dplyr::select(kingdom, taxon, lon, lat, unc, unc.cat) %>%
    group_by(kingdom, taxon, lon, lat, unc, unc.cat) %>%
    summarise(n=n(), .groups = "drop")

df2.res <- tibble()
for(kgd in c("Animalia", "Plantae", "Fungi")){
    df21 <- df2 %>% filter(kingdom == kgd) %>%
        group_by(kingdom, unc.cat) %>%
        summarise(n=n(), .groups = "drop") %>%
        mutate(unc.cat = factor(unc.cat,
                                levels = c("<= 1m", "10m", "100m", "250m", "1km", "5km", "10km", "50km", "100km", "> 100km"))) %>%
        arrange(kingdom, unc.cat) %>%
        mutate(cs = cumsum(n), cum.perc=round(cs/sum(n)*100,1), .groups = 'drop')
    df2.res <- rbind(df2.res, df21)
}

df3 <- df1.res %>% mutate(subset = "All records") %>%
    bind_rows(df2.res %>% mutate(subset = "Unique taxon-location combinations")) %>%
    dplyr::select(-.groups) %>%
    mutate(kingdom = factor(kingdom, levels = c("Animalia", "Plantae", "Fungi")))

p.bar.resolution.research <- ggplot() +
    geom_col(data=df3, aes(x = unc.cat, y = cs, fill = subset), width= 0.75,
             position = position_dodge2(width=0.9, reverse = T)) +
    geom_text(data=df3, aes(x = unc.cat, y = cs, colour = subset,
                            label = paste0(format(round(cs/1000000, 1), digits = 1,
                                                  scientific=F, big.mark = " "), paste0(", ", cum.perc, " %")),
                            angle = 90),
              position = position_dodge2(width=0.7, reverse = T), hjust = -0.05, vjust = 0.5, size = 3) +
    scale_y_continuous("Nr of specimens (in millions of records)",
                       limits=c(0, 43000000),
                       breaks=seq(0, 43000000, 5000000),
                       labels = paste0(format(seq(0, 43000000, 5000000)/1000000,
                                              scientific = F, big.mark = " "), "M"),
                       expand = c(0,0)) +
    scale_x_discrete("Target research resolution") +
    scale_fill_manual(values = c("grey40", "firebrick")) +
    scale_colour_manual(values = c("grey40", "firebrick")) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 14),
          axis.title.x = element_text(margin = margin(t = -7, r = 0, b = 0, l = 0)),
          panel.grid.major = element_blank(),
          axis.line = element_line(size = 0.3, colour = "grey30"),
          panel.background = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank()) +
    facet_wrap(~kingdom, ncol=1)
p.bar.resolution.research
f.out <- paste0(manuscript.dir, "figures/figure_S6.png")
ggsave(plot=p.bar.resolution.research, device="png", filename = f.out, width = 16, height = 16, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

#### FIGURE S7 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S7 ...\n")

for(kgd in c("Animalia", "Plantae", "Fungi")){
    cat("  plotting", kgd, "...\n")
    for(i in 1:length(grid.sizes)){
        gs <- grid.sizes[i]
        cat("    grid size", gs, "...\n")
        dfg <- df.res.kgd %>% filter(kingdom == kgd & grid.size == gs)
        min.n <- formatC(min(dfg$n), format = "f", digits = 0, big.mark = " ")
        max.n <- formatC(max(dfg$n), format = "f", digits = 0, big.mark = " ")
        sum.n <- formatC(sum(dfg$n), format = "f", digits = 0, big.mark = " ")
        br <- br[br< max(dfg$n)]
        br <- c(br, br[length(br)] * 10)
        p <- ggplot(dfg) +
            geom_tile(aes(x=xc, y=yc, fill =n)) +
            geom_sf(data = wc, fill = NA, colour = 'white', size = 0.2) +
            scale_x_continuous(breaks = c(-180, 180)) +
            scale_y_continuous(breaks = c(-89, 89)) +
            scale_fill_gradient("N", trans = "log", low = "yellow", high = "red", breaks = br, limits = c(br[1], br[length(br)])) +
            theme(panel.background = element_rect(fill = "grey20"),
                  axis.title = element_blank(),
                  legend.key.width = unit(0.2, "cm"),
                  legend.key.height = unit(0.45, "cm"),
                  legend.title = element_text(colour = "ivory", size = 6),
                  legend.text = element_text(size = 4.5, colour = "ivory"),
                  legend.position = c(0.965, 0.18),
                  panel.grid = element_line(size = 0.1, colour = "ivory"),
                  legend.background = element_rect(fill="transparent")) +
            guides(fill=guide_colourbar(title.vjust = 2, title.hjust = 0.7)) +
            annotate("text", x = -15000000, y = 9000000,
                     label = paste0("Resolution: ",
                                    ifelse(gs < 10000, formatC(gs, digits=0, format="f", big.mark = ""),
                                           formatC(gs, digits=0, format="f", big.mark = " ")),
                                    "m"),
                     colour = "ivory", size = 3) +
            annotate("text", x = 16500000, y = 9000000,
                     label = paste0("N=", formatC(sum.n, digits=0, format="f", big.mark = " ")),
                     colour = "ivory", size = 3)

        f.out <- paste0(manuscript.dir, "figures/figS7_", kgd, "_", formatC(gs, digits=0, format="f"), ".png")
        ggsave(plot=p, device="png", filename = f.out, width = 18, height = 9.3, units = "cm", dpi = 600)
        cat("      ==>", f.out, "\n")
    }
}

#### FIGURE S8 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S8 ...\n")

plotRegion <- function(df.res.kgd, region.name, region.extent, region.coastline, region.plot.width,
                       region.plot.height){
    br <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
    dfg <- df.res.kgd %>% filter(xc >= region.extent["xmin"] & xc <= region.extent["xmax"] &
                                     yc >= region.extent["ymin"] & yc <= region.extent["ymax"])
    for(kgd in c("Animalia", "Plantae", "Fungi")){
        cat("  plotting", kgd, "...\n")
        for(i in 1:length(grid.sizes)){
            gs <- grid.sizes[i]
            cat("    grid size", gs, "...\n")
            dfg.kg <- dfg %>% filter(kingdom == kgd & grid.size == gs)
            min.n <- formatC(min(dfg.kg$n), format = "f", digits = 0, big.mark = " ")
            max.n <- formatC(max(dfg.kg$n), format = "f", digits = 0, big.mark = " ")
            sum.n <- formatC(sum(dfg.kg$n), format = "f", digits = 0, big.mark = " ")
            br <- br[br< max(dfg.kg$n)]
            br <- c(br, br[length(br)] * 10)
            p <- ggplot(dfg.kg) +
                geom_tile(aes(x=xc, y=yc, fill =n)) +
                geom_sf(data = region.coastline, fill = NA, colour = 'white', size = 0.2) +
                scale_fill_gradient("N", trans = "log", low = "yellow", high = "red", breaks = br, limits = c(br[1], br[length(br)])) +
                ggtitle(paste0("Resolution: ",
                               ifelse(gs < 10000, formatC(gs, digits=0, format="f", big.mark = ""),
                                      formatC(gs, digits=0, format="f", big.mark = " ")),
                               "m")) +
                theme(panel.background = element_rect(fill = "grey20"),
                      axis.title = element_blank(),
                      legend.key.width = unit(0.1, "cm"),
                      legend.key.height = unit(0.20, "cm"),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      legend.title = element_text(colour = "ivory", size = 6),
                      legend.text = element_text(size = 4.5, colour = "ivory"),
                      legend.position = c(0.89, 0.25),
                      panel.grid = element_blank(),
                      legend.background = element_rect(fill="grey20"),
                      plot.title = element_text(size = 6, hjust = 0, vjust = -2)) +
                guides(fill=guide_colourbar(title.vjust = 2, title.hjust = 0.7))
            p
            f.out <- paste0(manuscript.dir, "figures/figS8_", region.name, "_", tolower(kgd), "_", formatC(gs, digits=0, format="f"), ".png")
            ggsave(plot=p, device="png", filename = f.out, width = region.plot.width, height = region.plot.height, units = "cm", dpi = 600)
            cat("      ==>", f.out, "\n")
        }
    }
}

# Australia
region.name <- "australia"
region.extent <- c("xmin"=9500000, "xmax"=14500000, "ymin"=-5500000, "ymax"=-1000000)
region.coastline <- st_crop(st_transform(world_coastline, 8857), region.extent)
region.plot.width <- 6
region.plot.height <- 5.7
cat("Plotting Australia ...\n")
plotRegion(df.res.kgd, region.name, region.extent, region.coastline, region.plot.width, region.plot.height)
# Europe
region.name <- "europe"
region.extent <- c("xmin"=-1000000, "xmax"=2500000, "ymin"=4240000, "ymax"=7800000)
region.coastline <- st_crop(st_transform(world_coastline, 8857), region.extent)
region.plot.width <- 6
region.plot.height <- 5.7
cat("Plotting Europe ...\n")
plotRegion(df.res.kgd, region.name, region.extent, region.coastline, region.plot.width, region.plot.height)
# North America
region.name <- "north_america"
region.extent <- c("xmin"=-10700000, "xmax"=-5000000, "ymin"=1900000, "ymax"=6200000)
region.coastline <- st_crop(st_transform(world_coastline, 8857), region.extent)
region.plot.width <- 6
region.plot.height <- 4.6
cat("Plotting North America ...\n")
plotRegion(df.res.kgd, region.name, region.extent, region.coastline, region.plot.width, region.plot.height)

#### FIGURE S9 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S9 ...\n")
df.res.diff <- data.frame()
taxons <- list.filter(taxons.to.model, max.uncertainty == Inf)
for(i in 1:length(taxons)){
    model.taxon <- taxons[[i]]
    modelling.dir <- paste0(output.dir, "maxent/", model.taxon$taxon.abbrv, "-uncertainty")
    f.out <- paste0(modelling.dir, "/data/", model.taxon$taxon.abbrv, "_", model.taxon$region.name, "_presences_env_data.csv")
    env.df <- vroom::vroom(f.out, col_types = cols()) %>% filter(!is.na(unc))
    predictors <- env.df %>% dplyr::select(-id, -pa, -lon, -lat, -unc, -cell) %>% names()
    predictors <- paste0("bio", sort(na.omit(as.numeric(unlist(str_split(predictors, "bio"))))))
    if(model.taxon$taxon.abbrv == "guazuma_ulmifolia")
        predictors <- c(predictors, "treecover")
    env.df <- env.df %>%
        mutate(unc.cat = cut(unc, breaks=c(0, 1000, 5000, 10000, 50000, 100000, 1000000),
                             labels = c("<=1k", "(1k-5k]", "(5k-10k]", "(10k-50k]", "(50k-100k]", ">100k"))) %>%
        dplyr::select(id, unc.cat, all_of(predictors)) %>%
        pivot_longer(cols = all_of(predictors), names_to = "variable") %>%
        filter(!is.na(value))

    diff.df <- env.df %>%
        group_by(id, unc.cat, variable) %>%
        summarise(min = min(value), max = max(value), .groups = "drop") %>%
        mutate(diff=max - min, species=model.taxon$taxon.name)
    df.res.diff <- rbind(df.res.diff, diff.df)
}
df.res.diff <- df.res.diff %>% mutate(species = ifelse(str_starts(species, "Rho"), "Rhododendron groenlandicum", species),
                                      species = ifelse(str_starts(species, "Gua"), "Guazuma ulmifolia", species),
                                      species = ifelse(str_starts(species, "Euc"), "Eucalyptus gongylocarpa", species)) %>%
    mutate(species = factor(species,
                            levels = c("Rhododendron groenlandicum", "Guazuma ulmifolia", "Eucalyptus gongylocarpa")))

# Predictors common to all species
common.predictors <- c("bio2", "bio3", "bio8", "bio9", "bio18", "bio19")
for(i in 1:length(common.predictors)){
    predictor <- common.predictors[i]
    f.out <- paste0(manuscript.dir, "figures/figure_S9A", i, ".png")
    p <- ggplot(df.res.diff %>% filter(variable == predictor)) +
        geom_point(aes(x=unc.cat, y=diff, group = species, colour = species),
                   size = 1, alpha = 0.5, position = position_dodge(width = 0.5)) +
        scale_y_continuous("Variable range") +
        scale_x_discrete("Uncertainty interval") +
        scale_colour_manual(values = c("firebrick", "darkgreen", "blue")) +
        theme(panel.background = element_blank(),
              legend.position = "bottom",
              panel.grid.major = element_line(linetype = "dashed", colour = "grey70", size = 0.2),
              axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1.2),
              axis.text.y = element_text(size = 6),
              axis.title = element_text(size = 7),
              strip.text = element_text(size = 8),
              axis.ticks = element_blank(),
              plot.title = element_text(size = 8)) +
        guides(colour = guide_legend(title.position = "top", order = 1, title = "Species",
                                     override.aes = list(size = 3),
                                     label.theme = element_text(hjust=0.5, size=8),
                                     title.theme = element_text(size=9, hjust = 0.5))) +
        ggtitle(toupper(predictor))
    # print(p + theme(legend.position = "none"))
    ggsave(plot=p + theme(legend.position = "none"), device="png", filename = f.out, width = 6, height = 4, units = "cm", dpi = 300)
    cat("  ==>", f.out, "\n")
}

# Plot the legend separately
leg <- get_legend(p)
f.out <- paste0(manuscript.dir, "figures/figure_S9A_legend.png")
ggsave(plot=leg, device="png", filename = f.out, width = 12, height = 2, units = "cm", dpi = 600)
cat("  ==>", f.out, "\n")

# Plot the rest of predictors, exclusive to each species
l.exclusive.predictors = list(
    rg = c("bio4", "bio15"),
    gu = c("bio13", "bio14", "bio15", "treecover"),
    eg = c("bio14", "bio15")
)
sps <- c("rg", "gu", "eg")
for(i in 1:length(sps)){
    predictors <- l.exclusive.predictors[[sps[i]]]
    for(j in 1:length(predictors)){
        predictor <- predictors[j]
        f.out <- paste0(manuscript.dir, "figures/figure_S9", toupper(base::letters[i+1]), j, ".png")
        p <- ggplot(df.res.diff %>% filter(variable == predictor)) +
            geom_point(aes(x=unc.cat, y=diff, group = species),
                       size = 1, alpha = 0.3, position = position_dodge(width = 0.5)) +
            scale_y_continuous("Variable range") +
            scale_x_discrete("Uncertainty interval") +
            theme(panel.background = element_blank(),
                  legend.position = "none",
                  panel.grid.major = element_line(linetype = "dashed", colour = "grey70", size = 0.2),
                  axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1.2),
                  axis.text.y = element_text(size = 6),
                  axis.title = element_text(size = 7),
                  strip.text = element_text(size = 8),
                  axis.ticks = element_blank(),
                  plot.title = element_text(size = 8)) +
            ggtitle(toupper(predictor))
        # print(p)
        ggsave(plot=p, device="png", filename = f.out, width = 6, height = 4, units = "cm", dpi = 300)
        cat("  ==>", f.out, "\n")
    }
}

#### FIGURE S10 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S10 ...\n")

# Generate csv files for suitability sd for each taxon
taxons.s10 <- list.filter(taxons.to.model, max.uncertainty == Inf)
for(taxon in taxons.s10){
    cat("  processing", taxon$taxon.name, "suitability sd ...\n")
    f.out <- paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_Inf/", taxon$taxon.abbrv, "_prediction_sd")
    if(!file.exists(paste0(f.out, ".csv"))){
        f <- paste0(f.out, ".tif")
        r <- raster(f)
        if(str_starts(taxon$taxon.abbrv, "rhododendron")){ # No need to mask for Guazuma or Eucalyptus
            r <- r * raster(paste0(data.processed.dir, taxon$region.name, "_mask.tif"))
        }
        df <- r %>% as.data.frame(., xy=TRUE, cells=FALSE, na.rm=TRUE) %>% setNames(., c("x", "y", "sd"))
        write_csv(df, file = paste0(f.out, ".csv"))
        cat("    ==>", f.out, "\n")
    }
}

# Plot Rhododendron groenlandicum
taxon <- taxons.s10[[1]]
cat("  plotting", taxon$name, "standard deviation ...\n" )
df <- read_csv(paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_Inf/", taxon$taxon.abbrv, "_prediction_sd.csv"),
               col_types = cols())
p <- ggplot(df) +
    geom_tile(aes(x=x, y=y, fill=sd)) +
    scale_x_continuous(breaks = seq(-170, -10, 10)) +
    scale_y_continuous(breaks = seq(40, 90, 10)) +
    scale_fill_viridis(limits=c(0,0.3), breaks=seq(0, 0.3, 0.05)) +
    coord_equal() +
    theme(axis.title = element_blank(),
          panel.grid.major = element_line(colour = "grey70", linetype = 'dashed', size = 0.2),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Standard deviation",
                                  title.position = "top",
                                  label.theme = element_text(hjust=0.5, size=8),
                                  title.theme = element_text(size=9, hjust = 0.5),
                                  barwidth = unit(7, "cm"),
                                  barheight = unit(0.25, "cm"),
                                  direction = "horizontal"))
f.out <- paste0(manuscript.dir, "figures/figS10_", taxon$taxon.abbrv, "_suitability-sd.tif")
ggsave(plot=p, device="tiff", filename = f.out, width = 16, height = 10, units = "cm", dpi = 600)
cat("    ==>", f.out, "\n")

# Plot Guazuma ulmifolia
taxon <- taxons.s10[[2]]
cat("  plotting", taxon$name, "standard deviation ...\n" )
df <- read_csv(paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_Inf/", taxon$taxon.abbrv, "_prediction_sd.csv"),
               col_types = cols())
p <- ggplot(df) +
    geom_tile(aes(x=x, y=y, fill=sd)) +
    scale_x_continuous(breaks = seq(-120, -30, 10)) +
    scale_y_continuous(breaks = seq(-40, 40, 10)) +
    scale_fill_viridis(limits=c(0,0.3), breaks=seq(0, 0.3, 0.05)) +
    coord_equal() +
    theme(axis.title = element_blank(),
          panel.grid.major = element_line(colour = "grey70", linetype = 'dashed', size = 0.2),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Standard deviation",
                                  title.position = "top",
                                  label.theme = element_text(hjust=0.5, size=8),
                                  title.theme = element_text(size=9, hjust = 0.5),
                                  barwidth = unit(7, "cm"),
                                  barheight = unit(0.25, "cm"),
                                  direction = "horizontal"))
f.out <- paste0(manuscript.dir, "figures/figS10_", taxon$taxon.abbrv, "_suitability-sd.tif")
ggsave(plot=p, device="tiff", filename = f.out, width = 16, height = 13, units = "cm", dpi = 600)
cat("    ==>", f.out, "\n")

# Plot Eucalyptus gongylocarpa
taxon <- taxons.s10[[3]]
cat("  plotting", taxon$name, "standard deviation ...\n" )
df <- read_csv(paste0(output.dir, "maxent/", taxon$taxon.abbrv, "-uncertainty/predictions/maxunc_Inf/", taxon$taxon.abbrv, "_prediction_sd.csv"),
               col_types = cols())
p <- ggplot(df) +
    geom_tile(aes(x=x, y=y, fill=sd)) +
    scale_x_continuous(breaks = seq(110, 155, 5)) +
    scale_y_continuous(breaks = seq(-40, -10, 5)) +
    scale_fill_viridis(limits=c(0,0.3), breaks=seq(0, 0.3, 0.05)) +
    coord_equal() +
    theme(axis.title = element_blank(),
          panel.grid.major = element_line(colour = "grey70", linetype = 'dashed', size = 0.2),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title = "Standard deviation",
                                  title.position = "top",
                                  label.theme = element_text(hjust=0.5, size=8),
                                  title.theme = element_text(size=9, hjust = 0.5),
                                  barwidth = unit(7, "cm"),
                                  barheight = unit(0.25, "cm"),
                                  direction = "horizontal"))
f.out <- paste0(manuscript.dir, "figures/figS10_", taxon$taxon.abbrv, "_suitability-sd.tif")
ggsave(plot=p, device="tiff", filename = f.out, width = 16, height = 14, units = "cm", dpi = 600)
cat("    ==>", f.out, "\n")

#### FIGURE S11 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S12 ...\n")

# Just need to call functions done for FIGURE 7
unc.thresholds <- c(3536, 7071, 14142, 28284)
for(i in 1:length(unc.thresholds)){
    tic()
    generateTaxonFigures("eucalyptus_gongylocarpa", unc.thresholds[i], "S11")
    generateTaxonFigures("rhododendron_groenlandicum", unc.thresholds[i], "S11")
    generateTaxonFigures("guazuma_ulmifolia", unc.thresholds[i], "S11")
    toc()
}

#### FIGURE S12 -------------------------------------------------------------------------------- ####
cat("\n", rep("-", 80), "\n", sep = "")
cat("Plotting figure S11 ...\n")

unc.thresholds <- c("3536", "7071", "14142", "28284", "Inf")
taxons <- c("eucalyptus_gongylocarpa", "guazuma_ulmifolia", "rhododendron_groenlandicum")
taxon.names <- c("Eucalyptus gongylocarpa", "Guazuma ulmifolia", "Rhododendron groenlandicum")
names(taxons) <- taxon.names
res <- tibble()
for(i in 1:length(taxons)){
    taxon <- taxons[i]
    for(unc in unc.thresholds){
        f.range.results <- paste0("outputs/maxent/", taxon, "-uncertainty/", taxon, "_maxunc_", unc, "_suitability_range_results.csv")
        range.results.df <- vroom(f.range.results, col_types = cols())
        m.min <- getQuantileModel(range.results.df, 0.05)
        m.max <- getQuantileModel(range.results.df, 0.95)
        min.area <- range.results.df %>% filter(model.id == m.min$model.id) %>% pull(range.area)
        max.area <- range.results.df %>% filter(model.id == m.max$model.id) %>% pull(range.area)
        res <- bind_rows(res, bind_cols(taxon.name = names(taxon), taxon.abbrv = taxon,
                                        unc = unc, min.area = min.area, max.area = max.area)) %>%
            mutate(delta = abs(max.area - min.area),
                   perc.increment = round(delta / min.area * 100, 1),
                   unc = factor(unc, levels = unc.thresholds),
                   taxon.abbrv = factor(taxon.abbrv, levels = taxons),
                   taxon.name = factor(taxon.name, levels = taxon.names))

    }
}

fig.theme <- theme(axis.text.x = element_text(size = 6, hjust = 1, vjust = 1.2),
                   axis.text.y = element_text(size = 6),
                   axis.title = element_text(size = 7),
                   strip.text = element_text(size = 8),
                   axis.ticks = element_blank(),
                   legend.text = element_text(size = 6),
                   legend.title = element_text(size = 6),
                   plot.title = element_text(size = 8))
p <- ggplot(res) +
    geom_point(aes(x = unc, y = delta, colour = taxon.name, shape = taxon.name), size = 3) +
    scale_shape_manual("Species", values = c(0, 5, 6)) +
    scale_colour_manual("Species", values = c("firebrick", "darkblue", "darkgreen")) +
    scale_x_discrete("Maximum uncertainty (km)") +
    scale_y_continuous("Number of 30 arc-second grid cells (~ 1km)") +
    fig.theme
p
f.out <- paste0(manuscript.dir, "figures/figS12_range_differences.tif")
ggsave(plot=p, device="tiff", filename = f.out, width = 12, height = 6, units = "cm", dpi = 300)
cat("    ==>", f.out, "\n")

p <- ggplot(res) +
    geom_point(aes(x = unc, y = perc.increment, colour = taxon.name, shape = taxon.name), size = 3) +
    scale_shape_manual("Species", values = c(0, 5, 6)) +
    scale_colour_manual("Species", values = c("firebrick", "darkblue", "darkgreen")) +
    scale_x_discrete("Maximum uncertainty (km)") +
    scale_y_continuous("Percentage") +
    fig.theme
p
f.out <- paste0(manuscript.dir, "figures/figS12_percentages_increment.tif")
ggsave(plot=p, device="tiff", filename = f.out, width = 12, height = 6, units = "cm", dpi = 300)
cat("    ==>", f.out, "\n")






