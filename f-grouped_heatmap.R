# Create a grouped heatmap

library(KRSA)
library(tidyverse)

pep_data <- list.files("datastore", "data", full.names = TRUE) |>
  set_names(~ str_extract(.x, "run\\d+")) |>
  map(readRDS) |>
  map(~ pluck(.x, "grouped")) |>
  map(~ filter(.x, r.seq >= 0.8, str_detect(Peptide, "REF", negate = TRUE), str_detect(Peptide, "^p", negate = TRUE))) |>
  list_rbind(names_to = "run") |>
  select(-r.seq) |>
  summarise(slope = mean(slope), .by = c(Peptide, Group))

peps <- pep_data |>
  pull(Peptide) |>
  unique()

krsa_heatmap_grouped(pep_data, peps, scale = "row")

svg(filename = "figures/combined_heatmap.svg")
krsa_heatmap_grouped(pep_data, peps, scale = "row")
dev.off()
