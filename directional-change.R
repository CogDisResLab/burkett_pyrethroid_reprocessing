# Get Direction of Top Kinases

library(tidyverse)
library(KRSA)

sanitize_names <- function(names) {
  names |>
    str_extract("(run.*)-dpp_.*-STK\\.csv", 1) |>
    str_replace_all("_", "-")
}

threshold <- 0.2

chip_coverage <- KRSA::KRSA_coverage_STK_PamChip_87102_v2 |>
  rename(Peptide = Substrates, Kinase = Kin)

dpp_files <- list.files("results", "dpp") |>
  set_names(~ sanitize_names(.x)) |>
  map(~ read_csv(file.path("results", .x), show_col_types = FALSE)) |>
  map(~ inner_join(.x, chip_coverage, by = "Peptide")) |>
  imap(~ mutate(.x, Comparison = .y)) |>
  map(~ select(.x, Barcode, Peptide, LFC, Kinase, Comparison)) |>
  map(~ mutate(.x, Score = case_when(
    LFC < -threshold ~ -1,
    LFC > threshold ~ 1,
    .default = 0
  ))) |>
  bind_rows() |>
  nest(.by = c(Comparison, Kinase)) |>
  filter(Kinase %in% c("BRSK", "ERK", "HAL", "MLCK", "PAKB", "PIM", "PKCH", "PKN", "RAD53", "WNK"), Comparison == "run-all") |>
  mutate(Direction = map_dbl(data, ~ mean(.x$Score))) |>
  select(-data) |>
  write_csv("results/selected_kinase_activity_directions.csv")

