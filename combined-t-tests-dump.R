# Dump all data for one-sample t-tests

library(tidyverse)
library(KRSA)

pep_data <- readRDS("datastore/run-all-data_modeled-STK.RDS") |>
  pluck("grouped") |>
  filter(r.seq >= 0.8, str_detect(Peptide, "REF", negate = TRUE), str_detect(Peptide, "^p", negate = TRUE)) |>
  mutate(Group = if_else(Group == "A", "Exposed", "Control")) |>
  select(-r.seq)


mapping <- KRSA_Mapping_STK_PamChip_87102_v1 |>
  rename(Peptide = Substrates, Kinase = Kinases) |>
  separate_longer_delim(Kinase, " ") |>
  unique()

all_pep <- pep_data |>
  left_join(mapping, by = "Peptide") |>
  write_csv("results/all_kinase_peptide_Intensity.csv")


all_dpp <- read_csv("results/run-all-dpp_A_B-STK.csv") |>
  select(-totalMeanLFC, -LFC_SD) |>
  left_join(mapping, by = "Peptide") |>
  write_csv("results/all_kinase_peptide_LFC.csv")
