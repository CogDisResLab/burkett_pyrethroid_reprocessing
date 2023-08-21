# Create barplot for Kinases of Interest

library(tidyverse)
library(ggprism)

pep_data <- list.files("datastore", "data", full.names = TRUE) |>
  set_names(~ str_extract(.x, "run\\d+")) |>
  map(readRDS) |>
  map(~ pluck(.x, "normalized")) |>
  map(~ filter(.x, r.seq >= 0.8, str_detect(Peptide, "REF", negate = TRUE))) |>
  list_rbind(names_to = "run") |>
  mutate(Group = if_else(Group == "A", "Exposed", "Control"))

mapping <- KRSA::KRSA_Mapping_STK_PamChip_87102_v1

kinases_of_interest <- read_csv("results/combined_quartile_ranked_KRSA.csv") |>
  select(run, Kinase, quartile) |>
  pivot_wider(names_from = run, values_from = quartile) |>
  mutate(mean_rank = (run01 + run02 + run03) / 3) |>
  slice_min(mean_rank, n = 10) |>
  pull(Kinase)

combined_data <- mapping |>
  separate_longer_delim(Kinases, delim = " ") |>
  filter(Kinases %in% kinases_of_interest) |>
  rename(Peptide = Substrates, Kinase = Kinases) |>
  inner_join(pep_data) |>
  select(SampleName, Group, Kinase, Peptide, Run = run, Slope = slope)

tested <- combined_data |>
  nest(.by = c(Run, Kinase)) |>
  mutate(ttest = map(data, ~ t.test(formula = Slope ~ Group, data = .x)),
         glanced = map(ttest, broom::glance),
         mean_val = map(data, ~ summarise(.x, mean_value = mean(Slope), .by = Group)),
         mean_val = map(mean_val, ~ pivot_wider(.x, names_from = Group, values_from = mean_value, names_prefix = "MeanVal"))) |>
  select(-ttest, -data) |>
  unnest_wider(glanced) |>
  unnest_wider(mean_val) |>
  select(-method, -alternative, -parameter, -estimate1, -estimate2) |>
  mutate(Significant = p.value < 0.05) |>
  pivot_longer(cols = starts_with("MeanVal"), names_to = "Group", values_to = "Mean") |>
  mutate(Group = str_remove(Group, "MeanVal")) |>
  write_csv("results/t-test_all_kinases.csv")

tested_filtered <- tested |>
  filter(Kinase %in% c("PKCH", "AKT", "PAKB", "PIM"), Run == "run01") |>
  write_csv("results/t-test_selected_kinases.csv")

g <- ggplot(tested_filtered, aes(x = Kinase, y = Mean, fill = Group, ymin = Mean - conf.low, ymax = Mean + conf.high))

p <- g +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(position = position_dodge(0.8), width = 0.7) +
  scale_fill_prism("shades_of_gray") +
  theme_prism()

ggsave("figures/barplot-comparison.svg", plot = p)
