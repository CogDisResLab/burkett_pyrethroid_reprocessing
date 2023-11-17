# Make Quartile Plots

library(tidyverse)

quartile_figure <- function(df, grouping = "KinaseFamily") {
    df %>%
        dplyr::select(hgnc_symbol, one_of(grouping), Qrt, Method) %>%
        tidyr::pivot_wider(
            names_from = Method,
            values_from = Qrt,
            values_fn = unique
        ) %>%
        tidyr::pivot_longer(3:ncol(.), names_to = "Method", values_to = "Qrt") %>%
        dplyr::mutate(
            present = ifelse(is.na(Qrt), "No", "Yes"),
            Qrt = ifelse(present == "No", 2, Qrt),
            present = factor(present, levels = c("Yes", "No")),
            Qrt = factor(Qrt, levels = c(1, 2, 3, 4)),
            Method = factor(Method, levels = c("UKA", "PTM-SEA", "KEA3", "KRSA"))
        ) %>%
        ggplot2::ggplot(ggplot2::aes(hgnc_symbol, Method)) +
        ggplot2::geom_point(ggplot2::aes(size = Qrt, shape = present)) +
        ggplot2::scale_size_manual(values = c(
            `4` = 4,
            `3` = 3,
            `2` = 2,
            `1` = 1
        )) +
        ggplot2::theme_bw() +
        {
            if (grouping == "subfamily") {
                ggplot2::facet_grid(. ~ subfamily, scales = "free", space = "free")
            } else if (grouping == "group") {
                ggplot2::facet_grid(. ~ group, scales = "free", space = "free")
            } else {
                ggplot2::facet_grid(. ~ KinaseFamily,
                    scales = "free",
                    space = "free"
                )
            }
        } +
        ggplot2::scale_shape_manual(values = c(Yes = 19, No = 1)) +
        ggplot2::theme(axis.text.x = element_text(
            angle = 30,
            size = 7.5,
            vjust = 0.7
        )) +
        ggplot2::labs(x = "", y = "")
}

generate_quartile_plot <- function(datafile) {
    creeden_data <-
        read_csv(file.path("results", datafile), show_col_types = F)

    sig_kinases <- creeden_data |>
        filter(Method == "KRSA", Qrt == 4) |>
        pull(hgnc_symbol) |>
        unique()

    creeden_data |>
        filter(hgnc_symbol %in% sig_kinases) |>
        quartile_figure()
}

creedenzymatic_files <- list.files("results", "creedenzymatic") |>
    set_names(~ str_remove(.x, "_.*")) |>
    map(generate_quartile_plot) |>
    imap(~ ggsave(
        str_glue("{.y}-creedenzymatic.png"),
        path = "figures",
        plot = .x,
        width = 20,
        height = 5
    ))
