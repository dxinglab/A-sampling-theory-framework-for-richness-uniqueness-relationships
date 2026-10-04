# Plot complete BCR-specific expectations for Figures S2 and S3.
# Run after S21-prepare-bbs-full-expectations.R.

# 1. Paths and packages --------------------------------------------------

input_file <- "process data/FiguresS2_S3_BBS_full_expectations.rds"

library(ggplot2)

# 2. Read the plotting data ---------------------------------------------

plot_data <- readRDS(input_file)

site_source <- plot_data$site_data
site <- rbind(
    data.frame(
        BCR = site_source$BCR, richness = site_source$richness,
        lcbd = site_source$coverage_lcbd, framework = "Coverage"
    ),
    data.frame(
        BCR = site_source$BCR, richness = site_source$richness,
        lcbd = site_source$sorensen_lcbd, framework = "Sorensen"
    )
)
curve <- plot_data$curves
curve$model[curve$model == "Annealed expectation"] <-
    "Conditional expectation"
curve$model <- factor(
    curve$model,
    levels = c("Parameterized expectation", "Conditional expectation")
)

# 3. Prepare independent log-richness axes -------------------------------

site$x_plot <- log10(site$richness + 1)
curve$x_plot <- log10(curve$richness + 1)
bcr_levels <- paste("BCR", sort(unique(site$BCR)))
site$BCR_label <- factor(paste("BCR", site$BCR), levels = bcr_levels)
curve$BCR_label <- factor(paste("BCR", curve$BCR), levels = bcr_levels)

# 4. One 28-BCR block per framework -------------------------------------

# Draw all BCR-specific annealed fits for one LCBD framework.
make_framework_plot <- function(framework) {
    site_now <- site[site$framework == framework, ]
    curve_now <- curve[curve$framework == framework, ]

    ggplot(site_now, aes(x_plot, lcbd)) +
        geom_hex(
            aes(fill = after_stat((log1p(count) / log1p(max(count)))^2.15)),
            bins = 18, colour = NA
        ) +
        scale_fill_gradientn(
            colours = c("#EEEAF2", "#C8BED5", "#9A89B0", "#6D5A8E"),
            limits = c(0, 1), guide = "none"
        ) +
        geom_line(
            data = curve_now,
            aes(
                x_plot, expected, group = model,
                linetype = model, linewidth = model
            ),
            colour = "#6D5A8E", lineend = "round",
            inherit.aes = FALSE
        ) +
        facet_wrap(
            vars(BCR_label), ncol = 4, scales = "free_x",
            axes = "all_x", axis.labels = "all_x"
        ) +
        scale_x_continuous(
            breaks = function(limits) {
                seq(limits[1L], limits[2L], length.out = 5L)[2:4]
            },
            labels = function(x) format(round(10^x - 1), trim = TRUE),
            expand = expansion(mult = c(0, 0.015))
        ) +
        scale_y_continuous(
            breaks = seq(0, 1, 0.2),
            expand = expansion(mult = c(0, 0))
        ) +
        coord_cartesian(ylim = c(0, 1.05)) +
        scale_linetype_manual(
            values = c(
                "Parameterized expectation" = "11",
                "Conditional expectation" = "solid"
            ),
            breaks = c(
                "Parameterized expectation", "Conditional expectation"
            ),
            name = NULL
        ) +
        scale_linewidth_manual(
            values = c(
                "Parameterized expectation" = 0.42,
                "Conditional expectation" = 0.78
            ),
            guide = "none"
        ) +
        labs(
            x = "Local richness", y = expression(LCBD~"\u00d7"~M)
        ) +
        theme_classic(base_size = 15, base_family = "sans") +
        theme(
            panel.background = element_rect(fill = "white", colour = NA),
            panel.border = element_rect(
                fill = NA, colour = "black", linewidth = 0.98
            ),
            axis.line = element_blank(),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 11, colour = "grey25"),
            strip.background = element_blank(),
            strip.text = element_text(size = 13, face = "bold"),
            legend.position = "bottom",
            legend.key.width = grid::unit(1.7, "cm"),
            legend.text = element_text(size = 13),
            panel.spacing.x = grid::unit(0.55, "lines"),
            panel.spacing.y = grid::unit(0.32, "lines"),
            plot.margin = margin(6, 8, 4, 7)
        )
}

# 5. Export publication and editable formats ----------------------------

for (framework in c("Coverage", "Sorensen")) {
    plot <- make_framework_plot(framework)
    figure_number <- if (framework == "Coverage") "S2" else "S3"
    framework_stem <- paste0(
        "figure/Figure", figure_number,
        "_BBS_annealed_all_BCR_", framework
    )
    ggsave(
        paste0(framework_stem, ".png"), plot,
        width = 10, height = 15.5, dpi = 400, bg = "white"
    )
    ragg::agg_tiff(
        paste0(framework_stem, ".tiff"), width = 10, height = 15.5,
        units = "in", res = 600, background = "white", compression = "lzw"
    )
    print(plot)
    grDevices::dev.off()
    grDevices::pdf(
        paste0(framework_stem, ".pdf"),
        width = 10, height = 15.5,
        family = "Helvetica", useDingbats = FALSE, bg = "white"
    )
    print(plot)
    grDevices::dev.off()
}
