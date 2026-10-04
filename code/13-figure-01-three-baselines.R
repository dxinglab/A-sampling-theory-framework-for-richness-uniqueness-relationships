# Compare the three LCBD baselines at a specified N parameter.

# 1. Settings --------------------------------------------------------------

N <- 211250
M <- 1250L
x_values <- c(0.999, 0.9999, 0.99999)
conditional_x <- 0.9999
xprime_values <- c(0.5, 0.99)
coverage_result_file <- "process data/coverage_fixed_N_parameter_space.rds"
sorensen_result_file <- "process data/sorensen_fixed_N_parameter_space.rds"
sad_file <- "process data/fixed_N_fisher_sads.rds"
output_base <- "figure/Figure1_three_baselines"

figure_size <- list(width = 12, height = 16.7, dpi = 400L)
figure_style <- list(
    text = 18,
    axis_text = 16,
    legend = 14,
    tag = 16,
    border = 0.55,
    line = 0.70,
    point = 2.65,
    diamond = 4.50
)

draw_key_xprime_point <- function(data, params, size) {
    data$fill <- data$colour
    data$colour <- "white"
    data$shape <- 23
    data$size <- figure_style$diamond
    data$stroke <- 0.75
    ggplot2::draw_key_point(data, params, size)
}

font_family <- "Helvetica"
xprime_colours <- c("0.5" = "#65549A", "0.99" = "#267AA2")
xprime_marker_colours <- c("0.5" = "#8E80B5", "0.99" = "#5D9CB8")
heatmap_colours <- c("#5E3C99", "#F7F7F7", "#2B8CBE")
x_ticks <- c(1, 3, 10, 30, 100, 300, 1000)

library(ggplot2)
library(patchwork)

source("code/01-lcbd-models.R")


# 2. Functions -------------------------------------------------------------

xprime_to_c <- function(xprime) (1 - xprime) / xprime

number_label <- function(x) {
    sub("[.]?0+$", "", formatC(x, digits = 6, format = "f"))
}

regional_richness <- function(N, x) {
    as.integer(round(-N * (1 - x) / x * log1p(-x)))
}

richness_summary <- function(S, M, x, xprime) {
    y <- neutral_y(M, xprime_to_c(xprime))
    q <- 1 - log1p(-x / y) / log1p(-x)
    probability <- stats::dbinom(0:S, S, q)
    probability[1L] <- 0
    probability <- probability / sum(probability)
    cumulative <- cumsum(probability)
    richness <- 0:S

    c(
        typical = sum(richness * probability),
        lower = richness[which(cumulative >= 0.025)[1L]],
        upper = richness[which(cumulative >= 0.975)[1L]]
    )
}

parameterized_curve <- function(framework, S, M, x, richness) {
    if (framework == "Coverage deficit") {
        return(coverage_parameterized(S, M, x, richness))
    }
    pmax(0, sorensen_mean_other(
        sorensen_parameterized(S, M, x, richness), M
    ))
}

annealed_curve <- function(framework, S, M, x, c, richness) {
    if (framework == "Coverage deficit") {
        return(coverage_annealed(S, M, x, c, richness))
    }
    sorensen_mean_other(sorensen_annealed(S, M, x, c, richness), M)
}

quenched_curve <- function(framework, sad, M, c) {
    if (framework == "Coverage deficit") {
        return(coverage_quenched(sad, M, c)[, c("alpha", "Quenched")])
    }
    result <- sorensen_quenched_log_domain(
        sad, M, c, alpha = 0:length(sad)
    )
    data.frame(alpha = result$alpha, Quenched = result$expected_Sorensen)
}


# 3. Regional species pools and representative SAD ------------------------

regional_settings <- data.frame(
    x = x_values,
    theta = N * (1 - x_values) / x_values,
    S = vapply(x_values, regional_richness, integer(1), N = N)
)
middle <- regional_settings[regional_settings$x == conditional_x, ]

sad_input <- readRDS(sad_file)
middle_index <- which.min(abs(sad_input$settings$x_grid - conditional_x))
candidate_sads <- sad_input$sad_ensembles[[middle_index]]
candidate_totals <- vapply(candidate_sads, sum, numeric(1))
representative_sad <- candidate_sads[[which.min(abs(candidate_totals - N))]]


# 4. Calculate curves and typical-richness points -------------------------

frameworks <- c("Coverage deficit", "Sørensen")
x_labels <- number_label(x_values)
xprime_labels <- number_label(xprime_values)
curve_rows <- list()
point_rows <- list()
range_rows <- list()

for (framework in frameworks) {
    for (i in seq_along(x_values)) {
        x <- x_values[i]
        S <- regional_settings$S[i]
        richness <- exp(seq(log(1), log(S), length.out = 500L))

        curve_rows[[length(curve_rows) + 1L]] <- data.frame(
            framework = framework,
            richness = richness,
            expected = parameterized_curve(framework, S, M, x, richness),
            baseline = "Parameterized",
            x = x_labels[i],
            xprime = NA_character_
        )

        for (j in seq_along(xprime_values)) {
            summary <- richness_summary(S, M, x, xprime_values[j])
            point_rows[[length(point_rows) + 1L]] <- data.frame(
                framework = framework,
                richness = unname(summary["typical"]),
                expected = parameterized_curve(
                    framework, S, M, x, unname(summary["typical"])
                ),
                baseline = "Parameterized",
                x = x_labels[i],
                xprime = xprime_labels[j]
            )
        }
    }

    richness <- seq_len(middle$S)
    for (j in seq_along(xprime_values)) {
        xprime <- xprime_values[j]
        c_value <- xprime_to_c(xprime)
        summary <- richness_summary(middle$S, M, conditional_x, xprime)
        annealed <- annealed_curve(
            framework, middle$S, M, conditional_x, c_value, richness
        )
        quenched <- quenched_curve(framework, representative_sad, M, c_value)
        quenched <- quenched$Quenched[match(richness, quenched$alpha)]

        curve_rows[[length(curve_rows) + 1L]] <- data.frame(
            framework = framework,
            richness = rep(richness, 2L),
            expected = c(quenched, annealed),
            baseline = rep(c("Quenched", "Annealed"), each = length(richness)),
            x = number_label(conditional_x),
            xprime = xprime_labels[j]
        )
        range_rows[[length(range_rows) + 1L]] <- data.frame(
            framework = framework,
            xprime = xprime_labels[j],
            typical = unname(summary["typical"]),
            lower = unname(summary["lower"]),
            upper = unname(summary["upper"])
        )
    }
}

curves <- do.call(rbind, curve_rows)
curves <- curves[is.finite(curves$expected), ]
parameterized_points <- do.call(rbind, point_rows)
typical_ranges <- unique(do.call(rbind, range_rows))


# 5. Plot -----------------------------------------------------------------

manuscript_theme <- theme_classic(
    base_size = figure_style$text,
    base_family = font_family
) +
    theme(
        panel.border = element_rect(
            fill = NA, linewidth = figure_style$border
        ),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = figure_style$border),
        axis.title = element_text(size = figure_style$text),
        axis.text = element_text(size = figure_style$axis_text),
        legend.title = element_text(size = figure_style$legend),
        legend.text = element_text(size = figure_style$legend),
        plot.title = element_text(
            size = figure_style$text, face = "bold", hjust = 0.5,
            margin = margin(b = 2)
        ),
        plot.margin = margin(2, 8, 0, 6)
    )

# Draw one three-baseline curve panel.
make_panel <- function(framework, tag, show_y, legend_type) {
    data <- curves[curves$framework == framework, ]
    parameterized <- data[data$baseline == "Parameterized", ]
    conditional <- data[data$baseline != "Parameterized", ]
    bands <- typical_ranges[typical_ranges$framework == framework, ]
    index <- match(conditional$xprime, bands$xprime)
    conditional$lower <- bands$lower[index]
    conditional$upper <- bands$upper[index]
    highlighted <- conditional[
        conditional$richness >= conditional$lower &
            conditional$richness <= conditional$upper,
    ]

    conditional$line_colour <- xprime_colours[conditional$xprime]
    highlighted$line_colour <- xprime_colours[highlighted$xprime]
    parameter_points <- parameterized_points[
        parameterized_points$framework == framework,
    ]
    parameter_points$point_colour <- xprime_marker_colours[
        parameter_points$xprime
    ]

    conditional_points <- do.call(rbind, lapply(seq_len(nrow(bands)), function(i) {
        rows <- conditional[conditional$xprime == bands$xprime[i], ]
        do.call(rbind, lapply(c("Quenched", "Annealed"), function(model) {
            curve <- rows[rows$baseline == model, ]
            data.frame(
                richness = bands$typical[i],
                expected = stats::approx(
                    curve$richness, curve$expected, xout = bands$typical[i]
                )$y,
                baseline = model,
                xprime = bands$xprime[i],
                point_colour = xprime_marker_colours[bands$xprime[i]]
            )
        }))
    }))
    legend_points <- data.frame(
        richness = NA_real_,
        expected = NA_real_,
        line_colour = unname(xprime_colours),
        point_colour = unname(xprime_colours)
    )

    label_targets <- data.frame(
        x = x_labels,
        target = c(350, 15, 3),
        label_x = c(285, 12, 4.2),
        segment_x = c(310, 13, 3.8),
        hjust = c(1, 1, 0)
    )
    labels <- do.call(rbind, lapply(seq_len(nrow(label_targets)), function(i) {
        curve <- parameterized[parameterized$x == label_targets$x[i], ]
        target_y <- stats::approx(
            curve$richness, curve$expected, xout = label_targets$target[i]
        )$y
        data.frame(
            x = label_targets$x[i],
            target = label_targets$target[i],
            target_y = target_y,
            label_x = label_targets$label_x[i],
            segment_x = label_targets$segment_x[i],
            hjust = label_targets$hjust[i],
            label_y = target_y
        )
    }))
    labels$label <- paste0("x = ", labels$x)
    labels$label_colour <- ifelse(
        labels$x == number_label(conditional_x), "#000000", "#777777"
    )

    ggplot() +
        geom_line(
            data = parameterized,
            aes(richness, expected, linetype = baseline, group = x),
            colour = "#858585", linewidth = figure_style$line * 0.62
        ) +
        geom_line(
            data = conditional,
            aes(
                richness, expected, colour = line_colour,
                linetype = baseline, group = interaction(xprime, baseline)
            ),
            linewidth = figure_style$line * 0.48, alpha = 0.72,
            na.rm = TRUE
        ) +
        geom_line(
            data = highlighted[highlighted$baseline == "Quenched", ],
            aes(
                richness, expected, colour = line_colour,
                linetype = baseline, group = xprime
            ),
            linewidth = figure_style$line * 1.10, na.rm = TRUE
        ) +
        geom_line(
            data = highlighted[highlighted$baseline == "Annealed", ],
            aes(
                richness, expected, colour = line_colour,
                linetype = baseline, group = xprime
            ),
            linewidth = figure_style$line * 1.45, na.rm = TRUE
        ) +
        geom_point(
            data = parameter_points,
            aes(richness, expected, fill = point_colour),
            shape = 23, colour = "white", stroke = 0.75,
            size = figure_style$diamond, show.legend = FALSE
        ) +
        geom_point(
            data = conditional_points,
            aes(richness, expected, fill = point_colour),
            shape = 23, colour = "white", stroke = 0.75,
            size = figure_style$diamond, show.legend = FALSE
        ) +
        geom_point(
            data = transform(
                legend_points,
                richness = NA_real_, expected = NA_real_
            ),
            aes(
                richness, expected,
                colour = line_colour
            ),
            shape = 23,
            size = figure_style$diamond,
            stroke = 0.75,
            key_glyph = draw_key_xprime_point,
            inherit.aes = FALSE,
            show.legend = legend_type == "colour",
            na.rm = TRUE
        ) +
        geom_segment(
            data = labels,
            aes(
                x = segment_x, y = label_y,
                xend = target, yend = target_y,
                colour = label_colour
            ),
            linewidth = figure_style$line * 0.72,
            show.legend = FALSE
        ) +
        geom_text(
            data = labels,
            aes(label_x, label_y, label = label, colour = label_colour),
            hjust = labels$hjust, vjust = 0.5,
            size = (figure_style$legend - 1) / ggplot2::.pt,
            show.legend = FALSE
        ) +
        annotate(
            "text", x = Inf, y = Inf, label = tag,
            hjust = 1.45, vjust = 1.35, fontface = "bold",
            size = figure_style$tag / ggplot2::.pt
        ) +
        scale_x_log10(
            breaks = x_ticks, limits = c(1, 1600),
            labels = scales::label_comma(),
            expand = expansion(mult = c(0.02, 0.03))
        ) +
        scale_y_continuous(
            breaks = seq(0, 1, 0.2),
            labels = if (show_y) waiver() else NULL,
            limits = c(0, 1.01), expand = expansion(mult = c(0.01, 0.01))
        ) +
        scale_colour_identity(
            name = expression(italic(x)*"'"),
            breaks = unname(xprime_colours), labels = names(xprime_colours),
            guide = if (legend_type == "colour") "legend" else "none"
        ) +
        scale_fill_identity(guide = "none") +
        scale_linetype_manual(
            name = NULL,
            values = c(
                Parameterized = "dotted",
                Quenched = "longdash",
                Annealed = "solid"
            ),
            breaks = c("Parameterized", "Quenched", "Annealed"),
            guide = if (legend_type == "linetype") "legend" else "none"
        ) +
        guides(
            linetype = if (legend_type == "linetype") guide_legend(
                override.aes = list(
                    colour = c("#858585", "#383838", "#383838"),
                    linewidth = c(
                        figure_style$line * 0.62,
                        figure_style$line,
                        figure_style$line
                    )
                )
            ) else "none",
            colour = if (legend_type == "colour") guide_legend(
                title.position = "left",
                override.aes = list(
                    linetype = "solid",
                    linewidth = figure_style$line
                )
            ) else "none"
        ) +
        labs(
            x = "Local richness",
            y = if (show_y) expression(LCBD~"\u00d7"~M) else NULL,
            title = framework
        ) +
        manuscript_theme +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.96, 0.90),
            legend.justification = c(1, 1),
            legend.box = "horizontal",
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.key.width = grid::unit(2.2, "lines"),
            legend.key.height = grid::unit(0.62, "lines"),
            legend.title = if (legend_type == "colour") {
                element_text(
                    size = figure_style$legend, vjust = 0.5,
                    margin = margin(r = 4)
                )
            } else {
                element_text(size = figure_style$legend)
            },
            legend.margin = margin(0),
            aspect.ratio = 0.75
        )
}

panel_a <- make_panel("Coverage deficit", "a", TRUE, "linetype")
panel_b <- make_panel("Sørensen", "b", FALSE, "colour")
curve_panel_a <- panel_a
curve_panel_b <- panel_b
figure_plot <- panel_a + panel_b + plot_layout(widths = c(1, 1)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))


# 6. Parameter-space results ---------------------------------------------

results <- list(
    coverage = readRDS(coverage_result_file),
    sorensen = readRDS(sorensen_result_file)
)

parameter_label <- function(x) {
    sub("[.]?0+$", "", formatC(x, digits = 6, format = "f"))
}

cell_edges <- function(x) {
    x <- sort(unique(x))
    c(
        x[1L] - diff(x[1:2]) / 2,
        (x[-1L] + x[-length(x)]) / 2,
        x[length(x)] + diff(tail(x, 2L)) / 2
    )
}

# Convert parameter-space output to the two plotted relative differences.
prepare_parameter_data <- function(result, framework) {
    summary <- result$slope_summary
    summary$annealed_parameterized <- ifelse(
        abs(summary$slope_model1_mean) > 1e-10,
        100 * (summary$slope_model3_mean -
            summary$slope_model1_mean) /
            abs(summary$slope_model1_mean),
        NA_real_
    )

    replicate <- result$replicate_slopes
    replicate$annealed_quenched <- ifelse(
        abs(replicate$slope_model3_pair_at_alpha2) > 1e-10,
        100 * abs(
            replicate$slope_model3_pair_at_alpha2 -
                replicate$slope_model2_at_alpha2
        ) / abs(replicate$slope_model3_pair_at_alpha2),
        NA_real_
    )
    quenched_summary <- aggregate(
        annealed_quenched ~ x + x_prime,
        data = replicate,
        FUN = stats::median
    )

    data <- merge(
        summary,
        quenched_summary,
        by = c("x", "x_prime"),
        sort = FALSE
    )
    data$framework <- framework
    data
}

parameter_data <- rbind(
    prepare_parameter_data(results$coverage, "Coverage deficit"),
    prepare_parameter_data(results$sorensen, "Sørensen")
)

parameter_data$x_plot <- -log10(1 - parameter_data$x)
parameter_data$xprime_plot <- -log10(1 - parameter_data$x_prime)

x_values <- sort(unique(parameter_data$x_plot))
xprime_values <- sort(unique(parameter_data$xprime_plot))
x_edges <- cell_edges(x_values)
xprime_edges <- cell_edges(xprime_values)
xprime_display_lower <- 0

parameter_data$xmin <- x_edges[
    match(parameter_data$x_plot, x_values)
]
parameter_data$xmax <- x_edges[
    match(parameter_data$x_plot, x_values) + 1L
]
parameter_data$ymin <- xprime_edges[
    match(parameter_data$xprime_plot, xprime_values)
]
parameter_data$ymax <- xprime_edges[
    match(parameter_data$xprime_plot, xprime_values) + 1L
]

shared_breaks <- seq(0, 120, by = 20)
shared_limits <- c(0, 120)

x_anchor <- results$coverage$settings$x_anchor
x_prime_anchor <- results$coverage$settings$x_prime_anchor
x_breaks <- -log10(1 - x_anchor)
x_break_labels <- parameter_label(x_anchor)
x_prime_ticks <- x_prime_anchor[x_prime_anchor >= 0.5]
xprime_breaks <- -log10(1 - x_prime_ticks)
xprime_break_labels <- parameter_label(x_prime_ticks)

marker_data <- expand.grid(
    framework = c("Coverage deficit", "Sørensen"),
    xprime = c(0.5, 0.99),
    stringsAsFactors = FALSE
)
marker_data$x_plot <- -log10(1 - 0.9999)
marker_data$xprime_plot <- -log10(1 - marker_data$xprime)
marker_data$colour <- unname(
    xprime_marker_colours[parameter_label(marker_data$xprime)]
)


# 7. Parameter-space panels ----------------------------------------------

manuscript_theme <- theme_classic(
    base_size = figure_style$text,
    base_family = font_family
) +
    theme(
        panel.border = element_rect(
            fill = NA,
            linewidth = figure_style$border
        ),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = figure_style$border),
        axis.title = element_text(size = figure_style$text),
        axis.text = element_text(size = figure_style$axis_text),
        legend.title = element_text(size = figure_style$legend),
        legend.text = element_text(size = figure_style$legend),
        plot.margin = margin(2, 8, 0, 6)
    )

# Draw one parameter-space comparison panel.
make_heatmap <- function(
        framework, field, tag, show_y, show_x, show_legend
) {
    data <- parameter_data[parameter_data$framework == framework, ]
    data$value <- data[[field]]
    markers <- marker_data[marker_data$framework == framework, ]

    ggplot(data) +
        geom_rect(
            aes(
                xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = value
            ),
            colour = NA
        ) +
        geom_abline(
            slope = 1,
            intercept = 0,
            colour = "#4D4D4D",
            linewidth = 0.55
        ) +
        geom_point(
            data = markers,
            aes(x_plot, xprime_plot),
            shape = 21,
            size = figure_style$point * 1.55,
            stroke = 0.9,
            colour = "white",
            fill = markers$colour,
            inherit.aes = FALSE,
            show.legend = FALSE
        ) +
        annotate(
            "text",
            x = Inf,
            y = Inf,
            label = tag,
            hjust = 1.45,
            vjust = 1.35,
            fontface = "bold",
            size = figure_style$tag / ggplot2::.pt
        ) +
        scale_x_continuous(
            breaks = x_breaks,
            labels = x_break_labels,
            expand = expansion(mult = 0)
        ) +
        scale_y_continuous(
            breaks = xprime_breaks,
            labels = if (show_y) xprime_break_labels else NULL,
            expand = expansion(mult = 0)
        ) +
        scale_fill_gradient2(
            low = heatmap_colours[1L],
            mid = heatmap_colours[2L],
            high = heatmap_colours[3L],
            midpoint = 0,
            limits = shared_limits,
            breaks = shared_breaks,
            oob = scales::squish,
            name = "Relative difference (%)"
        ) +
        coord_cartesian(
            xlim = range(x_edges),
            ylim = c(xprime_display_lower, tail(xprime_edges, 1L)),
            expand = FALSE
        ) +
        labs(
            x = if (show_x) expression(italic(x)) else NULL,
            y = if (show_y) expression(italic(x)*"'") else NULL
        ) +
        manuscript_theme +
        theme(
            axis.text.x = if (show_x) {
                element_text(angle = 25, hjust = 1, vjust = 1)
            } else {
                element_blank()
            },
            axis.ticks.x = if (show_x) {
                element_line(linewidth = figure_style$border)
            } else {
                element_blank()
            },
            legend.position = if (show_legend) "bottom" else "none",
            legend.key.width = grid::unit(5.5, "lines"),
            legend.key.height = grid::unit(0.55, "lines"),
            legend.margin = margin(-4, 0, 0, 0),
            axis.title.y = if (show_y) {
                element_text(margin = margin(r = -12))
            } else {
                element_blank()
            },
            aspect.ratio = 0.75
        ) +
        guides(fill = guide_colourbar(
            title.position = "top",
            title.hjust = 0.5,
            ticks.colour = "black",
            frame.colour = "black",
            barwidth = grid::unit(12, "lines"),
            barheight = grid::unit(0.55, "lines")
        ))
}

panel_c <- make_heatmap(
    "Coverage deficit", "annealed_parameterized",
    "c", TRUE, TRUE, TRUE
)
panel_d <- make_heatmap(
    "Sørensen", "annealed_parameterized",
    "d", FALSE, TRUE, TRUE
)
panel_e <- make_heatmap(
    "Coverage deficit", "annealed_quenched",
    "e", TRUE, TRUE, TRUE
)
panel_f <- make_heatmap(
    "Sørensen", "annealed_quenched",
    "f", FALSE, TRUE, TRUE
)


# 8. Assemble and export --------------------------------------------------

shared_legend <- cowplot::get_legend(
    panel_c + theme(legend.position = "bottom")
)
panel_c <- panel_c + theme(legend.position = "none")
panel_d <- panel_d + theme(legend.position = "none")
panel_e <- panel_e + theme(legend.position = "none")
panel_f <- panel_f + theme(legend.position = "none")

legend_panel <- cowplot::ggdraw() +
    cowplot::draw_grob(shared_legend, x = 0.005, width = 1)

curve_row <- curve_panel_a + curve_panel_b +
    plot_layout(widths = c(1, 1), guides = "keep")
parameterized_row <- panel_c + panel_d +
    plot_layout(widths = c(1, 1), guides = "keep")
quenched_row <- panel_e + panel_f +
    plot_layout(widths = c(1, 1), guides = "keep")

figure_plot <- curve_row / parameterized_row / quenched_row / legend_panel +
    plot_layout(heights = c(1, 1, 1, 0.13)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))

ragg::agg_png(
    paste0(output_base, ".png"),
    width = figure_size$width,
    height = figure_size$height,
    units = "in",
    res = figure_size$dpi,
    background = "white"
)
print(figure_plot)
grDevices::dev.off()

ragg::agg_tiff(
    paste0(output_base, ".tiff"),
    width = figure_size$width,
    height = figure_size$height,
    units = "in",
    res = 600,
    background = "white",
    compression = "lzw"
)
print(figure_plot)
grDevices::dev.off()

grDevices::pdf(
    paste0(output_base, ".pdf"),
    width = figure_size$width,
    height = figure_size$height,
    family = "Helvetica",
    useDingbats = FALSE,
    bg = "white"
)
print(figure_plot)
grDevices::dev.off()
