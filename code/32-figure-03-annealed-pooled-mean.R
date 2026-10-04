# Plot the pooled-mean annealed validation used in Figure 3.
# Run after script 31 has completed both Figure 3 simulations.

# 1. Settings -------------------------------------------------------------

S <- 300L
M <- 1250L
x <- 0.9999
x_prime_values <- c(0.5, 0.99)
seed <- 20260922L
min_errorbar_replicates <- 50L

point_colours <- c("#756BB1", "#5A8F79")
point_shapes <- c(23, 24)
figure_width <- 8.1
figure_height <- 4.6
font_family <- "Helvetica"

# 2. Input and output paths ----------------------------------------------

number_tag <- function(value) {
    gsub("[^0-9A-Za-z]+", "p", format(value, digits = 10, scientific = FALSE))
}

run_id <- function(x_prime) {
    paste0(
        "annealed_S", S, "_M", M, "_x", number_tag(x),
        "_xp", number_tag(x_prime), "_seed", seed
    )
}

result_dir <- "process data/Figure3_annealed_validation"
result_files <- file.path(
    result_dir, paste0(vapply(x_prime_values, run_id, character(1)), ".rds")
)
formal_figure_base <- paste0(
    "figure/Figure3_simulation_annealed_",
    "pooled_mean_1.96pooledSD"
)

# 3. Read the completed runs ---------------------------------------------

results <- lapply(result_files, readRDS)

# 4. Summarize all matrices by integer richness ---------------------------

frameworks <- c("Coverage", "Sørensen")

for (i in seq_along(results)) {
    summaries <- results[[i]]$replicate_summaries
    count <- summaries$count
    sample_size <- colSums(count)
    contributors <- colSums(count > 0L)
    tables <- vector("list", length(frameworks))

    for (j in seq_along(frameworks)) {
        framework <- frameworks[j]
        stored <- results[[i]]$point_results
        stored <- stored[stored$framework == framework, ]
        tables[[j]] <- data.frame(
            framework = framework,
            richness = 0:S,
            sample_size = sample_size,
            contributing_replicates = contributors,
            pooled_mean = stored$mc_mean,
            pooled_sd = stored$mc_sd,
            analytic = stored$analytic
        )
    }
    results[[i]]$point_results <- do.call(rbind, tables)
}

# 5. Two LCBD frameworks in one figure ----------------------------------

x_ticks <- sort(unique(c(1, 2, 5, 10, 20, 50, 100, S)))
x_ticks <- x_ticks[x_ticks <= S]
# Draw the pooled conditional means for both LCBD frameworks.
draw_figure <- function() {
all_values <- unlist(lapply(results, function(result) {
    z <- result$point_results
    z <- z[z$richness > 0L, ]
    supported <- z$contributing_replicates >= min_errorbar_replicates
    c(
        z$analytic, z$pooled_mean,
        (z$pooled_mean - 1.96 * z$pooled_sd)[supported],
        (z$pooled_mean + 1.96 * z$pooled_sd)[supported]
    )
}))
common_y_limits <- extendrange(all_values[is.finite(all_values)], f = 0.03)
common_y_limits[1L] <- max(0, common_y_limits[1L])

par(
    mfrow = c(1, 2), mar = c(3.3, 2.5, 2.6, 0.4),
    oma = c(1.00, 1.75, 0.2, 0), family = font_family,
    mgp = c(1.45, 0.30, 0), tcl = -0.20, las = 1
)

for (framework in frameworks) {
    if (framework == "Coverage") {
        par(mar = c(2.35, 2.10, 2.6, 0.2))
    } else {
        par(mar = c(2.35, 2.10, 2.6, 0.5))
    }
    plot(
        1, 1, type = "n", log = "x", axes = FALSE,
        xlim = c(1, S), ylim = common_y_limits,
        xlab = "", ylab = "", xaxs = "i", yaxs = "i"
    )
    axis(1, at = x_ticks, labels = FALSE, lwd = 0, lwd.ticks = 0.75)
    axis(1, at = x_ticks[-c(1L, length(x_ticks))],
        labels = x_ticks[-c(1L, length(x_ticks))],
        cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
    axis(1, at = x_ticks[1L], labels = x_ticks[1L], hadj = 0,
        cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
    axis(1, at = x_ticks[length(x_ticks)],
        labels = x_ticks[length(x_ticks)], hadj = 1,
        cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
    axis(
        2, at = seq(0, 1, 0.2), labels = format(seq(0, 1, 0.2), nsmall = 1),
        las = 1, cex.axis = 0.88, lwd = 0, lwd.ticks = 0.75
    )
    box(col = "grey20", lwd = 0.85)
    title(
        main = if (framework == "Coverage") "Coverage deficit" else framework,
        cex.main = 0.92, line = 0.25
    )
    mtext(
        if (framework == "Coverage") "a" else "b",
        side = 3, line = -1.35, adj = 0.97,
        cex = 1.08, font = 2
    )

    for (i in seq_along(results)) {
        z <- results[[i]]$point_results
        z <- z[z$framework == framework & z$richness > 0L, ]
        centre <- z$pooled_mean
        lower <- z$pooled_mean - 1.96 * z$pooled_sd
        upper <- z$pooled_mean + 1.96 * z$pooled_sd
        shown_points <-
            z$contributing_replicates > 0L & is.finite(centre)
        shown_intervals <- shown_points &
            z$contributing_replicates >= min_errorbar_replicates &
            is.finite(lower) & is.finite(upper) & upper > lower
        sparse_points <- shown_points &
            z$contributing_replicates < min_errorbar_replicates
        lines(
            z$richness, z$analytic,
            col = "black", lty = 1, lwd = 1.25
        )
        error_colour <- adjustcolor(point_colours[i], 0.72)
        segments(
            z$richness[shown_intervals], lower[shown_intervals],
            z$richness[shown_intervals], upper[shown_intervals],
            col = error_colour, lwd = 1.0
        )
        segments(
            z$richness[shown_intervals] / 1.025,
            lower[shown_intervals],
            z$richness[shown_intervals] * 1.025,
            lower[shown_intervals],
            col = error_colour, lwd = 1.0
        )
        segments(
            z$richness[shown_intervals] / 1.025,
            upper[shown_intervals],
            z$richness[shown_intervals] * 1.025,
            upper[shown_intervals],
            col = error_colour, lwd = 1.0
        )
        points(
            z$richness[shown_intervals], centre[shown_intervals],
            pch = point_shapes[i], col = "white",
            bg = point_colours[i], cex = 0.60, lwd = 0.55
        )
        points(
            z$richness[sparse_points], centre[sparse_points],
            pch = 21, col = "white", bg = "grey65",
            cex = 0.60, lwd = 0.55
        )
    }

    if (framework == "Coverage") {
        legend(
            "bottomleft",
            c(
                "Analytical annealed model",
                paste0("x' = ", trimws(formatC(
                    x_prime_values, format = "fg", digits = 6,
                    drop0trailing = TRUE
                ))),
                "Rarely realized richness"
            ),
            col = c(
                "black", rep("white", length(x_prime_values)), "white"
            ),
            pt.bg = c(NA, point_colours, "grey65"),
            lty = c(1, rep(NA, length(x_prime_values) + 1L)), lwd = 1.25,
            pch = c(NA, point_shapes, 21), pt.cex = 0.8,
            bty = "n", cex = 0.75, seg.len = 1.35, y.intersp = 0.90
        )
    }
}

mtext("Local richness", side = 1, outer = TRUE, line = 0.10, cex = 1.12)
mtext(
    "LCBD × M", side = 2, outer = TRUE,
    line = 0.25, cex = 1.14, las = 0
)
}

# Export Figure 3 in the manuscript formats.
export_figure <- function(output_base) {
    ragg::agg_png(
        paste0(output_base, ".png"), width = figure_width,
        height = figure_height, units = "in", res = 600,
        pointsize = 13, background = "white"
    )
    draw_figure()
    dev.off()

    ragg::agg_tiff(
        paste0(output_base, ".tiff"), width = figure_width,
        height = figure_height, units = "in", res = 600,
        pointsize = 13, background = "white", compression = "lzw"
    )
    draw_figure()
    dev.off()

    grDevices::pdf(
        paste0(output_base, ".pdf"), width = figure_width,
        height = figure_height, pointsize = 13,
        family = "Helvetica", useDingbats = FALSE, bg = "white"
    )
    draw_figure()
    dev.off()

}

export_figure(formal_figure_base)
