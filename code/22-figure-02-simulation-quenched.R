# Plot the quenched simulation validation prepared by script 21.

# 1. Paths and settings ---------------------------------------------------

model2_file <- "process data/Figure2_quenched_expectations.rds"
figure2_file <- "figure/Figure2_simulation_quenched.png"
panels <- c("implicit", "explicit")
frameworks <- c("coverage", "sorensen")
grains <- c("5", "20", "100")
grain_labels <- c(
    "5" = "5 m × 5 m",
    "20" = "20 m × 20 m",
    "100" = "100 m × 100 m"
)
grain_colours <- c("5" = "#3378A6", "20" = "#5D8B73", "100" = "#B27645")
residual_colour <- "#5D8B73"
font_family <- "Helvetica"
richness_coordinate <- function(x) log10(x + 1)

# 2. Read fitted objects -------------------------------------------------

model2_source <- readRDS(model2_file)$frameworks

# 3. Keep only the simulation panels ------------------------------------

# Combine the saved framework results into one plotting object.
assemble_result <- function() {
    observed <- fits <- list()
    for (framework in frameworks) {
        observed[[framework]] <- fits[[framework]] <- list()
        source <- model2_source[[framework]]
        for (panel in panels) {
            observed[[framework]][[panel]] <- source$observed[[panel]]
            fits[[framework]][[panel]] <- list()
            for (grain in grains) {
                fitted <- source$fits[[panel]][[grain]]
                fits[[framework]][[panel]][[grain]] <- list(
                    selected = fitted$model2,
                    known = fitted$known_parameter,
                    parameterized = fitted$fitted_x
                )
            }
        }
    }
    list(model = "Model 2", observed = observed, fits = fits)
}

model2_result <- assemble_result()

# 4. Prepare residuals and regressions -----------------------------------

# Calculate LCBD deviations and their richness regressions.
prepare_residuals <- function(result) {
    output <- diagnostics <- list()
    k <- 0L
    for (framework in frameworks) {
        output[[framework]] <- list()
        for (panel in panels) {
            for (grain in grains) {
                data <- result$observed[[framework]][[panel]][[grain]]
                curve <- result$fits[[framework]][[panel]][[grain]]$selected$curve
                row <- match(data$richness, curve$richness)
                data$expected <- curve$expected[row]
                data$residual <- data$lcbd - data$expected
                data <- data[data$richness > 0, , drop = FALSE]
                linear <- lm(residual ~ log(richness), data = data)
                k <- k + 1L
                diagnostics[[k]] <- data.frame(
                    model = result$model, framework = framework,
                    simulation = panel, grain = grain, n = nrow(data),
                    bias = mean(data$residual),
                    rmse = sqrt(mean(data$residual^2)),
                    slope = unname(coef(linear)[2L]),
                    slope_se = summary(linear)$coefficients[2L, 2L],
                    p = summary(linear)$coefficients[2L, 4L],
                    adjusted_r2 = summary(linear)$adj.r.squared
                )
                if (grain != "20") next

            richness <- exp(seq(
                log(min(data$richness)), log(max(data$richness)),
                length.out = 300L
            ))
            prediction <- predict(
                linear, data.frame(richness = richness), interval = "confidence"
            )
            output[[framework]][[panel]] <- list(
                data = data,
                line = data.frame(
                    richness = richness, fitted = prediction[, "fit"],
                    lower = prediction[, "lwr"], upper = prediction[, "upr"]
                ),
                significant = summary(linear)$coefficients[2L, 4L] < 0.05
            )
            }
        }
    }
    diagnostics <- do.call(rbind, diagnostics)
    diagnostics$significant <- diagnostics$p < 0.05
    list(
        panels = output,
        diagnostics = diagnostics
    )
}

model2_result$residuals <- prepare_residuals(model2_result)

# 5. Hexagon and curve functions -----------------------------------------

make_palette <- function(colour, opaque = FALSE) {
    palette <- colorRampPalette(c(
        adjustcolor(colour, alpha.f = 0.16),
        adjustcolor(colour, alpha.f = 0.42),
        adjustcolor(colour, alpha.f = 0.72), colour
    ), alpha = TRUE)(256L)
    if (!opaque) return(palette)
    rgba <- col2rgb(palette, alpha = TRUE) / 255
    rgb(
        rgba[1L, ] * rgba[4L, ] + 1 - rgba[4L, ],
        rgba[2L, ] * rgba[4L, ] + 1 - rgba[4L, ],
        rgba[3L, ] * rgba[4L, ] + 1 - rgba[4L, ]
    )
}

grain_palettes <- lapply(grain_colours, make_palette, opaque = TRUE)
residual_palette <- make_palette(residual_colour)

# Aggregate dense observations into hexagonal plotting cells.
make_hexagons <- function(x, y, x_range, y_range, nx = 30, ny = 20) {
    keep <- is.finite(x) & is.finite(y) &
        x >= x_range[1L] & x <= x_range[2L] &
        y >= y_range[1L] & y <= y_range[2L]
    x <- x[keep]
    y <- y[keep]
    panel_inches <- par("pin")
    x_scale <- panel_inches[1L] / diff(x_range)
    y_scale <- panel_inches[2L] / diff(y_range)
    cellsize <- min(panel_inches[1L] / nx, panel_inches[2L] / ny)
    extent <- sf::st_as_sfc(sf::st_bbox(c(
        xmin = -cellsize, ymin = -cellsize,
        xmax = panel_inches[1L] + cellsize,
        ymax = panel_inches[2L] + cellsize
    )))
    grid <- sf::st_make_grid(extent, cellsize = cellsize, square = FALSE)
    points <- data.frame(
        x = (x - x_range[1L]) * x_scale,
        y = (y - y_range[1L]) * y_scale
    )
    membership <- vapply(
        sf::st_intersects(sf::st_as_sf(points, coords = c("x", "y")), grid),
        min, integer(1)
    )
    counts <- tabulate(membership, nbins = length(grid))
    occupied <- which(counts > 0L)
    vertices <- do.call(rbind, lapply(occupied, function(cell) {
        xy <- sf::st_coordinates(grid[cell])[, 1:2, drop = FALSE]
        data.frame(
            cell = cell,
            x = xy[, 1L] / x_scale + x_range[1L],
            y = xy[, 2L] / y_scale + y_range[1L]
        )
    }))
    polygons <- sf::st_sfc(lapply(
        split(vertices[, c("x", "y")], vertices$cell),
        function(z) sf::st_polygon(list(as.matrix(z)))
    ))
    list(
        vertices = vertices, counts = counts, maximum = max(counts),
        mask = sf::st_union(polygons)
    )
}

draw_hexagons <- function(hexagons, palette) {
    for (cell in unique(hexagons$vertices$cell)) {
        vertices <- hexagons$vertices[hexagons$vertices$cell == cell, ]
        fraction <- (log1p(hexagons$counts[cell]) /
            log1p(hexagons$maximum))^2.15
        polygon(
            vertices$x, vertices$y,
            col = palette[1L + round(255 * fraction)], border = NA
        )
    }
}

smooth_curve <- function(curve, n = 1200L) {
    curve <- curve[curve$richness >= 1 & is.finite(curve$expected), ]
    x <- richness_coordinate(curve$richness)
    x_draw <- seq(min(x), max(x), length.out = n)
    data.frame(
        x = x_draw,
        expected = splinefun(x, curve$expected, method = "monoH.FC")(x_draw)
    )
}

draw_curve <- function(
        curve, colour, line_type, line_width = 1.25, hexagons = NULL,
        covering_hexagons = NULL
) {
    curve <- smooth_curve(curve)
    if (!is.null(hexagons)) {
        curve_sf <- sf::st_as_sf(curve, coords = c("x", "expected"))
        visible <- lengths(sf::st_intersects(curve_sf, hexagons$mask)) > 0L
        for (other in covering_hexagons) {
            visible <- visible &
                lengths(sf::st_intersects(curve_sf, other$mask)) == 0L
        }
        curve$expected[!visible] <- NA_real_
    }
    lines(
        curve$x, curve$expected, col = colour,
        lwd = line_width, lty = line_type
    )
}

richness_ticks <- function(S) {
    ticks <- sort(unique(c(1, 10, 50, 100, 200, S)))
    ticks[ticks <= S]
}

draw_x_axis <- function(ticks, labels, inset_ends = FALSE) {
    position <- richness_coordinate(ticks)
    axis(1, at = position, labels = FALSE, lwd = 0, lwd.ticks = 0.75)
    if (!labels) return(invisible(NULL))
    if (!inset_ends) {
        axis(1, at = position, labels = ticks,
            cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
        return(invisible(NULL))
    }
    inset <- 0.018 * diff(range(position))
    if (length(ticks) > 2L) {
        axis(1, at = position[-c(1L, length(ticks))],
            labels = ticks[-c(1L, length(ticks))],
            cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
    }
    axis(1, at = position[1L] + inset, labels = ticks[1L], hadj = 0,
        cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
    axis(1, at = position[length(position)] - inset,
        labels = ticks[length(ticks)], hadj = 1,
        cex.axis = 0.90, lwd = 0, lwd.ticks = 0)
}

draw_title_cell <- function(label, rotation = 0, cex = 1.0, font = 1, y = 0.5) {
    par(mar = rep(0, 4), family = font_family)
    plot.new()
    text(0.5, y, label, srt = rotation, cex = cex, font = font)
}

draw_tag <- function(tag) {
    usr <- par("usr")
    text(
        usr[2L] - 0.035 * diff(usr[1:2]),
        usr[4L] - 0.050 * diff(usr[3:4]),
        tag, adj = c(1, 1), cex = 1.08, font = 2
    )
}

# 6. Fitted and residual panels -----------------------------------------

panel_titles <- rep(c(
    "Spatially implicit simulation",
    "Spatially explicit simulation"
), 2L)
panel_keys <- expand.grid(
    panel = panels, framework = frameworks,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

# Draw one observed LCBD and expectation panel.
fit_panel <- function(result, framework, panel, tag, title, panel_index) {
    observed <- result$observed[[framework]][[panel]]
    fits <- result$fits[[framework]][[panel]]
    S <- max(vapply(fits, function(z) max(z$selected$curve$richness), numeric(1)))
    x_range <- richness_coordinate(c(1, S))
    y_range <- c(-0.03, 1.03)

    plot(NA, NA, type = "n", axes = FALSE, xlab = "", ylab = "",
        xlim = x_range, ylim = y_range, xaxs = "i", yaxs = "i"
    )
    hexagons <- lapply(grains, function(grain) {
        data <- observed[[grain]]
        make_hexagons(
            richness_coordinate(data$richness), data$lcbd,
            x_range, y_range, 34, 22
        )
    })
    names(hexagons) <- grains

    # Draw every expectation below the point clouds first.
    for (grain in grains) {
        draw_curve(
            fits[[grain]]$selected$curve,
            grain_colours[[grain]], "longdash"
        )
        draw_curve(
            fits[[grain]]$parameterized$curve,
            grain_colours[[grain]], "dotted", 1.0
        )
        if (!is.null(fits[[grain]]$known)) {
            draw_curve(
                fits[[grain]]$known$curve,
                grain_colours[[grain]], "longdash", 1.25
            )
        }
    }

    # Each point cloud covers expectations from the other two grains.
    for (grain in grains) {
        draw_hexagons(hexagons[[grain]], grain_palettes[[grain]])
    }

    # Redraw both expectations only on their corresponding point cloud.
    for (grain in grains) {
        draw_curve(
            fits[[grain]]$selected$curve, grain_colours[[grain]], "longdash",
            1.25, hexagons[[grain]], hexagons[setdiff(grains, grain)]
        )
        draw_curve(
            fits[[grain]]$parameterized$curve,
            grain_colours[[grain]], "dotted", 1.0,
            hexagons[[grain]], hexagons[setdiff(grains, grain)]
        )
        if (!is.null(fits[[grain]]$known)) {
            draw_curve(
                fits[[grain]]$known$curve,
                grain_colours[[grain]], "longdash", 1.25,
                hexagons[[grain]], hexagons[setdiff(grains, grain)]
            )
        }
    }

    draw_x_axis(
        richness_ticks(S), result$model == "Model 2", inset_ends = TRUE
    )
    axis(2, at = seq(0, 1, 0.2),
        labels = if (panel_index == 1L) {
            format(seq(0, 1, 0.2), nsmall = 1)
        } else FALSE,
        las = 1, cex.axis = 0.88, lwd = 0, lwd.ticks = 0.75
    )
    box(col = "grey20", lwd = 0.85)
    title(main = title, cex.main = 0.92, line = 0.25)
    draw_tag(tag)

    if (panel_index == 1L) {
        legend("bottomleft", legend = unname(grain_labels[grains]),
            col = grain_colours[grains], lty = 1, lwd = 1.25,
            bty = "n", cex = 0.78, seg.len = 1.0, y.intersp = 0.88
        )
    }
    if (panel_index == 2L) {
        legend(
            "bottomleft", legend = c("Parameterized", "Quenched"),
            col = "grey25", lty = c("dotted", "longdash"),
            lwd = c(1.0, 1.25), bty = "n", cex = 0.72,
            seg.len = 1.5, y.intersp = 0.88
        )
    }
}

residual_ranges <- function(result) {
    output <- list()
    for (framework in frameworks) {
        y <- unlist(lapply(panels, function(panel) {
            item <- result$residuals$panels[[framework]][[panel]]
            c(item$data$residual, item$line$lower, item$line$upper, 0)
        }))
        output[[framework]] <- range(y, finite = TRUE)
        output[[framework]] <- output[[framework]] +
            c(-1, 1) * 0.06 * diff(output[[framework]])
    }
    output
}

# Draw one LCBD-deviation panel and its fitted relationship.
residual_panel <- function(
        result, framework, panel, tag, panel_index, y_ranges
) {
    item <- result$residuals$panels[[framework]][[panel]]
    data <- item$data
    line <- item$line
    x_range <- range(richness_coordinate(data$richness))
    x_range <- x_range + c(-1, 1) * 0.04 * diff(x_range)
    y_range <- y_ranges[[framework]]

    plot(NA, NA, type = "n", axes = FALSE, xlab = "", ylab = "",
        xlim = x_range, ylim = y_range, xaxs = "i", yaxs = "i"
    )
    polygon(
        c(richness_coordinate(line$richness),
            rev(richness_coordinate(line$richness))),
        c(line$lower, rev(line$upper)),
        col = adjustcolor(residual_colour, alpha.f = 0.10), border = NA
    )
    draw_hexagons(
        make_hexagons(
            richness_coordinate(data$richness), data$residual,
            x_range, y_range, 28, 20
        ),
        residual_palette
    )
    abline(h = 0, col = "grey55", lwd = 0.65, lty = 3)
    lines(
        richness_coordinate(line$richness), line$fitted,
        col = residual_colour, lwd = 1.05,
        lty = if (item$significant) 1 else 2
    )

    ticks <- pretty(range(data$richness), n = 4)
    ticks <- unique(round(ticks[ticks >= min(data$richness) &
        ticks <= max(data$richness)]))
    draw_x_axis(ticks, TRUE)
    y_ticks <- pretty(y_range, n = 4)
    axis(2, at = y_ticks,
        labels = if (panel_index == 1L) {
            formatC(y_ticks, format = "f", digits = 3)
        } else FALSE,
        las = 1, cex.axis = 0.82, lwd = 0, lwd.ticks = 0.75
    )
    box(col = "grey20", lwd = 0.85)
    draw_tag(tag)
}

# 7. Draw one complete two-row figure -----------------------------------

# Assemble the fitted and deviation rows for one simulation scenario.
draw_simulation_figure <- function(result) {
    fit_row_height <- if (result$model == "Model 2") 2.13 else 1.89
    layout(
        matrix(c(
            0, 1, 1, 0, 2, 2,
            3, 4, 5, 0, 6, 7,
            8, 9, 10, 0, 11, 12,
            0, 13, 13, 13, 13, 13
        ), 4, byrow = TRUE),
        widths = c(0.19, 1, 1, 0, 1, 1),
        heights = c(0.18, fit_row_height, 1.95, 0.24)
    )
    par(oma = rep(0.05, 4), family = font_family)
    y_ranges <- residual_ranges(result)

    draw_title_cell("")
    draw_title_cell("")
    draw_title_cell(expression(LCBD~"\u00d7"~M), 90, 1.06)
    for (i in seq_len(nrow(panel_keys))) {
        par(mar = c(if (result$model == "Model 2") 1.45 else 0.68,
                1.10, 1.55, 0.08),
            mgp = c(1.45, 0.30, 0), tcl = -0.20
        )
        fit_panel(
            result, panel_keys$framework[i], panel_keys$panel[i],
            letters[i], panel_titles[i], i
        )
    }

    draw_title_cell(expression(Delta*LCBD~"\u00d7"~M), 90, 1.06)
    for (i in seq_len(nrow(panel_keys))) {
        par(mar = c(2.00, 1.10, 0.12, 0.08),
            mgp = c(1.45, 0.30, 0), tcl = -0.20
        )
        residual_panel(
            result, panel_keys$framework[i], panel_keys$panel[i],
            letters[i + 4L], i, y_ranges
        )
    }
    draw_title_cell("Local richness", cex = 1.02)
    grid::grid.text(
        "Coverage deficit", x = grid::unit(0.272, "npc"), y = grid::unit(0.975, "npc"),
        gp = grid::gpar(fontfamily = font_family, fontsize = 14.5, fontface = 2)
    )
    grid::grid.text(
        "Sørensen", x = grid::unit(0.772, "npc"), y = grid::unit(0.975, "npc"),
        gp = grid::gpar(fontfamily = font_family, fontsize = 14.5, fontface = 2)
    )
}

# Export one simulation figure in the manuscript formats.
save_figure <- function(result, png_file) {
    ragg::agg_png(
        png_file, width = 12.2, height = 6.6, units = "in", res = 600,
        pointsize = 18, background = "white"
    )
    draw_simulation_figure(result)
    dev.off()
    ragg::agg_tiff(
        sub("\\.png$", ".tiff", png_file),
        width = 12.2, height = 6.6, units = "in", res = 600,
        pointsize = 18, background = "white", compression = "lzw"
    )
    draw_simulation_figure(result)
    dev.off()
    grDevices::pdf(
        sub("\\.png$", ".pdf", png_file),
        width = 12.2, height = 6.6, pointsize = 18,
        family = "Helvetica", useDingbats = FALSE, bg = "white"
    )
    draw_simulation_figure(result)
    dev.off()
}

# 8. Save Figure 2 --------------------------------

save_figure(model2_result, figure2_file)
