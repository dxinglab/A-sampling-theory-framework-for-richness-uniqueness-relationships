# Figure 4: BCI expectation, species representation and priority map
# The expectation is the fitted annealed Model 3.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggnewscale)
  library(patchwork)
  library(cowplot)
})
source("code/01-lcbd-models.R")

# Calculate and draw all Figure 4 panels for one LCBD framework.
make_figure4 <- function(figure4_framework) {

# 1. Settings -----------------------------------------------------------------

map_proportion <- 0.05
bar_proportions <- c(0.025, 0.05, 0.10)
rare_abundance_max <- 50L
framework_label <- if (figure4_framework == "coverage") "Coverage" else "Sorensen"

model3_output_file <- "process data/Empirical_annealed_expectations.rds"
bci_matrix_file <- "data/BCI_2015_species_matrices.rds"
bci_habitat_file <- "data/bci_habitat_official_bciex.rda"
output_stem <- if (figure4_framework == "coverage") {
  "figure/Figure4_BCI_Coverage_top5pct"
} else {
  "figure/FigureS1_BCI_Sorensen_top5pct"
}
figure_width <- 12.4
figure_height <- 11.6
bar_x <- setNames(c(0.75, 1.55, 2.35), sprintf("%g%%", 100 * bar_proportions))
bar_x_limits <- c(0, 3.10)
bar_width <- 0.44

col_raw <- "#5E3C99"
col_shared <- "#806CA0"
col_deviation <- "#E66101"
col_neither <- "#D9D9D9"
panel_border_colour <- "grey20"
panel_border_width <- 0.6
fit_panel_border_width <- 1.2

top_strip_text <- element_text(
  size = 15, face = "bold", margin = margin(t = 0, b = 5)
)

# 2. Read BCI scores and annealed fits ----------------------------------------

model3_output <- readRDS(model3_output_file)
bci_data <- readRDS(bci_matrix_file)
bci_reference <- bci_data$site_data[["20"]]
bci_reference$richness <- rowSums(bci_data$communities[["20"]] > 0)
bci_source <- model3_output$frameworks[[figure4_framework]]
bci_observed_all <- bci_source$observed$bci
bci_fits_all <- lapply(bci_source$fits$bci, function(fit) {
  list(
    parameterized = fit$fitted_x,
    selected = fit$model3,
    selected_name = "Model 3"
  )
})
bci_observed <- bci_observed_all[["20"]]
selected_fit <- bci_fits_all[["20"]]
expected <- selected_fit$selected$curve$expected[
  match(bci_observed$richness, selected_fit$selected$curve$richness)
]

site_data <- bci_reference %>%
  transmute(
    site_id = as.character(site_id),
    richness,
    raw_lcbd = bci_observed$lcbd,
    lcbd_deviation = bci_observed$lcbd - expected,
    coord_x,
    coord_y
  )

# 3. Read the matching BCI community matrix -----------------------------------

community <- as.matrix(bci_data$communities[["20"]])
sad <- bci_data$regional_abundance
rare_species <- names(sad)[sad < rare_abundance_max]

# 4. Select quadrats in the highest 5% of each ranking -------------------------

rank_exact <- function(value, site_id) {
  ord <- order(-value, site_id, na.last = TRUE)
  out <- rep(NA_integer_, length(value))
  out[ord] <- seq_along(ord)
  out
}

top_n <- ceiling(map_proportion * nrow(site_data))
site_data <- site_data %>%
  mutate(
    raw_selected = rank_exact(raw_lcbd, site_id) <= top_n,
    deviation_selected = rank_exact(lcbd_deviation, site_id) <= top_n
  )

# 5. Summarise selected species ------------------------------------------------

covered_species <- function(selected) {
  colnames(community)[colSums(community[selected, , drop = FALSE]) > 0]
}

partition_counts <- function(pool, raw_species, deviation_species) {
  raw <- intersect(raw_species, pool)
  deviation <- intersect(deviation_species, pool)
  c(
    `Raw LCBD only` = length(setdiff(raw, deviation)),
    `Shared by both` = length(intersect(raw, deviation)),
    `LCBD deviation only` = length(setdiff(deviation, raw)),
    `Selected by neither` = length(setdiff(pool, union(raw, deviation)))
  )
}

bar_counts <- bind_rows(lapply(bar_proportions, function(proportion) {
  n_selected <- ceiling(proportion * nrow(site_data))
  raw_selected <- rank_exact(site_data$raw_lcbd, site_data$site_id) <= n_selected
  deviation_selected <- rank_exact(
    site_data$lcbd_deviation, site_data$site_id
  ) <= n_selected
  raw_species <- covered_species(raw_selected)
  deviation_species <- covered_species(deviation_selected)

  bind_rows(
    data.frame(
      threshold = sprintf("%g%%", 100 * proportion),
      panel = "Species pool",
      count = partition_counts(
        colnames(community), raw_species, deviation_species
      )
    ),
    data.frame(
      threshold = sprintf("%g%%", 100 * proportion),
      panel = "Rare species",
      count = partition_counts(rare_species, raw_species, deviation_species)
    )
  )
})) %>%
  mutate(
    component = rep(
      c(
        "Raw LCBD only", "Shared by both",
        "LCBD deviation only", "Selected by neither"
      ), 2 * length(bar_proportions)
    ),
    panel = factor(panel, levels = c("Species pool", "Rare species")),
    threshold = factor(
      threshold, levels = sprintf("%g%%", 100 * bar_proportions)
    )
  ) %>%
  group_by(panel, threshold) %>%
  mutate(
    percentage = 100 * count / sum(count),
    ymax = cumsum(percentage),
    ymin = ymax - percentage,
    ymid = (ymin + ymax) / 2
  ) %>%
  ungroup() %>%
  mutate(x_plot = unname(bar_x[as.character(threshold)]))

bar_counts$component <- factor(
  bar_counts$component,
  levels = c(
    "Raw LCBD only", "Shared by both",
    "LCBD deviation only", "Selected by neither"
  )
)

# 6. Read the official BCI habitat map -----------------------------------------

habitat_env <- new.env(parent = emptyenv())
load(bci_habitat_file, envir = habitat_env)
bci_habitat <- habitat_env$bci_habitat

habitat_colors <- c(
  low_plateau = "#377EB8", hi_plateau = "#E41A1C",
  slope = "#4DAF4A", swamp = "#984EA3", stream = "#00CED1",
  young = "#FF7F00", mixed = "#A65628"
)
habitat_labels <- c(
  low_plateau = "Low plateau",
  hi_plateau = "High plateau",
  slope = "Slope",
  swamp = "Swamp",
  stream = "Streamside",
  young = "Young forest",
  mixed = "Mixed habitats"
)

richness_limits <- as.numeric(
  quantile(site_data$richness, c(0.05, 0.95), na.rm = TRUE)
)
richness_breaks <- breaks_pretty(n = 4)(richness_limits)
richness_breaks <- richness_breaks[
  richness_breaks >= min(site_data$richness) &
    richness_breaks <= max(site_data$richness)
]

richness_to_radius <- function(richness) {
  richness <- oob_squish(richness, richness_limits)
  radius <- rescale(
    sqrt(richness), to = c(6, 12), from = sqrt(richness_limits)
  )
  radius[!is.finite(radius)] <- 5
  radius
}

make_wedge <- function(x0, y0, radius, angle0, angle1, n = 45,
                       include_center = TRUE) {
  angle <- seq(angle0, angle1, length.out = n)
  if (include_center) {
    return(data.frame(
      x = c(x0, x0 + radius * cos(angle)),
      y = c(y0, y0 + radius * sin(angle))
    ))
  }
  data.frame(x = x0 + radius * cos(angle), y = y0 + radius * sin(angle))
}

# 7. BCI fitted panel a ---------------------------------------------------

grain_colours <- c("5" = "#3378A6", "20" = "#5D8B73", "100" = "#B27645")
grain_labels <- c(
  "5" = "5 m × 5 m",
  "20" = "20 m × 20 m",
  "100" = "100 m × 100 m"
)
richness_coordinate <- function(x) log10(x + 1)

smooth_curve <- function(curve, grain, expectation, n = 1200L) {
  curve <- curve[curve$richness >= 1 & is.finite(curve$expected), ]
  x <- richness_coordinate(curve$richness)
  x_draw <- seq(min(x), max(x), length.out = n)
  data.frame(
    x = x_draw,
    expected = splinefun(x, curve$expected, method = "monoH.FC")(x_draw),
    grain = grain, expectation = expectation
  )
}

# Aggregate the BCI point cloud into hexagonal plotting cells.
make_hexagons <- function(data, x_range, y_range, bins = 31L) {
  built <- ggplot_build(
    ggplot(data, aes(x, lcbd)) +
      geom_hex(bins = bins) +
      scale_x_continuous(limits = x_range, expand = expansion(mult = 0)) +
      scale_y_continuous(limits = y_range, expand = expansion(mult = 0))
  )$data[[1L]]
  built <- built[
    is.finite(built$x) & is.finite(built$y) & is.finite(built$count),
  ]
  dx <- built$width[1L] / 2
  dy <- built$height[1L] / sqrt(3) / 2
  corners <- hexbin::hexcoords(dx, dy, n = 1L)
  maximum <- max(built$count)
  vertices <- bind_rows(lapply(seq_len(nrow(built)), function(i) {
    data.frame(
      cell = i,
      x = built$x[i] + corners$x,
      y = built$y[i] + corners$y,
      density = (log1p(built$count[i]) / log1p(maximum))^2.15
    )
  }))
  polygons <- lapply(split(vertices[, c("x", "y")], vertices$cell), function(z) {
    xy <- rbind(as.matrix(z), as.matrix(z[1L, ]))
    sf::st_polygon(list(xy))
  })
  list(vertices = vertices, mask = sf::st_union(sf::st_sfc(polygons)))
}

mask_curve <- function(curve, mask) {
  visible <- lengths(sf::st_intersects(
    sf::st_as_sf(curve, coords = c("x", "expected")), mask
  )) > 0L
  curve$expected[!visible] <- NA_real_
  curve
}

# Draw the BCI fitted relationship and grain-specific expectations.
make_fit_panel <- function() {
  grains <- c("5", "20", "100")
  curves <- bind_rows(lapply(c("5", "20", "100"), function(grain) {
    fit <- bci_fits_all[[grain]]
    bind_rows(
      smooth_curve(
        fit$selected$curve, grain, "Conditional expectation"
      ),
      smooth_curve(
        fit$parameterized$curve, grain, "Parameterized expectation"
      )
    )
  }))
  S <- max(vapply(bci_fits_all, function(x) {
    max(x$selected$curve$richness)
  }, numeric(1)))
  ticks <- sort(unique(c(1, 10, 50, 100, S)))
  ticks <- ticks[ticks <= S]
  x_range <- richness_coordinate(c(1, S))
  y_range <- c(0, 1.02)
  hexagons <- lapply(grains, function(grain) {
    data <- bci_observed_all[[grain]] %>%
      mutate(x = richness_coordinate(richness))
    make_hexagons(data, x_range, y_range)
  })
  names(hexagons) <- grains

  plot <- ggplot(data.frame(panel_title = "BCI tropical forest")) +
    geom_line(
      data = curves,
      aes(x, expected, colour = grain, linetype = expectation,
          linewidth = expectation,
          group = interaction(grain, expectation))
    )

  for (grain in grains) {
    alpha_palette <- colorRampPalette(c(
      adjustcolor(grain_colours[[grain]], alpha.f = 0.18),
      adjustcolor(grain_colours[[grain]], alpha.f = 0.45),
      adjustcolor(grain_colours[[grain]], alpha.f = 0.75),
      grain_colours[[grain]]
    ), alpha = TRUE)(256L)
    rgba <- col2rgb(alpha_palette, alpha = TRUE) / 255
    density_palette <- rgb(
      rgba[1, ] * rgba[4, ] + 1 - rgba[4, ],
      rgba[2, ] * rgba[4, ] + 1 - rgba[4, ],
      rgba[3, ] * rgba[4, ] + 1 - rgba[4, ]
    )
    grain_curves <- mask_curve(
      curves[curves$grain == grain, ], hexagons[[grain]]$mask
    )
    plot <- plot +
      geom_polygon(
        data = hexagons[[grain]]$vertices,
        aes(x, y, group = cell, fill = density),
        colour = NA, show.legend = FALSE
      ) +
      scale_fill_gradientn(
        colours = density_palette, limits = c(0, 1), guide = "none"
      ) +
      ggnewscale::new_scale_fill() +
      geom_line(
        data = grain_curves,
        aes(x, expected, colour = grain, linetype = expectation,
            linewidth = expectation,
            group = interaction(grain, expectation)),
        na.rm = TRUE
      )
  }

  plot +
    facet_wrap(~panel_title) +
    scale_colour_manual(
      values = grain_colours, breaks = c("5", "20", "100"),
      labels = grain_labels[c("5", "20", "100")], name = NULL
    ) +
    scale_linetype_manual(
      values = c(
        `Parameterized expectation` = "dotted",
        `Conditional expectation` = "solid"
      ),
      breaks = c(
        "Parameterized expectation", "Conditional expectation"
      ),
      name = NULL
    ) +
    scale_linewidth_manual(
      values = c(
        `Parameterized expectation` = 0.52,
        `Conditional expectation` = 0.68
      ),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = richness_coordinate(ticks), labels = ticks,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "Local richness", y = expression(LCBD~"\u00d7"~M),
      tag = "a"
    ) +
    guides(
      colour = guide_legend(
        order = 1, override.aes = list(linetype = "solid", linewidth = 0.8)
      ),
      linetype = "none"
    ) +
    theme_classic(base_size = 14, base_family = "Helvetica") +
    theme(
      panel.border = element_rect(
        fill = NA, colour = panel_border_colour,
        linewidth = fit_panel_border_width
      ),
      axis.line = element_blank(),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(size = 16, margin = margin(r = 21)),
      axis.text = element_text(size = 14, colour = "grey30"),
      strip.text.x = element_blank(),
      strip.background = element_blank(),
      plot.tag = element_text(size = 14.2, face = "bold"),
      plot.tag.position = c(0.97, 1.045),
      legend.position = c(0.035, 0.035),
      legend.justification = c(0, 0),
      legend.box = "vertical",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.width = grid::unit(0.62, "cm"),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.spacing.y = grid::unit(0.02, "cm"),
      plot.margin = margin(3, 12, 0, 12)
    )
}

# 8. Habitat and priority map panel c -----------------------------------

# Draw the BCI habitat map and the two priority selections.
make_map_panel <- function(data) {
  points <- data %>%
    mutate(
      priority_class = case_when(
        raw_selected & deviation_selected ~ "Both",
        raw_selected ~ "Raw LCBD",
        deviation_selected ~ "LCBD deviation",
        TRUE ~ "Not selected"
      ),
      radius = richness_to_radius(richness)
    ) %>%
    arrange(desc(priority_class == "Not selected"), desc(richness)) %>%
    mutate(polygon_id = seq_len(n()))

  polygons <- bind_rows(lapply(seq_len(nrow(points)), function(i) {
    row <- points[i, ]
    if (row$priority_class == "Both") {
      split_angle <- 135 * pi / 180
      raw_half <- make_wedge(
        row$coord_x, row$coord_y, row$radius,
        split_angle, split_angle + pi
      )
      deviation_half <- make_wedge(
        row$coord_x, row$coord_y, row$radius,
        split_angle - pi, split_angle
      )
      raw_half$group <- paste0(row$polygon_id, "_raw")
      raw_half$fill_group <- "Raw LCBD"
      raw_half$selected <- TRUE
      deviation_half$group <- paste0(row$polygon_id, "_deviation")
      deviation_half$fill_group <- "LCBD deviation"
      deviation_half$selected <- TRUE
      return(bind_rows(raw_half, deviation_half))
    }

    circle <- make_wedge(
      row$coord_x, row$coord_y, row$radius, 0, 2 * pi,
      include_center = FALSE
    )
    circle$group <- paste0(row$polygon_id, "_circle")
    circle$fill_group <- row$priority_class
    circle$selected <- row$priority_class != "Not selected"
    circle
  }))

  background <- filter(polygons, !selected)
  selected <- filter(polygons, selected)
  legend_x <- min(data$coord_x)
  legend_y <- min(data$coord_y)
  legend_sizes <- 0.52 * richness_to_radius(richness_breaks)

  ggplot() +
    geom_tile(
      data = bci_habitat,
      aes(x = x, y = y, fill = habitat),
      alpha = 0.45
    ) +
    scale_fill_manual(
      values = habitat_colors, labels = habitat_labels, name = "Habitat",
      guide = guide_legend(
        order = 1, ncol = 3, byrow = TRUE,
        title.position = "top", title.hjust = 0.5
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_polygon(
      data = background,
      aes(x, y, group = group, fill = fill_group),
      color = "white", linewidth = 0.08, alpha = 0.55
    ) +
    geom_polygon(
      data = selected,
      aes(x, y, group = group, fill = fill_group),
      color = "white", linewidth = 0.35, alpha = 0.95
    ) +
    scale_fill_manual(
      values = c(
        `Raw LCBD` = col_raw,
        `LCBD deviation` = col_deviation,
        `Not selected` = "grey75"
      ),
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = data.frame(
        x = rep(legend_x, 3), y = rep(legend_y, 3),
        selection = factor(
          c("Raw LCBD", "LCBD deviation", "Not selected"),
          levels = c("Raw LCBD", "LCBD deviation", "Not selected")
        )
      ),
      aes(x, y, fill = selection),
      shape = 21, size = 3, color = "white", alpha = 0
    ) +
    scale_fill_manual(
      values = c(
        `Raw LCBD` = col_raw,
        `LCBD deviation` = col_deviation,
        `Not selected` = "grey75"
      ),
      labels = c(
        `Raw LCBD` = "Top 5% by raw LCBD",
        `LCBD deviation` = "Top 5% by LCBD deviation",
        `Not selected` = "Not in Top 5%"
      ),
      name = "Quadrat ranking",
      guide = guide_legend(
        order = 2, ncol = 1, byrow = TRUE,
        title.position = "top", title.hjust = 0.5,
        override.aes = list(size = 4.8, alpha = c(1, 1, 0.5)),
        theme = theme(legend.key.width = grid::unit(0.28, "cm"))
      )
    ) +
    geom_point(
      data = data.frame(
        x = legend_x, y = legend_y,
        richness = richness_breaks, point_size = legend_sizes
      ),
      aes(x, y, size = point_size),
      shape = 21, fill = "grey60", color = "white", alpha = 0
    ) +
    scale_size_identity(
      name = "Richness", breaks = legend_sizes, labels = richness_breaks,
      guide = guide_legend(
        order = 3, ncol = 3, byrow = TRUE,
        title.position = "top", title.hjust = 0.5,
        override.aes = list(
          shape = 21, fill = "grey55", color = "white", alpha = 0.85
        )
      )
    ) +
    scale_x_continuous(breaks = seq(0, 1000, 200)) +
    scale_y_continuous(breaks = seq(0, 500, 100)) +
    coord_fixed(
      ratio = 1, xlim = c(-15, 1015), ylim = c(-15, 515), expand = FALSE
    ) +
    labs(x = "X coordinate (m)", y = "Y coordinate (m)", tag = "c") +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 16),
      axis.title.x = element_text(size = 16, margin = margin(t = 2)),
      axis.title.y = element_text(size = 16, margin = margin(r = 18)),
      axis.text = element_text(size = 14, color = "grey30"),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.justification = "center",
      legend.direction = "vertical",
      legend.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 14),
      legend.key.width = grid::unit(0.60, "cm"),
      legend.key.height = grid::unit(0.48, "cm"),
      legend.spacing.x = grid::unit(0.12, "cm"),
      legend.spacing.y = grid::unit(0.15, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.spacing = grid::unit(1, "pt"),
      plot.tag = element_text(size = 14.2, face = "bold"),
      plot.tag.position = c(0.965, 0.98),
      plot.margin = margin(0, 5, 2, 7)
    )
}

# 9. Species representation panel b -------------------------------------

bar_palette <- c(
  `Raw LCBD only` = col_raw,
  `Shared by both` = col_shared,
  `LCBD deviation only` = col_deviation,
  `Selected by neither` = col_neither
)

# Draw species representation by selection overlap.
make_bar_plot <- function(data) {
  inside <- filter(data, percentage >= 4)
  outside <- filter(data, percentage < 4) %>%
    mutate(label_side = if_else(panel == "Species pool", -1, 1))
  panel_labels <- c(
    `Species pool` = paste0("Species pool (", ncol(community), ")"),
    `Rare species` = paste0("Rare species (", length(rare_species), ")")
  )

  ggplot(data, aes(x = x_plot, y = percentage, fill = component)) +
    geom_col(width = bar_width, position = position_stack(reverse = TRUE)) +
    geom_text(
      data = inside,
      aes(y = ymid, label = sprintf("%.0f%%", percentage)),
      color = ifelse(inside$component == "Selected by neither", "grey25", "white"),
      fontface = "bold", size = 5
    ) +
    geom_segment(
      data = outside,
      aes(
        x = x_plot + 0.20 * label_side,
        xend = x_plot + 0.24 * label_side,
        y = ymid, yend = ymid + 4.5
      ),
      inherit.aes = FALSE, color = col_raw, linewidth = 0.7
    ) +
    geom_text(
      data = outside,
      aes(
        x = x_plot + 0.28 * label_side, y = ymid + 4.5,
        label = sprintf("%.0f%%", percentage)
      ),
      inherit.aes = FALSE,
      hjust = ifelse(outside$label_side < 0, 1, 0),
      color = col_raw, fontface = "bold", size = 5
    ) +
    scale_fill_manual(
      values = bar_palette,
      labels = c(
        `Raw LCBD only` = "Raw LCBD only",
        `Shared by both` = "Shared",
        `LCBD deviation only` = "LCBD deviation only",
        `Selected by neither` = "Neither"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_continuous(
      breaks = unname(bar_x), labels = names(bar_x),
      limits = bar_x_limits, expand = expansion(mult = 0)
    ) +
    facet_grid(~panel, labeller = as_labeller(panel_labels)) +
    labs(
      x = "Selection threshold",
      y = "Percentage of species",
      fill = NULL,
      tag = "b"
    ) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE, clip = "off") +
    theme_classic(base_size = 14, base_family = "Helvetica") +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(
        color = panel_border_colour, fill = NA,
        linewidth = panel_border_width
      ),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(size = 16, margin = margin(r = 2)),
      axis.text = element_text(size = 14, color = "grey30"),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.position = "bottom",
      legend.box.just = "center",
      legend.justification = "center",
      legend.text = element_text(size = 11.5),
      legend.key.width = grid::unit(0.46, "cm"),
      legend.spacing.x = grid::unit(0.04, "cm"),
      legend.margin = margin(t = 2, r = 0, b = 0, l = -8),
      legend.box.spacing = grid::unit(0, "pt"),
      panel.spacing.x = grid::unit(0, "pt"),
      strip.text.x = top_strip_text,
      strip.background = element_blank(),
      plot.tag = element_text(size = 14.2, face = "bold"),
      plot.tag.position = c(0.98, 0.98),
      plot.margin = margin(3, 25, 0, 7)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

panel_a <- make_fit_panel()
panel_b <- make_bar_plot(bar_counts)
panel_c <- make_map_panel(site_data)

top_row <- cowplot::plot_grid(
  panel_a, panel_b, nrow = 1, rel_widths = c(0.90, 1.55),
  align = "h", axis = "t"
)
figure4 <- cowplot::plot_grid(
  top_row, panel_c, ncol = 1, rel_heights = c(0.94, 1.62)
)

# 10. Export Figure 4 ----------------------------------------------------

ggsave(
  paste0(output_stem, ".png"), figure4,
  width = figure_width, height = figure_height, dpi = 600, bg = "white"
)

ragg::agg_tiff(
  paste0(output_stem, ".tiff"),
  width = figure_width, height = figure_height, units = "in", res = 600,
  background = "white", compression = "lzw"
)
print(figure4)
grDevices::dev.off()

if (capabilities("aqua")) {
  grDevices::quartz(
    type = "pdf", file = paste0(output_stem, ".pdf"),
    width = figure_width, height = figure_height, bg = "white"
  )
} else {
  grDevices::pdf(
    paste0(output_stem, ".pdf"),
    width = figure_width, height = figure_height,
    family = "Helvetica", useDingbats = FALSE, bg = "white"
  )
}
print(figure4)
grDevices::dev.off()

invisible(output_stem)
}

invisible(lapply(c("coverage", "sorensen"), make_figure4))
