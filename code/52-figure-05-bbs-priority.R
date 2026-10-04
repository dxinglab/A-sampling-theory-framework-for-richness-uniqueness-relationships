# Figure 5: BBS highest-5% route priorities and species representation.
# Run after 51-prepare-bbs-data.R.

suppressPackageStartupMessages({
    library(sf)
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(scales)
})
sf::sf_use_s2(FALSE)
source("code/01-lcbd-models.R")

# 1. Paths and settings ---------------------------------------------------

model_file <- "process data/Figure5_BBS_annealed_expectations.rds"
species_file <- "process data/Figure5_BBS_species_priority_two_frameworks.rds"
bcc_file <- "process data/Figure5_BBS_USFWS_BCC_2021_comparison.rds"
wdpa_polygon_files <- file.path(
    "data",
    paste0("WDPA_Mar2026_Public_shp_", 0:2),
    "WDPA_Mar2026_Public_shp-polygons.shp"
)
bcr_file <- "data/BCR_Terrestrial_master.shp"
output_dir <- "figure"

map_crs <- paste(
    "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96",
    "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
)
map_bbox <- c(xmin = -170, ymin = 15, xmax = -50, ymax = 75)
grid_cell_m <- 100000

focus_bcr <- 10L
figure_width <- 12
figure_height <- 10.4

raw_colour <- "#5E3C99"
deviation_colour <- "#E66101"
not_selected_colour <- "grey75"
coverage_colours <- c(
    "#F7FCF0", "#D9F0D3", "#A6DBA0", "#5AB4AC", "#2B8CBE", "#084081"
)

# 2. Read the BCR-specific Model 3 deviations and map data ---------------

model_result <- readRDS(model_file)
species_result <- readRDS(species_file)
bcc_result <- readRDS(bcc_file)

wide <- model_result$site_data %>%
    mutate(
        top_raw_coverage = coverage_raw_selected,
        top_delta_coverage = coverage_selected,
        top_raw_sorensen = sorensen_raw_selected,
        top_delta_sorensen = sorensen_selected,
        surveyID = site_id
    )

countries <- maps::map(
    "world", regions = c("USA$", "USA:Alaska", "Canada", "Mexico"),
    fill = TRUE, plot = FALSE
) %>%
    st_as_sf() %>%
    st_transform(4326) %>%
    st_make_valid() %>%
    st_crop(st_as_sfc(st_bbox(map_bbox, crs = st_crs(4326)))) %>%
    st_transform(map_crs)

wdpa_query <- function(path) {
    layer <- tools::file_path_sans_ext(basename(path))
    query <- paste0(
        'SELECT * FROM "', layer, '" WHERE ',
        "ISO3 IN ('USA','CAN','MEX') AND ",
        "STATUS IN ('Designated','Established') AND REALM <> 'Marine'"
    )
    part <- st_read(dirname(path), query = query, quiet = TRUE) %>%
        st_transform(map_crs) %>%
        st_make_valid()
    st_sf(geometry = st_geometry(part))
}

wdpa <- do.call(rbind, lapply(wdpa_polygon_files, wdpa_query))
land <- st_union(st_geometry(countries))
grid <- st_make_grid(land, cellsize = grid_cell_m, square = FALSE)
coverage <- st_sf(cell_id = seq_along(grid), geometry = grid) %>%
    st_intersection(st_sf(land_id = 1L, geometry = land)) %>%
    select(cell_id)
coverage$land_area_km2 <- as.numeric(st_area(coverage)) / 1e6

wdpa_hits <- st_intersects(coverage, wdpa)
progress <- txtProgressBar(min = 0, max = nrow(coverage), style = 3)
coverage$protected_area_km2 <- vapply(seq_len(nrow(coverage)), function(i) {
    setTxtProgressBar(progress, i)
    ids <- wdpa_hits[[i]]
    if (!length(ids)) return(0)
    clipped <- suppressWarnings(st_intersection(
        st_geometry(wdpa)[ids], st_geometry(coverage)[i]
    ))
    clipped <- clipped[!st_is_empty(clipped)]
    if (!length(clipped)) return(0)
    as.numeric(st_area(st_union(clipped))) / 1e6
}, numeric(1))
close(progress)
coverage$protected_percent <- pmin(
    100, 100 * coverage$protected_area_km2 / coverage$land_area_km2
)

great_lakes <- maps::map(
    "lakes", regions = c("Great Lakes", "Lake St. Clair"),
    fill = TRUE, plot = FALSE
) %>%
    st_as_sf() %>%
    st_make_valid() %>%
    st_transform(4326) %>%
    st_transform(map_crs)
bcr_boundaries <- st_read(bcr_file, quiet = TRUE) %>%
    mutate(
        BCR = as.integer(BCR),
        name_en = tools::toTitleCase(tolower(BCRNAME))
    ) %>%
    filter(BCR %in% model_result$settings$retained_bcr) %>%
    st_make_valid() %>%
    group_by(BCR, name_en) %>%
    summarise(do_union = TRUE, .groups = "drop") %>%
    st_transform(st_crs(coverage))
routes <- st_as_sf(
    wide, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE
) %>%
    st_transform(st_crs(coverage))
route_bounds <- st_bbox(routes)
map_xlim <- c(route_bounds[["xmin"]] - 250000,
    route_bounds[["xmax"]] + 250000)
map_ylim <- c(route_bounds[["ymin"]] - 150000,
    route_bounds[["ymax"]] + 150000)

# 3. Prepare the protected-area background -------------------------------

coverage_palette <- scales::colour_ramp(coverage_colours)
centres <- st_point_on_surface(st_geometry(coverage))
centre_sf <- st_sf(
    cell_id = seq_along(centres), geometry = centres, crs = st_crs(coverage)
)
neighbours <- st_is_within_distance(centre_sf, centre_sf, dist = 250000)
coverage_values <- coverage$protected_percent
coverage$protected_smooth <- vapply(seq_along(neighbours), function(i) {
    ids <- neighbours[[i]]
    weights <- rep(1, length(ids))
    weights[ids == i] <- 3
    weighted.mean(coverage_values[ids], weights, na.rm = TRUE)
}, numeric(1))
positive <- coverage$protected_smooth[coverage$protected_smooth > 0]
coverage_cap <- max(as.numeric(quantile(positive, 0.98, na.rm = TRUE)), 1)
coverage$protected_plot <- pmin(coverage$protected_smooth, coverage_cap)

richness_limits <- as.numeric(quantile(routes$richness, c(0.05, 0.95)))
richness_breaks <- c(30, 50, 70, 90)
reference_radius_range_m <- c(42000, 72000)
radius_at_50 <- rescale(
    sqrt(50), to = reference_radius_range_m, from = sqrt(richness_limits)
)
reference_legend_sizes <- rescale(
    rescale(
        sqrt(richness_breaks), to = reference_radius_range_m,
        from = sqrt(richness_limits)
    ),
    to = c(2.6, 4.8)
)
legend_max_size <- reference_legend_sizes[richness_breaks == 50]
radius_range_m <- c(30000, radius_at_50)

richness_to_radius <- function(values) {
    values <- oob_squish(values, richness_limits)
    rescale(sqrt(values), to = radius_range_m, from = sqrt(richness_limits))
}

# 4. Build split symbols for one framework -------------------------------

make_sector <- function(x, y, radius, start, end, group, fill_class) {
    angle <- seq(start, end, length.out = 49) * pi / 180
    data.frame(
        x = c(x, x + radius * cos(angle)),
        y = c(y, y + radius * sin(angle)),
        group = group, fill_class = fill_class
    )
}

make_circle <- function(x, y, radius, group, fill_class) {
    angle <- seq(0, 2 * pi, length.out = 65)
    data.frame(
        x = x + radius * cos(angle), y = y + radius * sin(angle),
        group = group, fill_class = fill_class
    )
}

# Convert route selections to circles and split symbols.
route_polygons <- function(route_sf, framework) {
    raw <- st_drop_geometry(route_sf)[[paste0("top_raw_", framework)]]
    deviation <- st_drop_geometry(route_sf)[[paste0("top_delta_", framework)]]
    xy <- st_coordinates(route_sf)
    points <- st_drop_geometry(route_sf) %>%
        mutate(
            x = xy[, 1L], y = xy[, 2L],
            raw_top = raw, deviation_top = deviation,
            selection = case_when(
                raw_top & deviation_top ~ "Both",
                raw_top ~ "Raw", deviation_top ~ "Deviation",
                TRUE ~ "Not selected"
            ),
            radius = richness_to_radius(richness),
            draw_layer = case_when(
                selection == "Not selected" ~ 1L,
                selection == "Both" ~ 3L,
                TRUE ~ 2L
            )
        ) %>%
        arrange(draw_layer, desc(richness)) %>%
        mutate(point_id = row_number())

    bind_rows(lapply(seq_len(nrow(points)), function(i) {
        row <- points[i, ]
        if (row$selection == "Both") {
            return(bind_rows(
                make_sector(
                    row$x, row$y, row$radius, 135, 315,
                    paste0(row$point_id, "_raw"), "Raw"
                ),
                make_sector(
                    row$x, row$y, row$radius, -45, 135,
                    paste0(row$point_id, "_deviation"), "Deviation"
                )
            ))
        }
        fill <- switch(
            row$selection, Raw = "Raw", Deviation = "Deviation",
            `Not selected` = "Not selected"
        )
        make_circle(
            row$x, row$y, row$radius,
            paste0(row$point_id, "_circle"), fill
        )
    })) %>%
        mutate(alpha = ifelse(fill_class == "Not selected", 0.22, 0.97))
}

# 5. Map and matching legend ---------------------------------------------

# Draw the BBS route-priority map for one LCBD framework.
make_map <- function(framework, title, highlight_boundary = NULL) {
    polygons <- route_polygons(routes, framework) %>%
        mutate(fill = case_when(
            fill_class == "Raw" ~ raw_colour,
            fill_class == "Deviation" ~ deviation_colour,
            TRUE ~ not_selected_colour
        ))
    background <- coverage %>%
        mutate(fill = coverage_palette(rescale(
            protected_plot, to = c(0, 1), from = c(0, coverage_cap)
        )))
    highlight_layer <- if (is.null(highlight_boundary)) NULL else geom_sf(
        data = highlight_boundary, fill = NA, colour = "grey18",
        linewidth = 0.50
    )

    ggplot() +
        geom_sf(data = countries, fill = "#F9F8F2", colour = NA) +
        geom_sf(data = background, aes(fill = fill), colour = NA, alpha = 0.90) +
        geom_sf(data = great_lakes, fill = "white", colour = NA) +
        geom_polygon(
            data = polygons,
            aes(x, y, group = group, fill = fill, alpha = I(alpha)),
            colour = "white", linewidth = 0.10, inherit.aes = FALSE
        ) +
        geom_sf(
            data = bcr_boundaries, fill = NA, colour = "grey48",
            linewidth = 0.18
        ) +
        highlight_layer +
        scale_fill_identity() +
        coord_sf(
            crs = st_crs(coverage), xlim = map_xlim, ylim = map_ylim,
            expand = FALSE, clip = "on"
        ) +
        labs(title = title) +
        theme_void(base_size = 14, base_family = "sans") +
        theme(
            plot.title = element_text(size = 16, hjust = 1),
            plot.margin = margin(2, 2, 2, 2),
            plot.background = element_rect(fill = "white", colour = NA)
        )
}

legend_theme <- theme_void(base_family = "sans") +
    theme(
        plot.margin = margin(0),
        plot.background = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = NA, colour = NA)
    )

# Assemble the map legends at the specified size.
make_legend <- function(text_size = 12, compact = FALSE) {
    gradient <- data.frame(x = seq(0.25, 49.75, length.out = 100), y = 1) %>%
        mutate(fill = coverage_palette(rescale(x, to = c(0, 1))))
    coverage_key <- ggplot(gradient, aes(x, y, fill = fill)) +
        geom_tile(height = 0.82) + scale_fill_identity() +
        scale_x_continuous(breaks = seq(0, 50, 10), limits = c(0, 50),
            expand = expansion(mult = c(0, 0))) +
        coord_cartesian(ylim = c(0.45, 1.55), clip = "off") +
        labs(title = "Protected land (%)") +
        theme_classic(base_size = text_size) +
        theme(
            axis.title = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.line = element_blank(),
            plot.title = element_text(size = text_size + 1, hjust = 0),
            plot.margin = margin(0),
            plot.background = element_rect(fill = NA, colour = NA),
            panel.background = element_rect(fill = NA, colour = NA)
        )

    selection <- data.frame(
        label = c(
            "Top 5% by raw LCBD",
            "Top 5% by LCBD deviation",
            "Not in Top 5%"
        ),
        colour = c(raw_colour, deviation_colour, not_selected_colour),
        y = c(0.61, 0.34, 0.07)
    )
    selection_key <- ggplot(selection) +
        geom_point(
            aes(0.06, y), colour = selection$colour,
            size = legend_max_size
        ) +
        geom_text(aes(0.14, y, label = label), hjust = 0,
            size = text_size / ggplot2::.pt) +
        annotate("text", 0, 0.98, label = "Route ranking within each BCR", hjust = 0,
            vjust = 1, size = (text_size + 1) / ggplot2::.pt) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
        legend_theme

    point_sizes <- legend_max_size *
        richness_to_radius(richness_breaks) / max(radius_range_m)
    richness_key <- ggplot() +
        annotate("text", 0, 0.98, label = "Richness", hjust = 0,
            vjust = 1, size = (text_size + 1) / ggplot2::.pt) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
        legend_theme
    for (i in seq_along(richness_breaks)) {
        point_x <- 0.05 + 0.24 * (i - 1L)
        point_y <- 0.30
        richness_key <- richness_key +
            annotate("point", point_x, point_y, size = point_sizes[i],
                colour = "grey50") +
            annotate("text", point_x + 0.06, point_y,
                label = richness_breaks[i], hjust = 0,
                size = text_size / ggplot2::.pt)
    }

    if (compact) {
        ggdraw() +
            draw_plot(coverage_key, 0.02, 0.50, 0.42, 0.45) +
            draw_plot(richness_key, 0.02, 0.04, 0.42, 0.43) +
            draw_plot(selection_key, 0.50, 0.05, 0.48, 0.84)
    } else {
        ggdraw() +
            draw_plot(coverage_key, 0.04, 0.76, 0.88, 0.22) +
            draw_plot(selection_key, 0.04, 0.38, 0.96, 0.34) +
            draw_plot(richness_key, 0.04, 0.02, 0.96, 0.30)
    }
}

smooth_curve <- function(curve, n = 1200L) {
    curve <- curve[curve$richness >= 1 & is.finite(curve$expected), ]
    x <- log10(curve$richness + 1)
    x_draw <- seq(min(x), max(x), length.out = n)
    data.frame(
        x = x_draw,
        expected = splinefun(x, curve$expected, method = "monoH.FC")(x_draw)
    )
}

# Draw the fitted BBS relationship for one framework.
make_bbs_fit_panel <- function(framework) {
    framework_name <- if (framework == "coverage") "Coverage" else "Sorensen"
    bbs_fit <- model_result$bcr10_fit_panel[[framework_name]]
    bbs_observed <- bbs_fit$observed
    observed <- transform(bbs_observed, x = log10(richness + 1))
    parameterized <- smooth_curve(bbs_fit$parameterized$curve)
    parameter <- subset(
        model_result$parameters,
        BCR == focus_bcr & framework == framework_name
    )
    richness <- seq_len(parameter$S)
    conditional_expected <- annealed_expectation_complete(
        framework,
        parameter$S,
        parameter$M,
        parameter$x,
        parameter$c,
        richness
    )
    conditional <- smooth_curve(data.frame(
        richness = richness,
        expected = conditional_expected
    ))
    ticks <- c(1, 10, 50, 100, 259)
    ggplot(observed, aes(x, lcbd)) +
        geom_hex(
            aes(fill = after_stat(
                (log1p(count) / log1p(max(count)))^2.15
            )), bins = 27, colour = NA
        ) +
        scale_fill_gradientn(
            colours = c("#EEEAF2", "#C8BED5", "#9A89B0", "#6D5A8E"),
            limits = c(0, 1), guide = "none"
        ) +
        geom_line(
            data = parameterized,
            aes(x, expected, linetype = "Parameterized expectation"),
            colour = "#6D5A8E", linewidth = 0.60,
            inherit.aes = FALSE
        ) +
        geom_line(
            data = conditional,
            aes(x, expected, linetype = "Conditional expectation"),
            colour = "#6D5A8E", linewidth = 0.78,
            inherit.aes = FALSE
        ) +
        scale_linetype_manual(
            values = c(
                "Parameterized expectation" = "dotted",
                "Conditional expectation" = "solid"
            ), breaks = c(
                "Parameterized expectation", "Conditional expectation"
            ), name = NULL
        ) +
        scale_x_continuous(
            breaks = log10(ticks + 1), labels = ticks,
            expand = expansion(mult = c(0, 0.015))
        ) +
        scale_y_continuous(
            limits = c(0, 1), breaks = seq(0, 1, 0.2),
            expand = expansion(mult = c(0, 0))
        ) +
        labs(
            x = "Local richness", y = expression(LCBD~"\u00d7"~M)
        ) +
        theme_classic(base_size = 15, base_family = "sans") +
        theme(
            panel.background = element_rect(
                fill = "white", colour = NA
            ),
            panel.border = element_rect(
                fill = NA, colour = "black", linewidth = 0.98
            ),
            axis.line = element_blank(),
            axis.title = element_text(size = 16),
            axis.title.y = element_text(size = 16, margin = margin(r = 18)),
            axis.text = element_text(size = 14, colour = "grey25"),
            legend.position = c(0.04, 0.04),
            legend.justification = c(0, 0),
            legend.direction = "vertical",
            legend.text = element_text(size = 12),
            legend.key.width = grid::unit(0.70, "cm"),
            legend.background = element_blank(),
            plot.background = element_rect(fill = NA, colour = NA),
            plot.margin = margin(5, 6, 5, 5)
        )
}

# Draw the BCR species-representation panels for one framework.
make_partition_panel <- function(framework) {
    framework_label <- if (framework == "coverage") "Coverage" else "Sorensen"
    partition_levels <- c(
        "Raw LCBD only", "Shared by both", "LCBD deviation only",
        "Selected by neither"
    )
    panel_levels <- c("Species pool", "PIF RI = 1", "USFWS breeding BCC")

    species_partition <- species_result$partition_by_bcr %>%
        filter(
            framework == framework_label,
            species_pool %in% c("Species pool", "Regional priority")
        ) %>%
        mutate(panel = recode(
            species_pool,
            `Species pool` = "Species pool",
            `Regional priority` = "PIF RI = 1"
        )) %>%
        select(BCR, panel, partition, percent, pool_size)
    bcc_partition <- bcc_result$partition_by_bcr %>%
        filter(
            framework == framework_label,
            bcc_definition == "Breeding BCC", pool_size > 0
        ) %>%
        transmute(
            BCR, panel = "USFWS breeding BCC", partition,
            percent, pool_size
        )
    bcr_levels <- species_partition %>%
        filter(
            panel == "Species pool",
            partition == "LCBD deviation only"
        ) %>%
        group_by(BCR) %>%
        summarise(deviation_only_percent = sum(percent), .groups = "drop") %>%
        arrange(desc(deviation_only_percent), BCR) %>%
        pull(BCR)
    data <- bind_rows(species_partition, bcc_partition) %>%
        mutate(
            BCR = factor(BCR, levels = bcr_levels),
            panel = factor(panel, levels = panel_levels),
            partition = factor(partition, levels = partition_levels)
        )
    panel_tags <- data.frame(
        panel = factor(panel_levels, levels = panel_levels),
        label = c("b", "c", "d")
    )
    bcr10_position <- match(10L, bcr_levels)
    bcr10_frames <- data.frame(
        panel = factor(panel_levels, levels = panel_levels),
        xmin = bcr10_position - 0.34,
        xmax = bcr10_position + 0.34,
        ymin = 0,
        ymax = 100
    )
    partition_palette <- c(
        "Raw LCBD only" = raw_colour,
        "Shared by both" = "#806CA0",
        "LCBD deviation only" = deviation_colour,
        "Selected by neither" = not_selected_colour
    )

    ggplot(data, aes(BCR, percent, fill = partition)) +
        geom_col(width = 0.68, position = position_stack(reverse = TRUE)) +
        geom_text(
            data = panel_tags,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = -0.35, size = 5.6, fontface = "bold",
            inherit.aes = FALSE
        ) +
        geom_rect(
            data = bcr10_frames,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "grey48", linetype = "dashed",
            linewidth = 0.55,
            inherit.aes = FALSE
        ) +
        scale_fill_manual(
            values = partition_palette, breaks = partition_levels,
            labels = c(
                "Raw LCBD only", "Shared",
                "LCBD deviation only", "Neither"
            ),
            name = NULL, drop = FALSE
        ) +
        facet_grid(~panel, drop = FALSE) +
        scale_y_continuous(
            breaks = c(0, 50, 100),
            labels = function(x) paste0(x, "%"),
            expand = expansion(mult = c(0, 0))
        ) +
        scale_x_discrete(drop = FALSE) +
        coord_cartesian(ylim = c(0, 100), clip = "off") +
        labs(x = "Bird Conservation Region", y = "Percentage of species") +
        theme_classic(base_size = 15, base_family = "sans") +
        theme(
            panel.border = element_rect(
                fill = NA, colour = "black", linewidth = 0.48
            ),
            axis.line = element_blank(),
            axis.title = element_text(size = 15),
            axis.title.y = element_text(size = 15, margin = margin(r = 3)),
            axis.text.x = element_text(
                size = 8.4, angle = 90, vjust = 0.5, hjust = 1,
                colour = "grey25"
            ),
            axis.text.y = element_text(size = 12, colour = "grey25"),
            strip.background = element_blank(),
            strip.text = element_text(size = 14.5, face = "bold"),
            panel.spacing = grid::unit(0.16, "cm"),
            legend.position = "bottom",
            legend.text = element_text(size = 11.5),
            legend.key.width = grid::unit(0.45, "cm"),
            legend.key.height = grid::unit(0.38, "cm"),
            legend.spacing.x = grid::unit(0.10, "cm"),
            legend.box.spacing = grid::unit(2, "pt"),
            legend.margin = margin(t = -1),
            plot.margin = margin(5, 5, 2, 5)
        ) +
        guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

# 6. Assemble and export Figure 5 and Figure S4 --------------------------

# Export one Figure 5 framework in all manuscript formats.
save_figure <- function(plot, stem) {
    ggsave(
        paste0(stem, ".png"), plot, width = figure_width,
        height = figure_height, dpi = 600, bg = "white"
    )
    ragg::agg_tiff(
        paste0(stem, ".tiff"), width = figure_width,
        height = figure_height, units = "in", res = 600,
        background = "white", compression = "lzw"
    )
    print(plot)
    dev.off()
    grDevices::pdf(
        paste0(stem, ".pdf"), width = figure_width,
        height = figure_height,
        family = "Helvetica", useDingbats = FALSE, bg = "white"
    )
    print(plot)
    grDevices::dev.off()
}

legend <- make_legend(text_size = 12.5, compact = FALSE)
focus_boundary <- bcr_boundaries %>% filter(BCR == focus_bcr)

for (framework in c("coverage", "sorensen")) {
    map <- make_map(framework, NULL, focus_boundary)
    fit_plot <- make_bbs_fit_panel(framework)
    partition_plot <- make_partition_panel(framework)
    map_tag <- "a"
    top_panel <- ggdraw() +
        draw_plot(map, 0.39, 0.05, 0.60, 0.89) +
        draw_grob(
            grid::rectGrob(gp = grid::gpar(fill = "grey94", col = NA)),
            0.0005, 0.395, 0.3841, 0.545
        ) +
        draw_plot(fit_plot, 0.0005, 0.405, 0.372, 0.525) +
        draw_plot(legend, 0.07, 0.01, 0.29, 0.38) +
        draw_label(map_tag, 0.985, 0.975, hjust = 1, vjust = 1,
            size = 15, fontface = "bold")
    figure <- plot_grid(
        top_panel, partition_plot, ncol = 1,
        rel_heights = c(1.28, 0.72)
    )
    figure_id <- if (framework == "coverage") "Figure5" else "FigureS4"
    stem <- file.path(
        output_dir,
        paste0(
            figure_id, "_BBS_2024_", tools::toTitleCase(framework),
            "_top5_map_top"
        )
    )
    save_figure(figure, stem)
}
