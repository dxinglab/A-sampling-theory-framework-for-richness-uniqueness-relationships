# Prepare all fitted expectations and comparison data used by Figure 5.

# 1. Paths and settings ---------------------------------------------------

bbs_file <- "data/BBS_2024_species_matrix.rds"
output_file <- "process data/Figure5_BBS_annealed_expectations.rds"

wdpa_polygon_files <- file.path(
    "data",
    paste0("WDPA_Mar2026_Public_shp_", 0:2),
    "WDPA_Mar2026_Public_shp-polygons.shp"
)

top_proportion <- 0.05
x_bounds <- c(0.5, 1 - 1e-10)
c_bounds <- c(1e-12, 1e6)
x_starts <- c(0.9, 0.99, 0.999, 0.9999, 0.99999)
c_starts <- c(0.01, 0.1, 1, 10, 1e3, 1e5)

region_bcr <- list(
    Western = c(4, 5, 9, 10, 15, 16, 32, 33, 34),
    Central = c(6, 11, 17, 18, 19, 20, 21, 35, 36, 37),
    Eastern = c(7, 8, 12, 13, 14, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)
)

# 2. Load the common model and LCBD functions ----------------------------

source("code/01-lcbd-models.R")

# 3. Define the deterministic ranking rule -------------------------------

rank_exact <- function(value, id) {
    order_id <- order(-value, id, na.last = TRUE)
    rank <- rep(NA_integer_, length(value))
    rank[order_id] <- seq_along(order_id)
    rank
}

# 4. Read the 28-BCR sample ----------------------------------------------

bbs <- readRDS(bbs_file)
communities <- bbs$communities
routes_by_bcr <- bbs$site_data
retained_bcr <- as.integer(names(communities))

for (bcr in retained_bcr) {
    key <- as.character(bcr)
    region <- names(region_bcr)[vapply(
        region_bcr, function(x) bcr %in% x, logical(1)
    )]
    routes_by_bcr[[key]]$region <- region
}
route <- do.call(rbind, routes_by_bcr)
rownames(route) <- NULL

# 5. Fit both frameworks and rank routes within each BCR -----------------

parameter_rows <- vector("list", 2L * length(retained_bcr))
curve_rows <- vector("list", 2L * length(retained_bcr))
site_rows <- vector("list", length(retained_bcr))
bcr10_fit_panel <- list()

for (i in seq_along(retained_bcr)) {
    bcr <- retained_bcr[i]
    key <- as.character(bcr)
    matrix_bcr <- communities[[key]]
    route_bcr <- routes_by_bcr[[key]]
    matrix_bcr <- matrix_bcr[, colSums(matrix_bcr) > 0, drop = FALSE]
    M <- nrow(matrix_bcr)
    S <- ncol(matrix_bcr)
    richness <- rowSums(matrix_bcr > 0)
    site_id <- rownames(matrix_bcr)
    top_n <- ceiling(top_proportion * M)

    observed <- list(
        Coverage = coverage_lcbd(matrix_bcr),
        Sorensen = sorensen_lcbd(matrix_bcr)
    )
    site <- data.frame(
        BCR = bcr, region = route_bcr$region[1L],
        site_id = site_id, richness = richness, M = M, S = S,
        top_n_within_BCR = top_n
    )

    for (j in seq_along(observed)) {
        framework <- names(observed)[j]
        data <- data.frame(richness, lcbd = observed[[framework]])
        prefix <- tolower(framework)
        fit <- fit_annealed_xc(
            data, S, M, prefix,
            x_bounds, c_bounds, x_starts, c_starts,
            curve_richness = sort(unique(richness))
        )
        expected <- fit$curve$expected[match(richness, fit$curve$richness)]
        deviation <- observed[[framework]] - expected
        raw_rank <- rank_exact(observed[[framework]], site_id)
        deviation_rank <- rank_exact(deviation, site_id)

        site[[paste0(prefix, "_lcbd")]] <- observed[[framework]]
        site[[paste0(prefix, "_expected")]] <- expected
        site[[paste0(prefix, "_deviation")]] <- deviation
        site[[paste0(prefix, "_raw_rank")]] <- raw_rank
        site[[paste0(prefix, "_deviation_rank")]] <- deviation_rank
        site[[paste0(prefix, "_raw_selected")]] <- raw_rank <= top_n
        site[[paste0(prefix, "_selected")]] <- deviation_rank <= top_n

        row_id <- 2L * i - 2L + j
        parameter_rows[[row_id]] <- data.frame(
            BCR = bcr, region = site$region[1L], framework = framework,
            M = M, S = S, x = fit$x, c = fit$c, mse = fit$mse,
            convergence = fit$convergence, boundary = fit$boundary
        )
        curve_rows[[row_id]] <- transform(
            fit$curve, BCR = bcr, framework = framework
        )

        if (bcr == 10L) {
            bcr10_fit_panel[[framework]] <- list(
                observed = data,
                parameterized = fit_parameterized_x(
                    data, S, M, prefix, x_bounds
                ),
                annealed = fit
            )
        }
    }
    site$selection_group <- factor(
        ifelse(
            site$coverage_selected & site$sorensen_selected,
            "Shared by both",
            ifelse(
                site$coverage_selected, "Coverage only",
                ifelse(site$sorensen_selected, "Sorensen only", "Neither")
            )
        ),
        levels = c("Coverage only", "Shared by both", "Sorensen only", "Neither")
    )
    site_rows[[i]] <- site
}

parameters <- do.call(rbind, parameter_rows)
curves <- do.call(rbind, curve_rows)
site_data <- do.call(rbind, site_rows)
rownames(site_data) <- NULL

# 6. Match route starts to protected-area polygons -----------------------

route_metadata <- route[, c(
    "surveyID", "RouteDataID", "Year", "StateNum", "Route",
    "Latitude", "Longitude"
)]

route_points <- sf::st_as_sf(
    route_metadata,
    coords = c("Longitude", "Latitude"),
    crs = 4326,
    remove = FALSE
)
route_points <- sf::st_transform(route_points, 5070)
point_in_wdpa_polygon <- rep(FALSE, nrow(route_points))

for (path in wdpa_polygon_files) {
    polygons <- sf::st_read(path, quiet = TRUE)
    polygons <- polygons[
        !polygons$REALM %in% "Marine" &
            polygons$STATUS %in% c("Designated", "Established"),
    ]
    polygons <- sf::st_make_valid(sf::st_transform(polygons, 5070))
    point_in_wdpa_polygon <- point_in_wdpa_polygon |
        lengths(sf::st_intersects(route_points, polygons)) > 0L
}

route_fields <- data.frame(
    surveyID = route_metadata$surveyID,
    point_in_wdpa_polygon = as.integer(point_in_wdpa_polygon)
)

site_order <- site_data$site_id
site_data <- merge(
    site_data, route_metadata, by.x = "site_id", by.y = "surveyID",
    all.x = TRUE, sort = FALSE
)
site_data <- merge(
    site_data, route_fields, by.x = "site_id", by.y = "surveyID",
    all.x = TRUE, sort = FALSE
)
site_data <- site_data[match(site_order, site_data$site_id), ]
rownames(site_data) <- NULL

# 7. Compare framework rankings within BCR -------------------------------

bcr_comparison <- do.call(rbind, lapply(split(site_data, site_data$BCR), function(z) {
    shared <- sum(z$coverage_selected & z$sorensen_selected)
    top_n <- z$top_n_within_BCR[1L]
    data.frame(
        BCR = z$BCR[1L], region = z$region[1L], M = nrow(z), top_n = top_n,
        n_shared = shared, overlap_each_percent = 100 * shared / top_n,
        raw_spearman = cor(
            z$coverage_lcbd, z$sorensen_lcbd, method = "spearman"
        ),
        deviation_spearman = cor(
            z$coverage_deviation, z$sorensen_deviation, method = "spearman"
        )
    )
}))
rownames(bcr_comparison) <- NULL
bcr_comparison$spearman_difference <-
    bcr_comparison$deviation_spearman - bcr_comparison$raw_spearman
bcr_comparison$absolute_difference <- abs(bcr_comparison$spearman_difference)

paired_test <- wilcox.test(
    bcr_comparison$deviation_spearman,
    bcr_comparison$raw_spearman,
    paired = TRUE, exact = FALSE, conf.int = TRUE
)
sample_size_test <- cor.test(
    bcr_comparison$M,
    bcr_comparison$absolute_difference,
    method = "spearman", exact = FALSE
)
framework_agreement_test <- data.frame(
    n_BCR = nrow(bcr_comparison),
    median_difference = median(bcr_comparison$spearman_difference),
    confidence_lower = paired_test$conf.int[1L],
    confidence_upper = paired_test$conf.int[2L],
    paired_statistic = unname(paired_test$statistic),
    paired_p = paired_test$p.value,
    sample_size_rho = unname(sample_size_test$estimate),
    sample_size_p = sample_size_test$p.value
)

n_selected_coverage <- sum(site_data$coverage_selected)
n_selected_sorensen <- sum(site_data$sorensen_selected)
n_shared <- sum(site_data$coverage_selected & site_data$sorensen_selected)
comparison <- data.frame(
    n_routes = nrow(site_data), n_bcr = length(retained_bcr),
    top_proportion = top_proportion,
    n_selected_coverage = n_selected_coverage,
    n_selected_sorensen = n_selected_sorensen,
    n_shared = n_shared,
    overlap_coverage_percent = 100 * n_shared / n_selected_coverage,
    overlap_sorensen_percent = 100 * n_shared / n_selected_sorensen
)

selection_summary <- aggregate(
    cbind(richness, coverage_deviation, sorensen_deviation) ~ selection_group,
    site_data, mean
)
selection_summary$n_routes <- as.integer(table(site_data$selection_group)[
    as.character(selection_summary$selection_group)
])

raw_deviation_overlap <- do.call(rbind, lapply(c("coverage", "sorensen"), function(framework) {
    raw_name <- paste0(framework, "_raw_selected")
    deviation_name <- paste0(framework, "_selected")
    do.call(rbind, lapply(split(site_data, site_data$BCR), function(z) {
        raw <- z[[raw_name]]
        deviation <- z[[deviation_name]]
        shared <- sum(raw & deviation)
        top_n <- z$top_n_within_BCR[1L]
        data.frame(
            framework = tools::toTitleCase(framework),
            BCR = z$BCR[1L], region = z$region[1L], M = nrow(z),
            top_n = top_n, n_shared = shared,
            overlap_percent = 100 * shared / top_n,
            expected_overlap = top_n^2 / nrow(z),
            adjusted_overlap_percent = 100 *
                (shared / top_n - top_n / nrow(z)) /
                (1 - top_n / nrow(z)),
            raw_only_n = sum(raw & !deviation),
            deviation_only_n = sum(!raw & deviation)
        )
    }))
}))
rownames(raw_deviation_overlap) <- NULL

# 8. Save one reusable result --------------------------------------------

result <- list(
    settings = list(
        model = "Model 3",
        fit_unit = "BCR",
        fit_objective = "site-level MSE",
        selection_unit = "within BCR",
        top_proportion = top_proportion,
        retained_bcr = retained_bcr,
        x_bounds = x_bounds,
        c_bounds = c_bounds,
        model_source = "code/01-lcbd-models.R",
        framework_rules = c(
            Coverage = "annealed Coverage-deficit expectation",
            Sorensen = "annealed mean-other Sorensen expectation"
        )
    ),
    parameters = parameters,
    curves = curves,
    site_data = site_data,
    comparison = comparison,
    bcr_comparison = bcr_comparison,
    framework_agreement_test = framework_agreement_test,
    raw_deviation_overlap = raw_deviation_overlap,
    selection_summary = selection_summary,
    bcr10_fit_panel = bcr10_fit_panel
)

saveRDS(result, output_file, compress = "gzip")


# Species representation -------------------------------------------------


# 9. Species-representation paths and settings ---------------------------

bbs_file <- "data/BBS_2024_species_matrix.rds"
priority_file <- "process data/Figure5_BBS_annealed_expectations.rds"
conservation_file <- "data/BBS_conservation_status_by_BCR.xlsx"
output_file <- "process data/Figure5_BBS_species_priority_two_frameworks.rds"

frameworks <- c("Coverage", "Sorensen")
partition_levels <- c(
    "Raw LCBD only", "Shared by both", "LCBD deviation only",
    "Selected by neither"
)

# 10. Define species-summary functions -----------------------------------

partition_counts <- function(raw, deviation, pool) {
    status <- ifelse(
        raw & deviation, "Shared by both",
        ifelse(
            raw, "Raw LCBD only",
            ifelse(deviation, "LCBD deviation only", "Selected by neither")
        )
    )
    total <- sum(pool)
    count <- table(factor(status[pool], levels = partition_levels))
    data.frame(
        partition = partition_levels,
        n_species = as.integer(count),
        pool_size = total,
        percent = if (total == 0L) NA_real_ else 100 * as.integer(count) / total
    )
}

# Summarise species represented by the two route rankings.
summarise_representation <- function(data, group_column) {
    result <- reshape(
        data[data$partition != "Selected by neither", ],
        idvar = c("BCR", "framework", group_column, "pool_size"),
        timevar = "partition", direction = "wide"
    )
    names(result) <- gsub("[ .]", "_", names(result))
    result$raw_selected_n <- result$n_species_Raw_LCBD_only +
        result$n_species_Shared_by_both
    result$deviation_selected_n <- result$n_species_LCBD_deviation_only +
        result$n_species_Shared_by_both
    result$raw_coverage_percent <- 100 * result$raw_selected_n / result$pool_size
    result$deviation_coverage_percent <-
        100 * result$deviation_selected_n / result$pool_size
    result$deviation_gain_percent <-
        result$deviation_coverage_percent - result$raw_coverage_percent
    result
}

# 11. Read route priorities and conservation status ----------------------

priority <- readRDS(priority_file)

bbs <- readRDS(bbs_file)
communities <- bbs$communities
species_id <- colnames(communities[[1L]])
site_data <- priority$site_data

conservation <- read_xlsx_data(conservation_file)
metadata_source <- conservation[!duplicated(conservation$AOU), ]
metadata_source <- metadata_source[
    match(as.integer(species_id), metadata_source$AOU),
]
species_metadata <- data.frame(
    species_id = species_id,
    AOU = as.integer(species_id),
    common_name = metadata_source$Common_Name,
    scientific_name = metadata_source$Scientific_Name,
    order = metadata_source$Order,
    family = metadata_source$Family,
    stringsAsFactors = FALSE
)

# 12. Classify species independently within each BCR ---------------------

partition_rows <- vector("list", 2L * length(priority$settings$retained_bcr))
species_rows <- vector("list", 2L * length(priority$settings$retained_bcr))
framework_rows <- vector("list", length(priority$settings$retained_bcr))
match_rows <- vector("list", length(priority$settings$retained_bcr))

for (i in seq_along(priority$settings$retained_bcr)) {
    bcr <- priority$settings$retained_bcr[i]
    routes <- site_data[site_data$BCR == bcr, ]
    matrix_bcr <- communities[[as.character(bcr)]][
        routes$site_id, , drop = FALSE
    ]
    local_pool <- colSums(matrix_bcr) > 0
    abundance <- colSums(matrix_bcr)
    occurrence <- colSums(matrix_bcr > 0)

    status_bcr <- conservation[conservation$BCR == bcr, ]
    status_index <- match(as.integer(species_id), status_bcr$AOU)
    regional_match <- status_bcr$PIF_Record_Matched[status_index] == 1
    regional_priority <- status_bcr$PIF_RI_1[status_index] == 1
    breeding_bcc <- status_bcr$USFWS_Breeding_BCC[status_index] == 1
    regional_match[is.na(regional_match)] <- FALSE
    regional_priority[is.na(regional_priority)] <- FALSE
    match_rows[[i]] <- data.frame(
        BCR = bcr,
        local_species_n = sum(local_pool),
        pif_matched_n = sum(local_pool & regional_match),
        pif_unmatched_n = sum(local_pool & !regional_match),
        regional_priority_n = sum(local_pool & regional_priority)
    )

    selected <- list()
    for (j in seq_along(frameworks)) {
        framework <- frameworks[j]
        prefix <- tolower(framework)
        raw_sites <- routes$site_id[routes[[paste0(prefix, "_raw_selected")]]]
        deviation_sites <- routes$site_id[routes[[paste0(prefix, "_selected")]]]
        raw <- colSums(matrix_bcr[raw_sites, , drop = FALSE]) > 0
        deviation <- colSums(matrix_bcr[deviation_sites, , drop = FALSE]) > 0
        selected[[framework]] <- deviation

        pools <- list(
            `Species pool` = local_pool,
            `Regional priority` = local_pool & regional_priority
        )
        partitions <- do.call(rbind, lapply(names(pools), function(pool_name) {
            data.frame(
                species_pool = pool_name,
                partition_counts(raw, deviation, pools[[pool_name]])
            )
        }))
        partitions$BCR <- bcr
        partitions$framework <- framework
        partition_rows[[2L * i - 2L + j]] <- partitions

        species_rows[[2L * i - 2L + j]] <- transform(
            species_metadata,
            BCR = bcr,
            framework = framework,
            local_abundance = abundance,
            local_occurrence = occurrence,
            local_occupancy = occurrence / nrow(matrix_bcr),
            regional_match = regional_match,
            regional_priority = regional_priority,
            usfws_breeding_bcc = breeding_bcc,
            raw_selected_species = raw,
            deviation_selected_species = deviation,
            selection_partition = ifelse(
                raw & deviation, "Shared by both",
                ifelse(
                    raw, "Raw LCBD only",
                    ifelse(deviation, "LCBD deviation only", "Selected by neither")
                )
            )
        )[local_pool, ]
    }

    coverage <- selected$Coverage
    sorensen <- selected$Sorensen
    pools <- list(
        `Species pool` = local_pool,
        `Regional priority` = local_pool & regional_priority
    )
    comparison <- do.call(rbind, lapply(names(pools), function(pool_name) {
        status <- ifelse(
            coverage & sorensen, "Shared by both",
            ifelse(
                coverage, "Coverage only",
                ifelse(sorensen, "Sorensen only", "Selected by neither")
            )
        )
        total <- sum(pools[[pool_name]])
        count <- table(factor(
            status[pools[[pool_name]]],
            levels = c(
                "Coverage only", "Shared by both", "Sorensen only",
                "Selected by neither"
            )
        ))
        data.frame(
            BCR = bcr, species_pool = pool_name,
            partition = names(count), n_species = as.integer(count),
            pool_size = total,
            percent = if (total == 0L) NA_real_ else 100 * as.integer(count) / total
        )
    }))
    framework_rows[[i]] <- comparison
}

partition_by_bcr <- do.call(rbind, partition_rows)
species_traits <- do.call(rbind, species_rows)
framework_comparison <- do.call(rbind, framework_rows)
pif_match_audit <- do.call(rbind, match_rows)
rownames(partition_by_bcr) <- NULL
rownames(species_traits) <- NULL
rownames(framework_comparison) <- NULL
rownames(pif_match_audit) <- NULL

# 13. Summarise species coverage and biological composition --------------

mean_partition <- aggregate(
    percent ~ framework + species_pool + partition,
    partition_by_bcr, mean, na.rm = TRUE
)
mean_partition$n_bcr <- aggregate(
    percent ~ framework + species_pool + partition,
    partition_by_bcr, function(x) sum(is.finite(x))
)$percent

coverage_by_bcr <- summarise_representation(partition_by_bcr, "species_pool")

metadata_coverage <- data.frame(
    n_bbs_taxa = nrow(species_metadata),
    n_taxonomy_matched = sum(
        !is.na(species_metadata$common_name) &
            nzchar(species_metadata$common_name)
    )
)

# 14. Save species-representation results --------------------------------

result <- list(
    settings = list(
        top_proportion = priority$settings$top_proportion,
        selection_unit = "within BCR",
        bcr_model_source = priority_file,
        figure2_source = priority$settings$figure2_source,
        priority_definition = "PIF Species of Regional Importance RI = 1",
        conservation_metadata = conservation_file
    ),
    metadata_coverage = metadata_coverage,
    pif_match_audit = pif_match_audit,
    species_metadata = species_metadata,
    partition_by_bcr = partition_by_bcr,
    mean_partition = mean_partition,
    coverage_by_bcr = coverage_by_bcr,
    framework_comparison = framework_comparison,
    species_traits = species_traits
)

saveRDS(result, output_file, compress = "gzip")


# USFWS Birds of Conservation Concern -----------------------------------


# 15. USFWS BCC paths and settings ---------------------------------------

species_file <- "process data/Figure5_BBS_species_priority_two_frameworks.rds"
output_file <- "process data/Figure5_BBS_USFWS_BCC_2021_comparison.rds"

bcc_definition <- "Breeding BCC"


# 16. Read species selections --------------------------------------------

species_result <- readRDS(species_file)
species <- species_result$species_traits
analysed_bcr <- sort(unique(species$BCR))
available_bcr <- sort(unique(species$BCR[
    !is.na(species$usfws_breeding_bcc)
]))
excluded_bcr <- setdiff(analysed_bcr, available_bcr)


# 17. Classify BCC taxa and compare route priorities ---------------------

rows <- vector(
    "list", length(frameworks) * length(analysed_bcr)
)
match_rows <- vector("list", length(analysed_bcr))

for (i in seq_along(analysed_bcr)) {
    bcr_now <- analysed_bcr[i]
    species_bcr <- species[
        species$BCR == bcr_now & species$framework == frameworks[1],
    ]
    match_rows[[i]] <- data.frame(
        BCR = bcr_now,
        observed_breeding_bcc_n = sum(
            species_bcr$usfws_breeding_bcc, na.rm = TRUE
        ),
        bcc_available = bcr_now %in% available_bcr
    )

    for (j in seq_along(frameworks)) {
        framework <- frameworks[j]
        data_now <- species[
            species$BCR == bcr_now & species$framework == framework,
        ]
        pool <- data_now$usfws_breeding_bcc
        pool[is.na(pool)] <- FALSE
        result_now <- partition_counts(
            data_now$raw_selected_species,
            data_now$deviation_selected_species,
            pool
        )
        result_now$BCR <- bcr_now
        result_now$framework <- framework
        result_now$bcc_definition <- bcc_definition
        rows[[length(frameworks) * (i - 1L) + j]] <- result_now
    }
}

partition_by_bcr <- do.call(rbind, rows)
match_audit <- do.call(rbind, match_rows)
rownames(partition_by_bcr) <- NULL
rownames(match_audit) <- NULL


# 18. Summarise BCC representation ---------------------------------------

valid <- partition_by_bcr$BCR %in% available_bcr &
    partition_by_bcr$pool_size > 0
valid_partition <- partition_by_bcr[valid, ]

coverage_by_bcr <- summarise_representation(valid_partition, "bcc_definition")

summary_grid <- data.frame(
    framework = frameworks,
    bcc_definition = bcc_definition,
    stringsAsFactors = FALSE
)
summary_by_framework <- do.call(rbind, lapply(seq_len(nrow(summary_grid)), function(i) {
    framework <- summary_grid$framework[i]
    definition <- summary_grid$bcc_definition[i]
    data_now <- coverage_by_bcr[
        coverage_by_bcr$framework == framework &
            coverage_by_bcr$bcc_definition == definition,
    ]
    data.frame(
        framework = framework,
        bcc_definition = definition,
        n_bcr = nrow(data_now),
        mean_raw_percent = mean(data_now$raw_coverage_percent),
        mean_deviation_percent = mean(data_now$deviation_coverage_percent),
        mean_gain_percent = mean(data_now$deviation_gain_percent),
        median_gain_percent = median(data_now$deviation_gain_percent),
        positive_gain_bcr = sum(data_now$deviation_gain_percent > 0),
        equal_gain_bcr = sum(data_now$deviation_gain_percent == 0),
        negative_gain_bcr = sum(data_now$deviation_gain_percent < 0)
    )
}))


# 19. Save the BCC comparison --------------------------------------------

result <- list(
    settings = list(
        definition = paste(
            "Taxon listed by USFWS Birds of Conservation Concern 2021",
            "for the corresponding Bird Conservation Region"
        ),
        source = conservation_file,
        source_scope = "United States BCR designations",
        population_entries_collapsed_to_species = TRUE,
        excluded_bcr = excluded_bcr
    ),
    match_audit = match_audit,
    partition_by_bcr = partition_by_bcr,
    coverage_by_bcr = coverage_by_bcr,
    summary_by_framework = summary_by_framework
)

saveRDS(result, output_file, compress = "gzip")


# Protected-area representation -----------------------------------------


# 20. Protected-route analysis settings ----------------------------------

input_file <- "process data/Figure5_BBS_annealed_expectations.rds"
n_permutations <- 99999L
seed <- 20260928L

frameworks <- c("coverage", "sorensen")

# 21. Define protected-route statistical functions -----------------------

exact_binomial_ci <- function(success, total) {
    unname(stats::binom.test(success, total)$conf.int)
}

paired_selection_test <- function(data, raw, deviation, protected, seed_now) {
    set.seed(seed_now)
    observed <- mean(protected[deviation]) - mean(protected[raw])
    null_count_difference <- numeric(n_permutations)
    for (bcr in sort(unique(data$BCR))) {
        keep <- data$BCR == bcr
        discordant <- xor(raw[keep], deviation[keep])
        if (!any(discordant)) next
        raw_only_n <- sum(raw[keep] & !deviation[keep])
        protected_n <- sum(protected[keep][discordant])
        raw_protected <- stats::rhyper(
            n_permutations,
            protected_n,
            sum(discordant) - protected_n,
            raw_only_n
        )
        null_count_difference <- null_count_difference +
            protected_n - 2 * raw_protected
    }
    null_difference <- null_count_difference / sum(raw)
    c(
        difference = observed,
        p = (1 + sum(abs(null_difference) >= abs(observed) - 1e-15)) /
            (n_permutations + 1)
    )
}

# 22. Calculate protected-route statistics -------------------------------

data <- readRDS(input_file)$site_data
protected <- as.logical(data$point_in_wdpa_polygon)
rows <- list()

for (framework_index in seq_along(frameworks)) {
    framework <- frameworks[framework_index]
    raw <- data[[paste0(framework, "_raw_selected")]]
    deviation <- data[[paste0(framework, "_selected")]]
    raw_success <- sum(protected & raw)
    deviation_success <- sum(protected & deviation)
    raw_ci <- exact_binomial_ci(raw_success, sum(raw))
    deviation_ci <- exact_binomial_ci(deviation_success, sum(deviation))
    offset <- 1000L + 10L * framework_index
    paired <- paired_selection_test(
        data, raw, deviation, protected, seed + offset + 3L
    )

    rows[[framework_index]] <- data.frame(
        protection_definition = "route_start",
        framework = if (framework == "coverage") "Coverage" else "Sørensen",
        n_BCR = length(unique(data$BCR)),
        n_routes = nrow(data),
        selected_n = sum(raw),
        raw_protected_n = raw_success,
        raw_percent = 100 * mean(protected[raw]),
        raw_ci_low_percent = 100 * raw_ci[1L],
        raw_ci_high_percent = 100 * raw_ci[2L],
        deviation_protected_n = deviation_success,
        deviation_percent = 100 * mean(protected[deviation]),
        deviation_ci_low_percent = 100 * deviation_ci[1L],
        deviation_ci_high_percent = 100 * deviation_ci[2L],
        deviation_minus_raw_percentage_points = 100 * paired["difference"],
        paired_p = paired["p"]
    )
}

statistics <- do.call(rbind, rows)
priority_result <- readRDS(input_file)
priority_result$protection_statistics <- statistics
saveRDS(priority_result, input_file, compress = "gzip")
