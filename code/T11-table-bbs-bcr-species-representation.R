# Supplementary table: BBS species representation at the top 5% threshold

suppressPackageStartupMessages(library(dplyr))
source("code/01-lcbd-models.R")

# 1. Input files ---------------------------------------------------------

matrix_file <- "data/BBS_2024_species_matrix.rds"
fit_file <- "process data/Figure5_BBS_annealed_expectations.rds"
species_file <- "process data/Figure5_BBS_species_priority_two_frameworks.rds"
conservation_file <- "data/BBS_conservation_status_by_BCR.xlsx"
output_file <- "process data/Table_BBS_BCR_species_representation_top5.csv"


# 2. Read BBS objects ----------------------------------------------------

bbs_matrix <- readRDS(matrix_file)
bbs_fit <- readRDS(fit_file)
bbs_species <- readRDS(species_file)
conservation <- read_xlsx_data(conservation_file)

communities <- bbs_matrix$communities
sites <- bbs_fit$site_data
traits <- bbs_species$species_traits %>%
    filter(framework == "Coverage")

bcr_names <- c(
    `5` = "Northern Pacific Rainforest",
    `6` = "Boreal Plains",
    `8` = "Boreal Softwood Shield",
    `9` = "Great Basin",
    `10` = "Northern Rockies",
    `11` = "Prairie Potholes",
    `12` = "Boreal Hardwood Transition",
    `13` = "Lower Great Lakes/St. Lawrence Plain",
    `14` = "Atlantic Northern Forest",
    `16` = "Southern Rockies/Colorado Plateau",
    `17` = "Badlands and Prairies",
    `18` = "Shortgrass Prairie",
    `19` = "Central Mixed Grass Prairie",
    `21` = "Oaks and Prairies",
    `22` = "Eastern Tallgrass Prairie",
    `23` = "Prairie Hardwood Transition",
    `24` = "Central Hardwoods",
    `25` = "West Gulf Coastal Plain/Ouachitas",
    `26` = "Mississippi Alluvial Valley",
    `27` = "Southeastern Coastal Plain",
    `28` = "Appalachian Mountains",
    `29` = "Piedmont",
    `30` = "New England/Mid-Atlantic Coast",
    `31` = "Peninsular Florida",
    `32` = "Coastal California",
    `33` = "Sonoran and Mojave Deserts",
    `35` = "Chihuahuan Desert",
    `37` = "Gulf Coastal Prairie"
)


# 3. Count each BCR-level species set -----------------------------------

count_species_sets <- function(bcr_now) {
    routes <- sites[sites$BCR == bcr_now, ]
    community <- communities[[as.character(bcr_now)]][
        routes$site_id, , drop = FALSE
    ]
    local_pool <- colSums(community) > 0

    traits_now <- traits[traits$BCR == bcr_now, ]
    traits_now <- traits_now[match(colnames(community), traits_now$species_id), ]
    pif_pool <- local_pool & traits_now$regional_priority

    status <- conservation[conservation$BCR == bcr_now, ]
    status <- status[
        match(as.integer(colnames(community)), status$AOU),
    ]
    bcc_available <- any(!is.na(status$USFWS_Breeding_BCC))
    bcc_pool <- local_pool & !is.na(status$USFWS_Breeding_BCC) &
        status$USFWS_Breeding_BCC == 1

    data.frame(
        BCR = bcr_now,
        regional_pool_n = sum(local_pool),
        pif_ri1_n = sum(pif_pool),
        usfws_bcc_n = if (bcc_available) sum(bcc_pool) else NA_integer_
    )
}

species_counts <- bind_rows(lapply(
    bbs_fit$settings$retained_bcr,
    count_species_sets
))


# 4. Calculate top 5% representation in both LCBD frameworks ------------

represented_percent <- function(community, selected, pool) {
    if (!any(pool)) return(NA_real_)
    represented <- colSums(community[selected, , drop = FALSE]) > 0
    100 * sum(represented & pool) / sum(pool)
}

# Summarize top-5% species representation for one BCR and framework.
summarise_framework <- function(bcr_now, framework_now) {
    routes <- sites[sites$BCR == bcr_now, ]
    community <- communities[[as.character(bcr_now)]][
        routes$site_id, , drop = FALSE
    ]

    traits_now <- traits[traits$BCR == bcr_now, ]
    traits_now <- traits_now[match(colnames(community), traits_now$species_id), ]
    local_pool <- colSums(community) > 0
    pif_pool <- local_pool & traits_now$regional_priority

    status <- conservation[conservation$BCR == bcr_now, ]
    status <- status[
        match(as.integer(colnames(community)), status$AOU),
    ]
    bcc_pool <- local_pool & !is.na(status$USFWS_Breeding_BCC) &
        status$USFWS_Breeding_BCC == 1

    raw_selected <- routes[[paste0(framework_now, "_raw_selected")]]
    deviation_column <- if (framework_now == "coverage") {
        "coverage_selected"
    } else {
        "sorensen_selected"
    }
    deviation_selected <- routes[[deviation_column]]

    data.frame(
        BCR = bcr_now,
        framework = if (framework_now == "coverage") "Coverage deficit" else "Sørensen",
        regional_pool_raw_percent = represented_percent(
            community, raw_selected, local_pool
        ),
        regional_pool_deviation_percent = represented_percent(
            community, deviation_selected, local_pool
        ),
        pif_ri1_raw_percent = represented_percent(
            community, raw_selected, pif_pool
        ),
        pif_ri1_deviation_percent = represented_percent(
            community, deviation_selected, pif_pool
        ),
        usfws_bcc_raw_percent = represented_percent(
            community, raw_selected, bcc_pool
        ),
        usfws_bcc_deviation_percent = represented_percent(
            community, deviation_selected, bcc_pool
        )
    )
}

top5_results <- bind_rows(lapply(
    bbs_fit$settings$retained_bcr,
    function(bcr_now) bind_rows(
        summarise_framework(bcr_now, "coverage"),
        summarise_framework(bcr_now, "sorensen")
    )
))

table_data <- top5_results %>%
    left_join(species_counts, by = "BCR") %>%
    mutate(BCR_name = unname(bcr_names[as.character(BCR)])) %>%
    select(
        BCR, BCR_name, framework,
        regional_pool_n, pif_ri1_n, usfws_bcc_n,
        regional_pool_raw_percent, regional_pool_deviation_percent,
        pif_ri1_raw_percent, pif_ri1_deviation_percent,
        usfws_bcc_raw_percent, usfws_bcc_deviation_percent
    ) %>%
    arrange(BCR, factor(framework, c("Coverage deficit", "Sørensen")))

write.csv(table_data, output_file, row.names = FALSE, na = "")
