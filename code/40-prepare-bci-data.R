# Prepare BCI community matrices from the eighth forest census.

suppressPackageStartupMessages(library(Matrix))

# 1. Files and grains ----------------------------------------------------

input_file <- "data/bci.tree8.rdata"
output_file <- "data/BCI_2015_species_matrices.rds"
grains <- c("5" = 5L, "20" = 20L, "100" = 100L)


# 2. Living trees --------------------------------------------------------

input <- new.env(parent = emptyenv())
load(input_file, envir = input)

trees <- input$bci.tree8
trees <- trees[
    trees$status == "A" &
        !is.na(trees$sp) &
        !is.na(trees$gx) &
        !is.na(trees$gy),
    c("sp", "gx", "gy")
]

species <- sort(unique(trees$sp), method = "radix")


# 3. Quadrat-by-species matrices ----------------------------------------

make_community <- function(grain) {
    cells_y <- 500L %/% grain
    cell_x <- floor(trees$gx / grain)
    cell_y <- floor(trees$gy / grain)

    if (grain == 5L) {
        site_label <- paste(cell_x, cell_y, sep = "_")
        site_levels <- sort(unique(site_label), method = "radix")
        parts <- do.call(rbind, strsplit(site_levels, "_", fixed = TRUE))
        source_levels <- as.character(
            as.integer(parts[, 1L]) * cells_y + as.integer(parts[, 2L])
        )
        site_index <- match(site_label, site_levels)
    } else {
        source_site_id <- cell_x * cells_y + cell_y
        source_levels <- as.character(sort(unique(source_site_id)))
        site_index <- match(as.character(source_site_id), source_levels)
    }

    sparseMatrix(
        i = site_index,
        j = match(trees$sp, species),
        x = 1,
        dims = c(length(source_levels), length(species)),
        dimnames = list(source_levels, species)
    )
}

communities <- lapply(grains, make_community)


# 4. Quadrat coordinates and regional abundances ------------------------

make_site_data <- function(community, grain) {
    source_site_id <- as.integer(rownames(community))
    cells_y <- 500L %/% grain

    data.frame(
        site_id = seq_len(nrow(community)),
        source_site_id = rownames(community),
        coord_x = floor(source_site_id / cells_y) * grain,
        coord_y = (source_site_id %% cells_y) * grain
    )
}

site_data <- Map(make_site_data, communities, grains)
regional_abundance <- as.numeric(colSums(communities[[1L]]))
names(regional_abundance) <- colnames(communities[[1L]])

bci <- list(
    census_year = 2015L,
    grains = grains,
    communities = communities,
    site_data = site_data,
    regional_abundance = regional_abundance
)


# 5. Save ----------------------------------------------------------------

saveRDS(bci, output_file, compress = "xz")
