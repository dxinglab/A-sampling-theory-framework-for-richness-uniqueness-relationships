# Generate the spatially explicit neutral community used in Figure 2.

suppressPackageStartupMessages(library(rcoalescence))

output_file <- paste0(
    "data/",
    "spatially_explicit_sigma26_nu1p77e4.rds"
)
output_directory <- tempdir()
landscape_width <- 650L
speciation_rate <- 1.77e-4
dispersal_scale <- 26

simulation <- SpatialTreeSimulation$new()
simulation$setSimulationParameters(
    task = 1,
    seed = dispersal_scale,
    output_directory = output_directory,
    min_speciation_rate = speciation_rate,
    fine_map_file = "null",
    fine_map_x_size = landscape_width,
    fine_map_y_size = landscape_width / 2,
    deme = 1,
    sigma = dispersal_scale,
    landscape_type = "closed"
)
simulation$runSimulation()
simulation$applySpeciationRates(
    speciation_rates = speciation_rate,
    use_spatial = TRUE
)
simulation$output()

community <- simulation$getSpeciesLocations(1L)
names(community) <- c("sp", "gx", "gy")
community$gx <- community$gx / landscape_width * 1000
community$gy <- community$gy / landscape_width * 1000

saveRDS(community, output_file, version = 2)
