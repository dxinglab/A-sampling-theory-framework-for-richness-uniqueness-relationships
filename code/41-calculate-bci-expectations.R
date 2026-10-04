# Fit the BCI annealed and parameterized expectations used by Figure 4.

suppressPackageStartupMessages(library(Matrix))

# 1. Paths and settings ---------------------------------------------------

input_file <- "data/BCI_2015_species_matrices.rds"
output_file <- "process data/Empirical_annealed_expectations.rds"

x_bounds <- c(0.5, 1 - 1e-10)
c_bounds <- c(1e-12, 1e6)
x_starts <- c(0.9, 0.99, 0.999, 0.9999, 0.99999)
c_starts <- c(0.01, 0.1, 1, 10, 1e3, 1e5)

source("code/01-lcbd-models.R")

# 2. Calculate the observed BCI scores -----------------------------------

bci <- readRDS(input_file)
grains <- c("5", "20", "100")
stopifnot(all(vapply(
    grains,
    function(grain) identical(
        names(bci$regional_abundance),
        colnames(bci$communities[[grain]])
    ),
    logical(1)
)))
observed_data <- setNames(lapply(grains, function(grain) {
    community <- bci$communities[[grain]]
    data.frame(
        richness = rowSums(community > 0),
        coverage = coverage_lcbd(community, bci$regional_abundance),
        sorensen = sorensen_lcbd(community)
    )
}), grains)

# 3. Fit the three BCI grains --------------------------------------------

frameworks <- list()
parameter_rows <- list()
k <- 0L

for (framework in c("coverage", "sorensen")) {
    observed <- setNames(lapply(grains, function(grain) {
        data <- observed_data[[grain]]
        data.frame(
            richness = data$richness,
            lcbd = data[[framework]]
        )
    }), grains)
    fits <- list()

    for (grain in names(observed)) {
        data <- observed[[grain]]
        S <- ncol(bci$communities[[grain]])
        M <- nrow(data)
        fitted_x_result <- fit_parameterized_x(
            data, S, M, framework, x_bounds
        )
        model3_result <- fit_annealed_xc(
            data, S, M, framework,
            x_bounds, c_bounds, x_starts, c_starts
        )
        fitted_x <- list(
            curve = fitted_x_result$curve,
            parameter = data.frame(
                model = "Parameterized",
                x = fitted_x_result$x,
                c = NA_real_,
                mse = fitted_x_result$mse
            )
        )
        model3 <- list(
            curve = model3_result$curve,
            parameter = data.frame(
                model = "Annealed",
                x = model3_result$x,
                c = model3_result$c,
                y = model3_result$y,
                mse = model3_result$mse,
                convergence = model3_result$convergence,
                boundary = model3_result$boundary
            )
        )
        fits[[grain]] <- list(model3 = model3, fitted_x = fitted_x)
        k <- k + 1L
        parameter_rows[[k]] <- transform(
            model3$parameter,
            framework = framework, panel = "bci", grain = grain,
            S = S, M = M, fit_objective = "site-level MSE"
        )
    }

    frameworks[[framework]] <- list(
        observed = list(bci = observed),
        fits = list(bci = fits)
    )
}

# 4. Compare framework rankings -----------------------------------------

framework_agreement <- do.call(rbind, lapply(grains, function(grain) {
    coverage <- frameworks$coverage$observed$bci[[grain]]
    sorensen <- frameworks$sorensen$observed$bci[[grain]]
    coverage_curve <- frameworks$coverage$fits$bci[[grain]]$model3$curve
    sorensen_curve <- frameworks$sorensen$fits$bci[[grain]]$model3$curve
    coverage_deviation <- coverage$lcbd - coverage_curve$expected[
        match(coverage$richness, coverage_curve$richness)
    ]
    sorensen_deviation <- sorensen$lcbd - sorensen_curve$expected[
        match(sorensen$richness, sorensen_curve$richness)
    ]
    data.frame(
        grain = grain,
        n_sites = nrow(coverage),
        raw_spearman = cor(coverage$lcbd, sorensen$lcbd, method = "spearman"),
        deviation_spearman = cor(
            coverage_deviation, sorensen_deviation, method = "spearman"
        )
    )
}))

# 5. Save the fitted BCI object ------------------------------------------

saveRDS(
    list(
        settings = list(
            x_bounds = x_bounds, c_bounds = c_bounds,
            fit_objective = "site-level MSE",
            sorensen_scale = "mean dissimilarity to the other M - 1 sites"
        ),
        frameworks = frameworks,
        parameters = do.call(rbind, parameter_rows),
        framework_agreement = framework_agreement
    ),
    output_file,
    version = 2
)
