# Prepare all fitted controls and expectations used by Figure 2.

suppressPackageStartupMessages(library(Matrix))

# Load LCBD and neutral-expectation functions.
source("code/01-lcbd-models.R")

# 1. Paths and settings ---------------------------------------------------

explicit_file <- "data/spatially_explicit_sigma26_nu1p77e4.rds"
output_dir <- "process data/spatially_implicit_grain_matched_c"
matrix_file <- file.path(output_dir, "maintext_explicit_grain_matrices.rds")
estimate_file <- file.path(output_dir, "grain_matched_c_estimates.rds")
log_c_bounds <- log(c(1e-8, 1e6))
grains <- c(`5` = 5L, `20` = 20L, `100` = 100L)

# 2. Reconstruct the three explicit abundance matrices ------------------

individuals <- readRDS(explicit_file)
species_levels <- sort(unique(individuals$sp))

make_community <- function(grain) {
    site <- paste(
        floor(individuals$gx / grain),
        floor(individuals$gy / grain),
        sep = "_"
    )
    site_levels <- sort(unique(site), method = "radix")
    sparseMatrix(
        i = match(site, site_levels),
        j = match(individuals$sp, species_levels),
        x = 1,
        dims = c(length(site_levels), length(species_levels)),
        dimnames = list(site_levels, as.character(species_levels))
    )
}

communities <- lapply(grains, make_community)
regional_sad <- as.numeric(colSums(communities[[1L]]))
names(regional_sad) <- colnames(communities[[1L]])

saveRDS(
    list(
        source = basename(explicit_file),
        seed = 26L,
        sigma = 26,
        speciation_rate = 1.77e-4,
        grains = grains,
        regional_sad = regional_sad,
        communities = communities
    ),
    matrix_file,
    version = 2
)

# 3. Fit c with d fixed exactly to one -----------------------------------

prepare_nb_counts <- function(community) {
    cells <- summary(community)
    nonzero <- split(
        cells$x,
        factor(cells$j, levels = seq_len(ncol(community)))
    )
    list(
        M = nrow(community),
        mu = as.numeric(colSums(community)) / nrow(community),
        nonzero = nonzero
    )
}

negative_log_likelihood <- function(log_c, counts) {
    c_value <- exp(log_c)
    log_likelihood <- 0
    for (i in seq_along(counts$mu)) {
        observed <- counts$nonzero[[i]]
        size <- c_value * counts$mu[i]
        log_likelihood <- log_likelihood +
            (counts$M - length(observed)) * dnbinom(
                0, mu = counts$mu[i], size = size, log = TRUE
            ) +
            sum(dnbinom(
                observed, mu = counts$mu[i], size = size, log = TRUE
            ))
    }
    -log_likelihood
}

fit_c <- function(community) {
    counts <- prepare_nb_counts(community)
    optimum <- optimize(
        negative_log_likelihood,
        interval = log_c_bounds,
        counts = counts,
        tol = 1e-9
    )
    candidates <- c(log_c_bounds, optimum$minimum)
    losses <- vapply(
        candidates,
        negative_log_likelihood,
        numeric(1),
        counts = counts
    )
    best <- which.min(losses)
    list(
        c = exp(candidates[best]),
        negative_log_likelihood = losses[best],
        boundary = best < 3L
    )
}

fits <- lapply(communities, fit_c)
regional_x <- uniroot(
    function(x) {
        -length(regional_sad) / log1p(-x) * x / (1 - x) -
            sum(regional_sad)
    },
    interval = c(0.5, 1 - 1e-12), tol = 1e-12
)$root

replicate_table <- do.call(rbind, lapply(names(communities), function(grain) {
    M <- nrow(communities[[grain]])
    c_value <- fits[[grain]]$c
    log_y <- (c_value / M) * log1p(1 / c_value)
    data.frame(
        seed = 26L,
        sigma = 26,
        speciation_rate = 1.77e-4,
        grain = as.integer(grain),
        M = M,
        c = c_value,
        xprime = 1 / (1 + c_value),
        y = exp(log_y),
        M_log_y = M * log_y,
        regional_x = regional_x,
        regional_S = length(regional_sad),
        regional_N = sum(regional_sad),
        negative_log_likelihood = fits[[grain]]$negative_log_likelihood,
        boundary = fits[[grain]]$boundary,
        stringsAsFactors = FALSE
    )
}))

# One realization gives identical point estimates and quartiles.
summary_table <- data.frame(
    grain = replicate_table$grain,
    M = replicate_table$M,
    n_replicates = 1L,
    median_c = replicate_table$c,
    c_q25 = replicate_table$c,
    c_q75 = replicate_table$c,
    median_xprime = replicate_table$xprime,
    xprime_q25 = replicate_table$xprime,
    xprime_q75 = replicate_table$xprime,
    median_y = replicate_table$y,
    median_M_log_y = replicate_table$M_log_y,
    regional_x = replicate_table$regional_x,
    stringsAsFactors = FALSE
)

saveRDS(
    list(
        settings = list(d = 1, log_c_bounds = log_c_bounds),
        replicate_table = replicate_table,
        summary_table = summary_table
    ),
    estimate_file,
    version = 2
)


# 4. Spatially implicit settings -----------------------------------------

matched_file <- file.path(
    output_dir,
    "spatially_implicit_sigma26_nu1p77e4_grain_matched_c_matrices.rds"
)
base_seed <- 20260917L

explicit <- readRDS(matrix_file)
estimates <- readRDS(estimate_file)$replicate_table
regional_sad <- explicit$regional_sad

# 5. Generate independent negative-binomial local abundances ------------

draw_implicit <- function(M, c_value, seed) {
    set.seed(seed)
    rows <- vector("list", length(regional_sad))
    values <- vector("list", length(regional_sad))

    for (i in seq_along(regional_sad)) {
        mu <- regional_sad[i] / M
        abundance <- rnbinom(M, mu = mu, size = c_value * mu)
        present <- which(abundance > 0)
        rows[[i]] <- present
        values[[i]] <- abundance[present]
    }

    sparseMatrix(
        i = unlist(rows, use.names = FALSE),
        j = rep(seq_along(regional_sad), lengths(rows)),
        x = unlist(values, use.names = FALSE),
        dims = c(M, length(regional_sad)),
        dimnames = list(
            paste0("site_", seq_len(M)),
            names(regional_sad)
        )
    )
}

matched_matrices <- setNames(lapply(seq_len(nrow(estimates)), function(i) {
    draw_implicit(
        M = estimates$M[i],
        c_value = estimates$c[i],
        seed = base_seed + estimates$grain[i]
    )
}), estimates$grain)

# 6. Save the grain-matched neutral control ------------------------------

saveRDS(
    list(
        scenario = "implicit_grain_matched_c",
        source_explicit = explicit$source,
        explicit_seed = explicit$seed,
        explicit_sigma = explicit$sigma,
        explicit_speciation_rate = explicit$speciation_rate,
        base_seed = base_seed,
        d = 1,
        regional_sad = explicit$regional_sad,
        estimates = estimates,
        matrices = matched_matrices
    ),
    matched_file,
    version = 2
)


# 7. Figure 2 input paths ------------------------------------------------

output_file <- "process data/Figure2_model_input.rds"

# 8. Calculate Coverage-deficit and Sorensen scores ---------------------

calculate_scores <- function(community, site_id) {
    community <- as(community, "dgCMatrix")
    incidence <- as(community > 0, "dMatrix")
    M <- nrow(community)
    richness <- as.numeric(rowSums(incidence))
    weight <- as.numeric(colSums(community)) / sum(community)
    richness_levels <- sort(unique(richness))
    group <- match(richness, richness_levels)
    occurrence <- vapply(seq_along(richness_levels), function(i) {
        colSums(incidence[group == i, , drop = FALSE])
    }, numeric(ncol(community)))
    similarity <- 2 * incidence %*% occurrence /
        outer(richness, richness_levels, "+")
    similarity[!is.finite(similarity)] <- 0

    data.frame(
        site_id = as.character(site_id),
        richness = richness,
        coverage = M / (M - 1) *
            (1 - as.numeric(incidence %*% weight)),
        sorensen = (M - rowSums(similarity)) / (M - 1)
    )
}

# 9. Read the spatially explicit simulation -----------------------------

simulation <- readRDS(matrix_file)
simulation_scores <- lapply(simulation$communities, function(community) {
    calculate_scores(community, rownames(community))
})

# 10. Save the Figure 2 model input -------------------------------------

groups <- lapply(names(grains), function(grain) list(
    group_id = paste0("Simulation_grain", grain),
    M = as.numeric(nrow(simulation_scores[[grain]])),
    observed_data = simulation_scores[[grain]]
))

saveRDS(groups, output_file, version = 2)


# 11. Read the Figure 2 inputs -------------------------------------------

grain <- names(grains)
panels <- c("implicit", "explicit")
frameworks <- c("coverage", "sorensen")
output_file <- "process data/Figure2_quenched_expectations.rds"

groups <- readRDS("process data/Figure2_model_input.rds")
groups <- setNames(groups, vapply(groups, `[[`, character(1), "group_id"))
implicit_source <- readRDS(matched_file)
regional_sad <- implicit_source$regional_sad
sad_by_grain <- setNames(rep(list(regional_sad), length(grain)), grain)
sad_by_panel <- list(implicit = sad_by_grain, explicit = sad_by_grain)

# 12. Prepare observed LCBD ----------------------------------------------

make_observed <- function(framework) {
    implicit <- setNames(lapply(grain, function(g) {
        community <- implicit_source$matrices[[g]]
        score <- if (framework == "coverage") {
            coverage_lcbd(community, regional_sad)
        } else {
            sorensen_lcbd(community)
        }
        data.frame(richness = rowSums(community > 0), lcbd = score)
    }), grain)

    explicit <- setNames(lapply(grain, function(g) {
        data <- groups[[paste0("Simulation_grain", g)]]$observed_data
        data.frame(richness = data$richness, lcbd = data[[framework]])
    }), grain)

    list(implicit = implicit, explicit = explicit)
}

# 13. Fit quenched expectations -----------------------------------------

quenched_curve <- function(framework, sad, M, c_value) {
    if (framework == "coverage") {
        curve <- coverage_quenched(sad, M, c_value)
        return(data.frame(
            richness = curve$alpha,
            expected = curve$Quenched
        ))
    }
    curve <- sorensen_quenched(sad, M, c_value)
    data.frame(
        richness = curve$alpha,
        expected = curve$expected_Sorensen
    )
}

final_quenched_curve <- function(framework, sad, M, c_value) {
    if (framework == "coverage") {
        curve <- coverage_quenched_high_precision(sad, M, c_value)
        return(data.frame(
            richness = curve$alpha,
            expected = curve$Quenched
        ))
    }
    curve <- sorensen_quenched_log_domain(
        sad, M, c_value, alpha = 0:length(sad)
    )
    data.frame(
        richness = curve$alpha,
        expected = curve$expected_Sorensen
    )
}

# Fit the quenched aggregation parameter to site-level LCBD values.
fit_quenched_c <- function(data, sad, M, framework, fixed_c = NULL) {
    objective <- function(log_c) {
        curve <- quenched_curve(framework, sad, M, exp(log_c))
        expected <- curve$expected[match(data$richness, curve$richness)]
        if (any(!is.finite(expected))) return(1e6)
        mean((data$lcbd - expected)^2)
    }

    candidates <- if (is.null(fixed_c)) {
        optimum <- optimize(objective, log(c(1e-12, 1e6)))
        c(1e-308, 1e-12, exp(optimum$minimum), 1e6)
    } else {
        fixed_c
    }
    losses <- vapply(log(candidates), objective, numeric(1))
    c_value <- candidates[which.min(losses)]

    list(
        curve = final_quenched_curve(framework, sad, M, c_value),
        parameter = data.frame(
            model = "Model 2",
            x = NA_real_,
            c = c_value,
            y = neutral_y(M, c_value),
            mse = min(losses),
            parameter_source = if (is.null(fixed_c)) {
                "LCBD fit"
            } else {
                "known generator parameter"
            }
        )
    )
}

format_parameterized_fit <- function(fit) {
    list(
        curve = fit$curve,
        parameter = data.frame(
            model = "Fitted x",
            x = fit$x,
            c = NA_real_,
            y = NA_real_,
            mse = fit$mse,
            parameter_source = "LCBD fit"
        )
    )
}

# 14. Calculate both frameworks -----------------------------------------

# Assemble observed values and fitted curves for one LCBD framework.
calculate_framework <- function(framework) {
    observed <- make_observed(framework)
    fits <- list()
    parameters <- list()
    k <- 0L

    for (panel in panels) {
        fits[[panel]] <- list()
        for (g in grain) {
            data <- observed[[panel]][[g]]
            M <- nrow(data)
            sad <- sad_by_panel[[panel]][[g]]
            S <- length(sad)
            model2 <- fit_quenched_c(
                data,
                sad,
                M,
                framework
            )
            known_parameter <- if (panel == "implicit") {
                fit_quenched_c(
                    data,
                    sad,
                    M,
                    framework,
                    fixed_c = estimates$c[match(as.integer(g), estimates$grain)]
                )
            }
            fitted_x <- format_parameterized_fit(
                fit_parameterized_x(data, S, M, framework)
            )
            fits[[panel]][[g]] <- list(
                model2 = model2,
                known_parameter = known_parameter,
                fitted_x = fitted_x
            )

            for (fit in Filter(Negate(is.null), list(
                model2, known_parameter, fitted_x
            ))) {
                k <- k + 1L
                parameters[[k]] <- transform(
                    fit$parameter,
                    framework = framework,
                    panel = panel,
                    grain = g,
                    S = S,
                    M = M,
                    fit_objective = "site-level MSE"
                )
            }
        }
    }

    list(
        observed = observed,
        sad = sad_by_panel,
        fits = fits,
        parameters = do.call(rbind, parameters)
    )
}

result <- setNames(lapply(frameworks, calculate_framework), frameworks)

# 15. Save the common Figure 2 object -----------------------------------

saveRDS(
    list(
        frameworks = result,
        fit_objective = "site-level MSE",
        implicit_seed = implicit_source$base_seed,
        implicit_scenario = implicit_source$scenario,
        implicit_variant = "independent_nb"
    ),
    output_file,
    version = 2
)
