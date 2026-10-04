# Calculate the x-x-prime parameter space at a specified N parameter.


# 1. Files and shared functions ------------------------------------------

input_file <- "process data/fixed_N_fisher_sads.rds"
output_files <- c(
    coverage = "process data/coverage_fixed_N_parameter_space.rds",
    sorensen = "process data/sorensen_fixed_N_parameter_space.rds"
)

source("code/01-lcbd-models.R")
sad_input <- readRDS(input_file)

N <- sad_input$settings$N
M <- sad_input$settings$M
n_sad <- sad_input$settings$n_sad
x_grid <- sad_input$settings$x_grid
x_prime_grid <- sad_input$settings$x_prime_grid
sad_ensembles <- sad_input$sad_ensembles
sad_metadata <- sad_input$sad_metadata
if (!"batch_id" %in% names(sad_metadata)) sad_metadata$batch_id <- 1L
if (!"replicate_in_batch" %in% names(sad_metadata)) {
    sad_metadata$replicate_in_batch <- sad_metadata$replicate
}


# 2. Richness and slope functions ----------------------------------------

slope_pair <- function(alpha, S) {
    if (abs(alpha - round(alpha)) < sqrt(.Machine$double.eps)) {
        pair <- c(round(alpha) - 1L, round(alpha) + 1L)
    } else {
        pair <- c(floor(alpha), ceiling(alpha))
    }
    pair <- pmax(1L, pmin(S, pair))
    if (pair[1L] == pair[2L]) {
        pair <- if (pair[1L] == 1L) c(1L, 2L) else c(S - 1L, S)
    }
    as.integer(pair)
}

point_slope <- function(alpha, expected) {
    unname(diff(expected) / diff(log(alpha)))
}

log_richness_derivative <- function(FUN, alpha, h = 0.01) {
    lower <- alpha * exp(-h)
    upper <- alpha * exp(h)
    (FUN(upper) - FUN(lower)) / (2 * h)
}

conditional_mean_from_q <- function(q) {
    sum(q) / (-expm1(sum(log1p(-q))))
}

annealed_typical_richness <- function(S, x, c) {
    y <- neutral_y(M, c)
    q_bar <- 1 - log1p(-x / y) / log1p(-x)
    S * q_bar / (-expm1(S * log1p(-q_bar)))
}


# 3. Targeted quenched expectations --------------------------------------

gauss_legendre <- function(n = 64L) {
    i <- seq_len(n - 1L)
    beta <- i / sqrt(4 * i^2 - 1)
    matrix_J <- matrix(0, n, n)
    matrix_J[cbind(i, i + 1L)] <- beta
    matrix_J[cbind(i + 1L, i)] <- beta
    decomposition <- eigen(matrix_J, symmetric = TRUE)
    order <- order(decomposition$values)
    list(
        nodes = (decomposition$values[order] + 1) / 2,
        weights = decomposition$vectors[1L, order]^2
    )
}

# Calculate species inclusion probabilities conditional on richness.
conditional_inclusion <- function(q, alpha) {
    S <- length(q)
    K <- as.integer(alpha)
    log_odds <- qlogis(pmin(1 - 1e-14, pmax(1e-14, q)))
    eta <- uniroot(
        function(value) sum(plogis(log_odds + value)) - K,
        c(-60, 60),
        tol = 1e-11
    )$root
    tilted <- plogis(log_odds + eta)

    forward <- matrix(0, S + 1L, K + 1L)
    forward[1L, 1L] <- 1
    for (i in seq_len(S)) {
        previous <- forward[i, ]
        shifted <- c(0, previous[-(K + 1L)])
        forward[i + 1L, ] <-
            previous * (1 - tilted[i]) + shifted * tilted[i]
    }

    total <- forward[S + 1L, K + 1L]
    adjoint <- numeric(K + 1L)
    adjoint[K + 1L] <- 1
    derivative <- numeric(S)
    for (i in S:1L) {
        previous <- forward[i, ]
        shifted <- c(0, previous[-(K + 1L)])
        derivative[i] <- sum(adjoint * (shifted - previous))
        adjoint <- adjoint * (1 - tilted[i]) +
            c(adjoint[-1L] * tilted[i], 0)
    }

    pmin(
        1,
        pmax(0, tilted + tilted * (1 - tilted) * derivative / total)
    )
}

sorensen_quenched_at <- function(q, inclusion, alpha, quadrature) {
    log_factor <- outer(
        q,
        1 - quadrature$nodes,
        function(q_i, one_minus_t) log1p(-q_i * one_minus_t)
    )
    base <- alpha * log(quadrature$nodes) +
        colSums(log_factor) + log(quadrature$weights)
    integral <- rowSums(exp(sweep(-log_factor, 2L, base, "+")))
    1 - 2 * sum(q * inclusion * integral)
}

quenched_pair <- function(sad, c_value, pair, quadrature) {
    q <- neutral_occurrence_probability(sad, M, c_value)
    regional_weight <- sad / sum(sad)

    values <- lapply(pair, function(alpha) {
        inclusion <- conditional_inclusion(q, alpha)
        c(
            coverage = M / (M - 1) *
                (1 - sum(regional_weight * inclusion)),
            sorensen = sorensen_quenched_at(
                q, inclusion, alpha, quadrature
            )
        )
    })
    do.call(rbind, values)
}


# 4. One x-x-prime combination -------------------------------------------

# Calculate all three model slopes for one x and x-prime combination.
calculate_point <- function(grid_row, quadrature) {
    x_index <- grid_row$x_index
    x_prime_index <- grid_row$x_prime_index
    x <- x_grid[x_index]
    x_prime <- x_prime_grid[x_prime_index]
    c_value <- (1 - x_prime) / x_prime
    sads <- sad_ensembles[[x_index]]
    metadata <- sad_metadata[sad_metadata$x_index == x_index, ]
    stopifnot(nrow(metadata) == length(sads))
    S <- length(sads[[1L]])

    alpha3 <- annealed_typical_richness(S, x, c_value)
    parameter_slopes <- vapply(
        c("coverage", "sorensen"),
        function(framework) log_richness_derivative(
            function(alpha) parameterized_expectation(
                framework, S, M, x, alpha
            ),
            alpha3
        ),
        numeric(1)
    )
    annealed_slopes <- vapply(
        c("coverage", "sorensen"),
        function(framework) log_richness_derivative(
            function(alpha) annealed_expectation(
                framework, S, M, x, c_value, alpha
            ),
            alpha3
        ),
        numeric(1)
    )

    quenched <- lapply(sads, function(sad) {
        q <- neutral_occurrence_probability(sad, M, c_value)
        alpha2 <- conditional_mean_from_q(q)
        pair <- slope_pair(alpha2, length(sad))
        values <- quenched_pair(sad, c_value, pair, quadrature)
        data.frame(
            alpha2 = alpha2,
            lower = pair[1L],
            upper = pair[2L],
            slope_model2_coverage = point_slope(
                pair, values[, "coverage"]
            ),
            slope_model2_sorensen = point_slope(
                pair, values[, "sorensen"]
            )
        )
    })
    quenched <- do.call(rbind, quenched)

    unique_pairs <- unique(quenched[, c("lower", "upper")])
    pair_slopes <- lapply(seq_len(nrow(unique_pairs)), function(i) {
        pair <- unlist(unique_pairs[i, ], use.names = FALSE)
        data.frame(
            lower = pair[1L],
            upper = pair[2L],
            slope_model3_coverage = point_slope(
                pair,
                annealed_expectation(
                    "coverage", S, M, x, c_value, pair
                )
            ),
            slope_model3_sorensen = point_slope(
                pair,
                annealed_expectation(
                    "sorensen", S, M, x, c_value, pair
                )
            )
        )
    })
    pair_slopes <- do.call(rbind, pair_slopes)
    quenched <- merge(
        quenched,
        pair_slopes,
        by = c("lower", "upper"),
        sort = FALSE
    )

    base <- data.frame(
        x_index = x_index,
        x_prime_index = x_prime_index,
        batch_id = metadata$batch_id,
        replicate = metadata$replicate,
        replicate_in_batch = metadata$replicate_in_batch,
        sad_id = metadata$sad_id,
        seed = metadata$seed,
        x = x,
        x_prime = x_prime,
        c = c_value,
        N = N,
        S = S,
        alpha_target = quenched$alpha2,
        alpha_lower = quenched$lower,
        alpha_upper = quenched$upper,
        alpha3_target = alpha3
    )

    replicate_rows <- lapply(c("coverage", "sorensen"), function(framework) {
        model2 <- quenched[[paste0("slope_model2_", framework)]]
        model3 <- quenched[[paste0("slope_model3_", framework)]]
        data.frame(
            base,
            framework = framework,
            point_type = "model2_typical_alpha",
            slope_model2 = model2,
            slope_model3 = model3,
            slope_model2_at_alpha2 = model2,
            slope_model3_pair_at_alpha2 = model3,
            calculation_status = ifelse(
                is.finite(model2) & is.finite(model3),
                "ok",
                "calculation_failure"
            )
        )
    })
    replicate_rows <- do.call(rbind, replicate_rows)

    parameter_rows <- lapply(c("coverage", "sorensen"), function(framework) {
        rows <- replicate_rows[replicate_rows$framework == framework, ]
        data.frame(
            x_index = x_index,
            x_prime_index = x_prime_index,
            x = x,
            x_prime = x_prime,
            c = c_value,
            N = N,
            S = S,
            framework = framework,
            point_type = "model3_typical_alpha",
            n_sad = n_sad,
            n_success = sum(rows$calculation_status == "ok"),
            alpha_target_mean = alpha3,
            slope_model1_mean = parameter_slopes[framework],
            slope_model3_mean = annealed_slopes[framework],
            alpha2_mean = mean(rows$alpha_target),
            alpha2_sd = stats::sd(rows$alpha_target),
            slope_model2_at_alpha2_mean =
                mean(rows$slope_model2_at_alpha2),
            slope_model2_at_alpha2_sd =
                stats::sd(rows$slope_model2_at_alpha2),
            slope_model3_pair_at_alpha2_mean =
                mean(rows$slope_model3_pair_at_alpha2),
            slope_model3_pair_at_alpha2_sd =
                stats::sd(rows$slope_model3_pair_at_alpha2)
        )
    })

    list(
        parameter = do.call(rbind, parameter_rows),
        replicate = replicate_rows
    )
}


# 5. Parameter space ------------------------------------------------------

grid <- expand.grid(
    x_prime_index = seq_along(x_prime_grid),
    x_index = seq_along(x_grid)
)
grid$x <- x_grid[grid$x_index]
grid$x_prime <- x_prime_grid[grid$x_prime_index]
grid <- grid[grid$x_prime < grid$x, ]

quadrature <- gauss_legendre(64L)
result <- lapply(seq_len(nrow(grid)), function(i) {
    message(sprintf("Parameter combination %d / %d", i, nrow(grid)))
    calculate_point(grid[i, ], quadrature)
})

parameter_rows <- do.call(rbind, lapply(result, `[[`, "parameter"))
replicate_rows <- do.call(rbind, lapply(result, `[[`, "replicate"))


# 6. Save one result for each framework ----------------------------------

# Save one compact parameter-space result for each LCBD framework.
save_framework_result <- function(framework) {
    output <- list(
        settings = c(
            sad_input$settings,
            list(
                framework = framework,
                N_parameter = N,
                parameter_space = "x_prime < x",
                quenched_slope = paste0(
                    "Two-point slope between the integer richness values ",
                    "bracketing the SAD-specific typical alpha."
                ),
                analytic_slope = paste0(
                    "Derivative with respect to log richness at the exact ",
                    "annealed typical alpha, using h = 0.01."
                ),
                matched_comparison = paste0(
                    "Annealed and quenched slopes use the same two integer ",
                    "richness values around each SAD-specific typical alpha."
                )
            )
        ),
        sad_metadata = sad_input$sad_metadata,
        slope_summary = parameter_rows[
            parameter_rows$framework == framework,
        ],
        replicate_slopes = replicate_rows[
            replicate_rows$framework == framework,
        ]
    )
    saveRDS(output, output_files[framework], compress = "xz")
    message("Saved: ", output_files[framework])
}

save_framework_result("coverage")
save_framework_result("sorensen")
