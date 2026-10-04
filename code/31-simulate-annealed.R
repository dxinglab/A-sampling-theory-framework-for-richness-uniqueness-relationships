# Simulate the annealed validation used by Figure 3.

# 1. Settings -------------------------------------------------------------

S <- 300L
M <- 1250L
x <- 0.9999
n_replicates <- 10000L
seed <- 20260922L
min_site_count <- 500L
min_valid_matrices <- 30L

arguments <- commandArgs(trailingOnly = TRUE)
x_prime_values <- if (length(arguments)) {
    as.numeric(arguments)
} else {
    c(0.5, 0.99)
}

source("code/01-lcbd-models.R")


# 2. Short simulation functions -----------------------------------------

number_tag <- function(value) {
    gsub("[^0-9A-Za-z]+", "p", format(value, digits = 10, scientific = FALSE))
}

draw_sad <- function(theta, expected_N) {
    as.numeric(sads::rls(n = S, N = expected_N, alpha = theta))
}

draw_community <- function(sad, c_value) {
    mu <- sad / M
    vapply(seq_len(S), function(i) {
        rnbinom(M, mu = mu[i], size = c_value * mu[i])
    }, numeric(M))
}

sum_by_richness <- function(value, richness) {
    output <- numeric(S + 1L)
    grouped <- rowsum(value, richness, reorder = FALSE)
    output[as.integer(rownames(grouped)) + 1L] <- grouped[, 1L]
    output
}

summarize_matrix <- function(community, regional_sad) {
    richness <- rowSums(community > 0)
    coverage <- coverage_lcbd(community, regional_sad)
    sorensen <- sorensen_lcbd(community)
    list(
        count = tabulate(richness + 1L, nbins = S + 1L),
        coverage_sum = sum_by_richness(coverage, richness),
        coverage_sum_sq = sum_by_richness(coverage^2, richness),
        sorensen_sum = sum_by_richness(sorensen, richness),
        sorensen_sum_sq = sum_by_richness(sorensen^2, richness)
    )
}

make_point_table <- function(
        framework, value_sum, value_sum_sq, analytic, count,
        richness_probability
) {
    sample_size <- colSums(count)
    contributors <- colSums(count > 0L)
    pooled_mean <- colSums(value_sum) / sample_size
    pooled_variance <- (
        colSums(value_sum_sq) - sample_size * pooled_mean^2
    ) / (sample_size - 1)
    pooled_sd <- sqrt(pmax(pooled_variance, 0))
    pooled_sd[sample_size < 2L] <- NA_real_
    qualified <- sample_size >= min_site_count &
        contributors >= min_valid_matrices
    if (framework == "Sørensen") qualified[1L] <- FALSE

    data.frame(
        framework = framework,
        richness = 0:S,
        sample_size = sample_size,
        contributing_replicates = contributors,
        mc_mean = pooled_mean,
        mc_sd = pooled_sd,
        analytic = analytic,
        bias = pooled_mean - analytic,
        richness_probability = richness_probability,
        empirical_probability = sample_size / (n_replicates * M),
        qualified = qualified
    )
}

metric_summary <- function(data) {
    data <- data[data$qualified & is.finite(data$bias), ]
    probability <- data$empirical_probability /
        sum(data$empirical_probability)
    data.frame(
        framework = data$framework[1L],
        probability_weighted_rmse = sqrt(sum(probability * data$bias^2)),
        probability_weighted_bias = sum(probability * data$bias)
    )
}


# 3. Simulate one x-prime setting ----------------------------------------

# Simulate and summarize all matrices for one local x-prime value.
simulate_x_prime <- function(x_prime) {
    theta <- -S / log1p(-x)
    expected_N <- theta * x / (1 - x)
    c_value <- (1 - x_prime) / x_prime
    log_y <- (c_value / M) * log1p(1 / c_value)
    y_value <- exp(log_y)

    q_bar <- 1 - log1p(-x / y_value) / log1p(-x)
    richness_probability <- dbinom(0:S, S, q_bar)
    coverage_analytic <- coverage_annealed(S, M, x, c_value, 0:S)
    sorensen_analytic <- sorensen_mean_other(
        sorensen_annealed(S, M, x, c_value, 0:S), M
    )

    set.seed(seed)
    count <- matrix(0L, n_replicates, S + 1L)
    coverage_sum <- matrix(0, n_replicates, S + 1L)
    coverage_sum_sq <- matrix(0, n_replicates, S + 1L)
    sorensen_sum <- matrix(0, n_replicates, S + 1L)
    sorensen_sum_sq <- matrix(0, n_replicates, S + 1L)
    regional_total <- numeric(n_replicates)

    for (replicate in seq_len(n_replicates)) {
        sad <- draw_sad(theta, expected_N)
        community <- draw_community(sad, c_value)
        matrix_summary <- summarize_matrix(community, sad)

        count[replicate, ] <- matrix_summary$count
        coverage_sum[replicate, ] <- matrix_summary$coverage_sum
        coverage_sum_sq[replicate, ] <- matrix_summary$coverage_sum_sq
        sorensen_sum[replicate, ] <- matrix_summary$sorensen_sum
        sorensen_sum_sq[replicate, ] <- matrix_summary$sorensen_sum_sq
        regional_total[replicate] <- sum(sad)

        if (replicate %% 100L == 0L) {
            message(
                "x' = ", x_prime, ": ", replicate,
                " / ", n_replicates
            )
        }
    }

    point_results <- rbind(
        make_point_table(
            "Coverage", coverage_sum, coverage_sum_sq, coverage_analytic,
            count, richness_probability
        ),
        make_point_table(
            "Sørensen", sorensen_sum, sorensen_sum_sq, sorensen_analytic,
            count, richness_probability
        )
    )
    metrics <- do.call(rbind, lapply(
        split(point_results, point_results$framework), metric_summary
    ))
    rownames(metrics) <- NULL

    run_id <- paste0(
        "annealed_S", S, "_M", M, "_x", number_tag(x),
        "_xp", number_tag(x_prime), "_seed", seed
    )
    result_file <- file.path(
        "process data/Figure3_annealed_validation",
        paste0(run_id, ".rds")
    )

    saveRDS(
        list(
            settings = list(
                S = S, M = M, x = x, x_prime = x_prime, c = c_value,
                n_replicates = n_replicates, seed = seed,
                local_generator = "independent negative binomial",
                metric_weight = "empirical richness frequency",
                sorensen_scale = paste(
                    "mean dissimilarity to the other M - 1 sites"
                )
            ),
            regional_total = regional_total,
            replicate_summaries = list(
                count = count,
                coverage_sum = coverage_sum,
                coverage_sum_sq = coverage_sum_sq,
                sorensen_sum = sorensen_sum,
                sorensen_sum_sq = sorensen_sum_sq
            ),
            point_results = point_results,
            metrics = metrics
        ),
        result_file,
        compress = "gzip"
    )
}


# 4. Run all settings ----------------------------------------------------

for (x_prime in x_prime_values) {
    simulate_x_prime(x_prime)
    gc()
}
