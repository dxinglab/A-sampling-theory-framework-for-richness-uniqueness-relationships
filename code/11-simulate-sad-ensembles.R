# Generate Fisher log-series SADs at a specified regional-abundance parameter.


# 1. Parameters -----------------------------------------------------------

N <- 211250
M <- 1250L
batch_seeds <- c(20261002L, 20251005L, 20241003L, 20231005L)
batch_sizes <- c(30L, 30L, 30L, 10L)
batch_ends <- cumsum(batch_sizes)
batch_starts <- c(1L, head(batch_ends, -1L) + 1L)
n_sad <- sum(batch_sizes)
n_x <- 30L
n_x_prime <- 30L

x_anchor <- c(0.9, 0.99, 0.999, 0.9999, 0.99999)
x_prime_anchor <- c(0.1, 0.5, 0.9, 0.99, 0.999)

output_file <- "process data/fixed_N_fisher_sads.rds"


# 2. Parameter grids ------------------------------------------------------

segmented_grid <- function(anchors, n) {
    transformed <- -log10(1 - anchors)
    intervals <- rep(
        (n - 1L) %/% (length(anchors) - 1L),
        length(anchors) - 1L
    )
    intervals[seq_len((n - 1L) %% (length(anchors) - 1L))] <-
        intervals[seq_len((n - 1L) %% (length(anchors) - 1L))] + 1L

    values <- unlist(lapply(seq_along(intervals), function(i) {
        part <- seq(
            transformed[i],
            transformed[i + 1L],
            length.out = intervals[i] + 1L
        )
        if (i > 1L) part <- part[-1L]
        part
    }), use.names = FALSE)

    grid <- 1 - 10^(-values)
    grid[c(1L, 1L + cumsum(intervals))] <- anchors
    grid
}

x_grid <- segmented_grid(x_anchor, n_x)
x_prime_grid <- segmented_grid(x_prime_anchor, n_x_prime)


# 3. Regional SADs --------------------------------------------------------

draw_sad <- function(x_index, replicate) {
    x <- x_grid[x_index]
    theta <- N * (1 - x) / x
    S <- as.integer(round(-theta * log1p(-x)))
    batch_id <- which(replicate <= batch_ends)[1L]
    replicate_in_batch <- replicate - batch_starts[batch_id] + 1L
    seed_used <- batch_seeds[batch_id] +
        (x_index - 1L) * batch_sizes[batch_id] + replicate_in_batch
    set.seed(seed_used)

    abundance <- as.numeric(sads::rls(
        n = S,
        N = N,
        alpha = theta
    ))

    list(
        abundance = abundance,
        metadata = data.frame(
            sad_id = sprintf("x%02d_sad%03d", x_index, replicate),
            batch_id = batch_id,
            x_index = x_index,
            replicate = replicate,
            replicate_in_batch = replicate_in_batch,
            seed = seed_used,
            x = x,
            theta = theta,
            S = S,
            target_N = N,
            realized_N = sum(abundance)
        )
    )
}

draw_ensemble <- function(x_index) {
    message(sprintf("Generating SADs for x grid %d / %d", x_index, n_x))
    result <- lapply(
        seq_len(n_sad),
        function(replicate) draw_sad(x_index, replicate)
    )
    list(
        abundance = lapply(result, `[[`, "abundance"),
        metadata = do.call(rbind, lapply(result, `[[`, "metadata"))
    )
}

ensembles <- lapply(seq_along(x_grid), draw_ensemble)

sad_input <- list(
    settings = list(
        N = N,
        M = M,
        n_sad = n_sad,
        batch_seeds = batch_seeds,
        batch_sizes = batch_sizes,
        x_anchor = x_anchor,
        x_prime_anchor = x_prime_anchor,
        x_grid = x_grid,
        x_prime_grid = x_prime_grid
    ),
    sad_ensembles = lapply(ensembles, `[[`, "abundance"),
    sad_metadata = do.call(rbind, lapply(ensembles, `[[`, "metadata"))
)


# 4. Save -----------------------------------------------------------------

saveRDS(sad_input, output_file, compress = "xz")
message("Saved: ", output_file)
