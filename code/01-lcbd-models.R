# ==============================================================================
# LCBD observations and neutral expectations
#
# Copyright (c) 2026 Xiaoning Wang and Dingliang Xing.
# Licensed under the MIT License. See ../LICENSE.
#
# 1. Parameterized:
#      B_CD^par(alpha)
#
# 2. Quenched conditional:
#      B_CD^Q(alpha | n_1,...,n_S)
#      (Fixed realized regional SAD)
#
# 3. Fisher-averaged conditional (Annealed Exact):
#      B_CD^A(alpha)
#      (Fisher log-series averaged after conditioning on R = alpha)
#
# 4. Fisher-averaged Large-S Approximation:
#      Approximate annealed curve for analytical tractability.
#
# LCBD convention used here:
#      LCBD_CD = M/(M-1) * (1 - C)
#
# where C = sum_i p_i I_i and p_i = n_i / sum_j n_j
# ==============================================================================


# ==============================================================================
# 1. Strict-Neutral Occurrence Model
# ==============================================================================

neutral_log_y <- function(M, c) {
    if (M <= 1) stop("M must be > 1.")
    if (c <= 0) stop("c must be > 0.")
    (c / M) * log1p(1 / c)
}

neutral_y <- function(M, c) {
    exp(neutral_log_y(M, c))
}

neutral_occurrence_probability <- function(n, M, c) {
    # q_n = 1 - y^(-n)
    -expm1(-n * neutral_log_y(M, c))
}


# ==============================================================================
# 2. Parameterized Richness-Matched Coverage Baseline
# ==============================================================================

coverage_parameterized <- function(S, M, x, alpha) {
    if (S <= 0) stop("S must be positive.")
    if (M <= 1) stop("M must be > 1.")
    if (x <= 0 || x >= 1) stop("x must be between 0 and 1.")
    
    out <- numeric(length(alpha))
    scale_factor <- M / (M - 1)
    
    # Theoretical boundaries
    out[alpha <= 0] <- scale_factor
    out[alpha >= S] <- 0
    
    inside <- alpha > 0 & alpha < S
    
    if (any(inside)) {
        a <- alpha[inside]
        log_1mx <- log1p(-x)
        
        # 1 - (1-x)^(alpha/S)
        one_minus_term <- -expm1((a / S) * log_1mx)
        
        out[inside] <- scale_factor * (1 - one_minus_term / x)
    }
    
    out
}


# ==============================================================================
# 3. Quenched Conditional Coverage Baseline
#
# Richness probabilities are updated in log space to avoid underflow. The
# recursion updates E[C | R = k] directly.
#
# Transition equation:
# E_new[k] = E_old[k] * (1 - rho_k) + (E_old[k-1] + w_i) * rho_k
# where rho_k = P(I_i = 1 | R_new = k)
# ==============================================================================

coverage_quenched <- function(n, M, c) {
    n <- as.numeric(n)
    
    if (length(n) == 0) stop("n must contain at least one species.")
    if (any(!is.finite(n)) || any(n <= 0)) stop("All regional abundances must be positive.")
    if (M <= 1) stop("M must be > 1.")
    
    S <- length(n)
    weight <- n / sum(n)
    q <- neutral_occurrence_probability(n = n, M = M, c = c)
    
    # State representation in log-probability space to prevent underflow.
    # log_rp[k] represents log(P(R = k-1))
    log_rp <- rep(-Inf, S + 1)
    log_rp[1] <- 0
    
    # E[C | R=k-1]
    expected_C <- rep(0, S + 1)
    
    for (i in seq_len(S)) {
        qi <- q[i]
        wi <- weight[i]
        
        log_1mq <- log1p(-qi)
        log_q <- log(qi)
        
        # Log-probabilities of transitioning to state k from k and from k-1
        term_0 <- log_rp + log_1mq
        term_1 <- c(-Inf, log_rp[-length(log_rp)]) + log_q
        
        # LogSumExp to find the updated log_rp
        max_term <- pmax(term_0, term_1)
        valid <- is.finite(max_term)
        
        log_rp_new <- rep(-Inf, S + 1)
        log_rp_new[valid] <- max_term[valid] + log(
            exp(term_0[valid] - max_term[valid]) +
                exp(term_1[valid] - max_term[valid])
        )
        
        # rho_k: Conditional probability that species i is present given new richness is k
        log_rho <- term_1 - log_rp_new
        rho <- exp(log_rho)
        rho[is.na(rho)] <- 0 # Handles -Inf - (-Inf) when states are mathematically unreachable
        
        # Direct DP update of the conditional expectations
        expected_C_prev <- c(0, expected_C[-length(expected_C)])
        expected_C <- expected_C * (1 - rho) + (expected_C_prev + wi) * rho
        
        log_rp <- log_rp_new
    }
    
    # Apply exact boundary values.
    expected_C[1] <- 0
    expected_C[S + 1] <- 1
    
    scale_factor <- M / (M - 1)
    expected_LCBD <- scale_factor * (1 - expected_C)
    
    data.frame(
        alpha = 0:S,
        richness_probability = exp(log_rp),
        expected_C = expected_C,
        Quenched = expected_LCBD
    )
}


# ==============================================================================
# 4. Fisher-Averaged Probability Generating Functions
# ==============================================================================

U_presence <- function(t, x, y) {
    numerator <- log1p(-x * t) - log1p(-(x / y) * t)
    denominator <- log1p(-x) - log1p(-x / y)
    numerator / denominator
}

V_absence <- function(t, x, y) {
    r <- x / y
    log1p(-r * t) / log1p(-r)
}


# ==============================================================================
# 5. Exact Fisher-Averaged Conditional Coverage Baseline (Annealed)
#
# Integration by parts evaluates 1 - E[C] directly and avoids endpoint
# singularities and loss of precision near E[C] = 1.
# ==============================================================================

coverage_annealed_one <- function(
        S, M, x, c, alpha,
        rel.tol = 1e-9, abs.tol = 1e-12
) {
    scale_factor <- M / (M - 1)
    
    # Theoretical boundaries
    if (alpha <= 0) return(scale_factor)
    if (alpha >= S) return(0)
    
    y <- neutral_y(M = M, c = c)
    
    # Integrand derived via integration by parts over z, where t = 1 - exp(-z)
    integrand_z <- function(z) {
        t <- -expm1(-z)
        
        U <- U_presence(t = t, x = x, y = y)
        V <- V_absence(t = t, x = x, y = y)
        
        # dV(t) / dt
        r <- x / y
        dV <- (-r / (1 - r * t)) / log1p(-r)
        
        U <- pmax(U, .Machine$double.xmin)
        V <- pmax(V, .Machine$double.xmin)
        
        # Log-scale summation avoids underflow when computing large powers of U and V
        log_val <- log(S - alpha) + 
            alpha * log(U) + 
            (S - alpha - 1) * log(V) + 
            log(dV) - z
        
        exp(log_val)
    }
    
    # Truncate the transformed integral at z = 40.
    one_minus_expected_C <- integrate(
        integrand_z,
        lower = 0,
        upper = 40, 
        subdivisions = 1000,
        rel.tol = rel.tol,
        abs.tol = abs.tol
    )$value
    
    scale_factor * one_minus_expected_C
}


# ==============================================================================
# 6. Exact Annealed Curve Vectorization
# ==============================================================================

coverage_annealed <- function(S, M, x, c, alpha) {
    vapply(
        alpha,
        FUN = function(a) {
            coverage_annealed_one(S = S, M = M, x = x, c = c, alpha = a)
        },
        FUN.VALUE = numeric(1)
    )
}


# ==============================================================================
# 7. Large-S Approximation for the Annealed Curve
# ==============================================================================

fisher_mean <- function(x) {
    -x / ((1 - x) * log1p(-x))
}

mu_plus_fisher <- function(x, y) {
    mu_reg <- fisher_mean(x)
    log_den <- log1p(-x)
    
    EN_y_minus_N <- -(x / y) / ((1 - x / y) * log_den)
    qbar <- 1 - log1p(-x / y) / log_den
    
    EN_q <- mu_reg - EN_y_minus_N
    EN_q / qbar
}

mu_minus_fisher <- function(x, y) {
    absence_probability <- log1p(-x / y) / log1p(-x)
    
    EN_absent <- -(x / y) / ((1 - x / y) * log1p(-x))
    EN_absent / absence_probability
}

coverage_annealed_largeS <- function(S, M, x, c, alpha) {
    if (M <= 1) stop("M must be > 1.")
    
    out <- numeric(length(alpha))
    scale_factor <- M / (M - 1)
    
    # Theoretical boundaries
    out[alpha <= 0] <- scale_factor
    out[alpha >= S] <- 0
    
    inside <- alpha > 0 & alpha < S
    
    if (any(inside)) {
        a <- alpha[inside]
        y <- neutral_y(M = M, c = c)
        
        mu_plus <- mu_plus_fisher(x = x, y = y)
        mu_minus <- mu_minus_fisher(x = x, y = y)
        
        rho <- a / S
        
        expected_C <- rho * mu_plus / (rho * mu_plus + (1 - rho) * mu_minus)
        
        out[inside] <- scale_factor * (1 - expected_C)
    }
    
    out
}


# ==============================================================================
# 8. Sørensen expectations on the manuscript response scale
# ==============================================================================

# The original Sørensen equations return M / (M - 1) times the mean
# dissimilarity to the other sites. These functions make that conversion once.
sorensen_scale <- function(M) {
    M / (M - 1)
}

sorensen_mean_other <- function(value, M) {
    value / sorensen_scale(M)
}

sorensen_parameterized_one <- function(S, x, alpha) {
    if (alpha <= 0) return(1)

    v <- -log1p(-x)
    nu <- exp(-v)
    growth_minus_one <- expm1(v * alpha / S)
    growth <- 1 + growth_minus_one
    ratio_minus_one <- growth_minus_one *
        (2 - nu * (growth + 1)) / (-expm1(-v))
    delta <- -alpha + (S / v) * log1p(ratio_minus_one)
    delta / alpha
}

sorensen_parameterized <- function(S, M, x, alpha) {
    sorensen_scale(M) * vapply(
        alpha,
        sorensen_parameterized_one,
        numeric(1),
        S = S,
        x = x
    )
}

poisson_binomial_coefficients <- function(q) {
    S <- length(q)
    probability <- numeric(S + 1L)
    probability[1L] <- 1

    for (q_i in q) {
        probability <- probability * (1 - q_i) +
            c(0, probability[-(S + 1L)]) * q_i
    }
    probability
}

sorensen_quenched <- function(n, M, c) {
    n <- as.numeric(n)
    S <- length(n)
    q <- neutral_occurrence_probability(n, M, c)
    richness_probability <- poisson_binomial_coefficients(q)

    leave_one_out <- matrix(0, nrow = S, ncol = S)
    for (i in seq_len(S)) {
        probability <- numeric(S)
        probability[1L] <- 1
        for (q_j in q[-i]) {
            probability <- probability * (1 - q_j) +
                c(0, probability[-S]) * q_j
        }
        leave_one_out[i, ] <- probability
    }

    expected <- numeric(S + 1L)
    expected[1L] <- 1
    for (alpha in seq_len(S)) {
        denominator <- richness_probability[alpha + 1L]
        if (!is.finite(denominator) ||
                denominator <= .Machine$double.xmin) {
            expected[alpha + 1L] <- NA_real_
            next
        }
        integral <- as.vector(
            leave_one_out %*% (1 / (alpha + seq_len(S)))
        )
        numerator <- sum(q^2 * leave_one_out[, alpha] * integral)
        expected[alpha + 1L] <- 1 - 2 * numerator / denominator
    }

    data.frame(
        alpha = 0:S,
        richness_probability = richness_probability,
        expected_Sorensen = expected,
        Quenched = sorensen_scale(M) * expected
    )
}

fisher_occurrence_moments <- function(x, y) {
    R1 <- log1p(-x / y) / log1p(-x)
    R2 <- log1p(-x / y^2) / log1p(-x)
    pbar <- 1 - R1
    q2 <- 1 - 2 * R1 + R2

    list(
        pbar = pbar,
        q2 = q2,
        q_plus = q2 / pbar,
        q_minus = (pbar - q2) / (1 - pbar),
        R1 = R1,
        R2 = R2
    )
}

sorensen_annealed_one <- function(
        S, M, x, c, alpha,
        rel.tol = 1e-9, abs.tol = 1e-12
) {
    scale_factor <- sorensen_scale(M)
    if (alpha <= 0) return(scale_factor)

    moments <- fisher_occurrence_moments(x, neutral_y(M, c))
    integrand <- function(t) {
        exp_minus_t <- exp(-t)
        z <- pmax(1 - exp_minus_t, .Machine$double.xmin)
        h_plus <- pmax(
            1 - moments$q_plus * exp_minus_t,
            .Machine$double.xmin
        )
        h_minus <- pmax(
            1 - moments$q_minus * exp_minus_t,
            .Machine$double.xmin
        )
        exp(
            alpha * log(z) +
                (alpha - 1) * log(h_plus) +
                (S - alpha) * log(h_minus) - t
        )
    }

    integral <- integrate(
        integrand, 0, 40,
        rel.tol = rel.tol,
        abs.tol = abs.tol
    )$value
    scale_factor * (1 - 2 * alpha * moments$q_plus * integral)
}

sorensen_annealed <- function(S, M, x, c, alpha) {
    vapply(
        alpha,
        sorensen_annealed_one,
        numeric(1),
        S = S,
        M = M,
        x = x,
        c = c
    )
}

sorensen_annealed_largeS <- function(S, M, x, c, alpha) {
    out <- numeric(length(alpha))
    out[alpha <= 0] <- sorensen_scale(M)
    inside <- alpha > 0

    if (any(inside)) {
        moments <- fisher_occurrence_moments(x, neutral_y(M, c))
        rho <- alpha[inside] / S
        dissimilarity <- 1 -
            2 * rho * moments$q_plus /
            (rho * (1 + moments$q_plus) +
                (1 - rho) * moments$q_minus)
        out[inside] <- sorensen_scale(M) * dissimilarity
    }
    out
}


# ==============================================================================
# 9. Log-domain quenched Sørensen evaluation
# ==============================================================================

logspace_add <- function(x, y) {
    maximum <- pmax(x, y)
    out <- maximum + log(exp(x - maximum) + exp(y - maximum))
    out[!is.finite(maximum)] <- -Inf
    out
}

logspace_sum <- function(x) {
    maximum <- max(x)
    if (!is.finite(maximum)) return(-Inf)
    maximum + log(sum(exp(x - maximum)))
}

sorensen_quenched_log_domain <- function(n, M, c, alpha) {
    n <- as.numeric(n)
    S <- length(n)
    log_absence <- -n * neutral_log_y(M, c)
    presence <- -expm1(log_absence)
    log_presence <- log(presence)

    log_probability <- rep(-Inf, S + 1L)
    log_probability[1L] <- 0
    for (i in seq_len(S)) {
        previous <- c(-Inf, log_probability[-length(log_probability)])
        log_probability <- logspace_add(
            log_probability + log_absence[i],
            previous + log_presence[i]
        )
    }

    requested <- sort(unique(as.integer(alpha)))
    terms <- matrix(-Inf, nrow = S, ncol = length(requested))
    for (i in seq_len(S)) {
        log_loo <- rep(-Inf, S)
        log_loo[1L] <- 0
        for (j in setdiff(seq_len(S), i)) {
            previous <- c(-Inf, log_loo[-S])
            log_loo <- logspace_add(
                log_loo + log_absence[j],
                previous + log_presence[j]
            )
        }
        for (j in seq_along(requested)) {
            a <- requested[j]
            if (a == 0L) next
            log_integral <- logspace_sum(
                log_loo - log(a + seq_len(S))
            )
            terms[i, j] <- 2 * log_presence[i] +
                log_loo[a] + log_integral
        }
    }

    expected <- numeric(length(requested))
    expected[requested == 0L] <- 1
    for (j in which(requested > 0L)) {
        a <- requested[j]
        expected[j] <- 1 - 2 * exp(
            logspace_sum(terms[, j]) - log_probability[a + 1L]
        )
    }

    list(
        alpha = requested,
        log_richness_probability = log_probability,
        expected_Sorensen = expected,
        Quenched = sorensen_scale(M) * expected
    )
}


# ==============================================================================
# 10. High-precision quenched Coverage evaluation
# ==============================================================================

coverage_quenched_high_precision <- function(n, M, c, precBits = 512L) {
    S <- length(n)
    n <- Rmpfr::mpfr(n, precBits)
    M <- Rmpfr::mpfr(M, precBits)
    c <- Rmpfr::mpfr(c, precBits)
    one <- Rmpfr::mpfr(1, precBits)
    q <- -expm1(-n * (c / M) * log1p(one / c))
    weight <- n / sum(n)
    probability <- coverage <- Rmpfr::mpfr(rep(0, S + 1L), precBits)
    probability[1L] <- one

    for (i in seq_len(S)) {
        previous_probability <- c(
            one - one,
            probability[-length(probability)]
        )
        previous_coverage <- c(
            one - one,
            coverage[-length(coverage)]
        )
        coverage <- coverage * (one - q[i]) +
            (previous_coverage + weight[i] * previous_probability) * q[i]
        probability <- probability * (one - q[i]) +
            previous_probability * q[i]
    }

    expected_C <- coverage / probability
    expected_C[c(1L, S + 1L)] <- c(0, 1)
    data.frame(
        alpha = 0:S,
        richness_probability = Rmpfr::asNumeric(probability),
        expected_C = Rmpfr::asNumeric(expected_C),
        Quenched = Rmpfr::asNumeric(
            (M / (M - one)) * (one - expected_C)
        )
    )
}


# ==============================================================================
# 11. Observed LCBD scores
# ==============================================================================

coverage_lcbd <- function(community, regional_sad = colSums(community)) {
    incidence <- community > 0
    M <- nrow(community)
    weight <- regional_sad / sum(regional_sad)
    M / (M - 1) * (1 - drop(incidence %*% weight))
}

sorensen_lcbd <- function(community) {
    incidence <- community > 0
    M <- nrow(incidence)
    richness <- rowSums(incidence)
    richness_levels <- sort(unique(richness))
    group <- match(richness, richness_levels)
    group_size <- tabulate(group)
    occurrence <- vapply(
        seq_along(richness_levels),
        function(i) colSums(incidence[group == i, , drop = FALSE]),
        numeric(ncol(incidence))
    )
    shared <- incidence %*% occurrence
    similarity <- 2 * shared / outer(richness, richness_levels, "+")
    if (richness_levels[1L] == 0L) {
        similarity[richness == 0L, 1L] <- group_size[1L]
    }
    similarity[!is.finite(similarity)] <- 0

    # The self-similarity of one is included in the sum. Dividing the
    # resulting dissimilarity sum by M - 1 gives the mean to other sites.
    (M - rowSums(similarity)) / (M - 1)
}


# ==============================================================================
# 12. Common model interfaces and fitting functions
# ==============================================================================

parameterized_expectation <- function(framework, S, M, x, alpha) {
    if (framework == "coverage") {
        return(coverage_parameterized(S, M, x, alpha))
    }
    sorensen_mean_other(
        sorensen_parameterized(S, M, x, alpha),
        M
    )
}

annealed_expectation <- function(framework, S, M, x, c, alpha) {
    if (framework == "coverage") {
        return(coverage_annealed(S, M, x, c, alpha))
    }
    sorensen_mean_other(
        sorensen_annealed(S, M, x, c, alpha),
        M
    )
}

annealed_expectation_complete <- function(
        framework, S, M, x, c, alpha
) {
    approximate <- if (framework == "coverage") {
        coverage_annealed_largeS(S, M, x, c, alpha)
    } else {
        sorensen_mean_other(
            sorensen_annealed_largeS(S, M, x, c, alpha),
            M
        )
    }
    exact <- vapply(alpha, function(a) {
        tryCatch(
            annealed_expectation(framework, S, M, x, c, a),
            error = function(error) NA_real_
        )
    }, numeric(1))
    exact[!is.finite(exact)] <- approximate[!is.finite(exact)]
    exact
}

fit_parameterized_x <- function(
        data, S, M, framework,
        x_bounds = c(0.5, 1 - 1e-10)
) {
    bounds <- -log1p(-x_bounds)
    objective <- function(v) {
        expected <- parameterized_expectation(
            framework, S, M, -expm1(-v), data$richness
        )
        mean((data$lcbd - expected)^2)
    }
    optimum <- optimize(objective, bounds)
    candidates <- c(bounds, optimum$minimum)
    v <- candidates[which.min(vapply(candidates, objective, numeric(1)))]
    richness <- 0:S

    list(
        curve = data.frame(
            richness = richness,
            expected = parameterized_expectation(
                framework, S, M, -expm1(-v), richness
            )
        ),
        x = -expm1(-v),
        mse = objective(v)
    )
}

fit_annealed_xc <- function(
        data, S, M, framework,
        x_bounds = c(0.5, 1 - 1e-10),
        c_bounds = c(1e-12, 1e6),
        x_starts = c(0.9, 0.99, 0.999, 0.9999, 0.99999),
        c_starts = c(0.01, 0.1, 1, 10, 1e3, 1e5),
        curve_richness = 0:S
) {
    richness <- sort(unique(data$richness))
    group <- match(data$richness, richness)
    count <- tabulate(group, nbins = length(richness))
    sum_y <- drop(rowsum(data$lcbd, group, reorder = TRUE))
    sum_y2 <- drop(rowsum(data$lcbd^2, group, reorder = TRUE))
    v_bounds <- -log1p(-x_bounds)
    lower <- c(v_bounds[1L], log(c_bounds[1L]))
    upper <- c(v_bounds[2L], log(c_bounds[2L]))

    objective <- function(theta) {
        expected <- tryCatch(
            suppressWarnings(annealed_expectation(
                framework, S, M,
                -expm1(-theta[1L]), exp(theta[2L]), richness
            )),
            error = function(error) rep(NA_real_, length(richness))
        )
        if (any(!is.finite(expected))) return(1e6)
        sum(sum_y2 - 2 * expected * sum_y + count * expected^2) /
            sum(count)
    }

    starts <- expand.grid(x = x_starts, c = c_starts)
    starts$loss <- mapply(
        function(x, c) objective(c(-log1p(-x), log(c))),
        starts$x,
        starts$c
    )
    starts <- starts[order(starts$loss)[seq_len(3L)], ]
    fits <- lapply(seq_len(nrow(starts)), function(i) {
        optim(
            c(-log1p(-starts$x[i]), log(starts$c[i])),
            objective,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            control = list(maxit = 250L)
        )
    })
    fit <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
    x <- -expm1(-fit$par[1L])
    c_value <- exp(fit$par[2L])
    list(
        curve = data.frame(
            richness = curve_richness,
            expected = annealed_expectation(
                framework, S, M, x, c_value, curve_richness
            )
        ),
        x = x,
        c = c_value,
        y = neutral_y(M, c_value),
        mse = fit$value,
        convergence = fit$convergence,
        boundary = any(abs(fit$par - lower) < 1e-6) ||
            any(abs(fit$par - upper) < 1e-6)
    )
}

# Read the first worksheet from a plain .xlsx analysis file.
read_xlsx_data <- function(path) {
    strings_document <- xml2::read_xml(unz(path, "xl/sharedStrings.xml"))
    strings <- xml2::xml_text(xml2::xml_find_all(
        strings_document, ".//*[local-name()='si']"
    ))
    document <- xml2::read_xml(unz(path, "xl/worksheets/sheet1.xml"))
    cells <- xml2::xml_find_all(
        document,
        ".//*[local-name()='sheetData']/*[local-name()='row']/*[local-name()='c']"
    )
    reference <- xml2::xml_attr(cells, "r")
    row <- as.integer(gsub("[A-Z]+", "", reference))
    column <- vapply(
        strsplit(gsub("[0-9]+", "", reference), ""),
        function(x) Reduce(
            function(value, letter) value * 26L + match(letter, LETTERS),
            x, init = 0L
        ),
        integer(1)
    )
    value <- xml2::xml_text(xml2::xml_find_first(
        cells, "./*[local-name()='v']"
    ))
    type <- xml2::xml_attr(cells, "t")
    shared <- !is.na(type) & type == "s" & !is.na(value)
    value[shared] <- strings[as.integer(value[shared]) + 1L]
    table <- matrix(NA_character_, max(row), max(column))
    table[cbind(row, column)] <- value
    result <- as.data.frame(table[-1L, , drop = FALSE], stringsAsFactors = FALSE)
    names(result) <- table[1L, ]
    type.convert(result, as.is = TRUE)
}
