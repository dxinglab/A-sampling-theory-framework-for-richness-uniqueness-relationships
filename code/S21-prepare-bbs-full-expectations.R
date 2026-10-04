# Prepare complete BCR-specific expectations for Figures S2 and S3.
# Run after 51-prepare-bbs-data.R.

# 1. Paths and settings ---------------------------------------------------

input_file <- "process data/Figure5_BBS_annealed_expectations.rds"
output_file <- "process data/FiguresS2_S3_BBS_full_expectations.rds"
x_bounds <- c(0.5, 1 - 1e-10)

# 2. Read fitted annealed parameters -------------------------------------

source("code/01-lcbd-models.R")
result <- readRDS(input_file)
site <- result$site_data
annealed_parameter <- result$parameters

model_functions <- list(
    Coverage = list(
        parameterized = coverage_parameterized,
        annealed = coverage_annealed,
        approximation = coverage_annealed_largeS
    ),
    Sorensen = list(
        parameterized = function(S, M, x, alpha) {
            parameterized_expectation("sorensen", S, M, x, alpha)
        },
        annealed = function(S, M, x, c, alpha) {
            annealed_expectation("sorensen", S, M, x, c, alpha)
        },
        approximation = function(S, M, x, c, alpha) {
            sorensen_mean_other(
                sorensen_annealed_largeS(S, M, x, c, alpha),
                M
            )
        }
    )
)

# 3. Complete the fitted expectations within each BCR --------------------

complete_annealed_curve <- function(
        S, M, x, c, alpha, model_function, approximation
) {
    approximate <- approximation(S, M, x, c, alpha)
    exact <- vapply(alpha, function(a) {
        tryCatch(
            model_function(S, M, x, c, a),
            error = function(e) NA_real_
        )
    }, numeric(1))
    fallback <- !is.finite(exact)
    exact[fallback] <- approximate[fallback]
    list(expected = exact, numerical_fallback = fallback)
}

parameterized_rows <- vector("list", nrow(annealed_parameter))
curve_rows <- vector("list", 2L * nrow(annealed_parameter))

for (i in seq_len(nrow(annealed_parameter))) {
    parameter <- annealed_parameter[i, ]
    framework <- parameter$framework
    functions <- model_functions[[framework]]
    response <- if (framework == "Coverage") {
        site$coverage_lcbd
    } else {
        site$sorensen_lcbd
    }
    data <- data.frame(
        richness = site$richness,
        lcbd = response
    )[site$BCR == parameter$BCR, ]
    alpha <- seq_len(parameter$S)
    fit <- fit_parameterized_x(
        data,
        parameter$S,
        parameter$M,
        tolower(framework),
        x_bounds
    )

    parameterized_rows[[i]] <- data.frame(
        BCR = parameter$BCR, framework = framework,
        M = parameter$M, S = parameter$S,
        x = fit$x, mse = fit$mse
    )
    curve_rows[[2L * i - 1L]] <- data.frame(
        BCR = parameter$BCR, framework = framework, richness = alpha,
        expected = functions$parameterized(
            parameter$S, parameter$M, fit$x, alpha
        ),
        model = "Parameterized expectation"
    )
    annealed_curve <- complete_annealed_curve(
        parameter$S, parameter$M, parameter$x, parameter$c, alpha,
        functions$annealed, functions$approximation
    )
    curve_rows[[2L * i]] <- data.frame(
        BCR = parameter$BCR, framework = framework, richness = alpha,
        expected = annealed_curve$expected,
        model = "Annealed expectation",
        numerical_fallback = annealed_curve$numerical_fallback
    )
    curve_rows[[2L * i - 1L]]$numerical_fallback <- FALSE
}

parameterized_parameter <- do.call(rbind, parameterized_rows)
curves <- do.call(rbind, curve_rows)
rownames(curves) <- NULL

# 4. Save plotting data ---------------------------------------------------

saveRDS(
    list(
        settings = list(
            frameworks = names(model_functions),
            richness_support = "complete 1:S",
            fit_objective = "site-level MSE",
            numerical_fallback = paste(
                "The large-S approximation is used only when exact",
                "integration fails in the low-probability tail."
            )
        ),
        site_data = site,
        curves = curves,
        annealed_parameters = annealed_parameter,
        parameterized_parameters = parameterized_parameter
    ),
    output_file
)
