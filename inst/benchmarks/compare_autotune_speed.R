# Run after installing CRAFT >= 0.4.0:
# Rscript compare_autotune_speed.R [output-directory]
library(CRAFT)
output_args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(output_args)) output_args[1L] else tempdir()
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(21)
series <- as.data.frame(replicate(3, cumprod(1 + rnorm(180, .001, .01))))
args <- list(series = series, n_models = 12, top_k = 2, n_testing = 5,
  search_space = list(transition_alpha = c(.25, .5, 1, 2),
                      sampler_weight = c("uniform", "bal_acc")),
  fixed_args = list(window = 10, n_draws = 300, n_factors = 2,
                    min_segment = 10, max_regimes_per_factor = 4),
  seed = 42, verbose = FALSE)
# Warm up namespace loading and byte compilation outside measured runs.
warm <- args
warm$n_models <- warm$top_k <- warm$n_testing <- 1L
invisible(do.call(craft_autotune, warm))
measure <- function(method, fun, extra = list()) {
  elapsed <- system.time(value <- do.call(fun, c(args, extra)))[["elapsed"]]
  origin_fits <- sum(value$performance$origins_computed)
  list(value = value, row = data.frame(method = method, elapsed = elapsed,
    origin_fits = origin_fits, winner = value$top_models[[1]]$model_id))
}
rows <- list()
for (repeat_id in seq_len(3)) {
  baseline <- measure("uncached_full_budget", craft_autotune, list(stage_budgets = 5, cache = FALSE))
  exact <- measure("cached_full_budget", craft_autotune, list(stage_budgets = 5))
  staged <- measure("staged", craft_autotune)
  stopifnot(isTRUE(all.equal(baseline$value$results$score, exact$value$results$score)))
  for (item in list(baseline, exact, staged)) {
    row <- item$row
    row$repeat_id <- repeat_id
    rows[[length(rows) + 1L]] <- row
  }
}
results <- do.call(rbind, rows)
print(results)
utils::write.csv(results, file.path(output_dir, "autotune_speed.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "autotune_session_info.txt"))
