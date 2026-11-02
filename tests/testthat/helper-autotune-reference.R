# Independent full-pipeline reference: every candidate builds distributions
# and runs craft_fit's ordinary backtests, without the tuning engine's caches
# or score-only path. Used only to verify the accelerated scoring path.
reference_autotune <- function(series, n_models, top_k, n_testing,
                                search_space, fixed_args, seed = 42,
                                verbose = FALSE, ...) {
  set.seed(seed)
  params <- lapply(seq_len(n_models), function(i) CRAFT:::.craft_autotune_draw_params(search_space))
  fits <- lapply(seq_len(n_models), function(i) {
    args <- utils::modifyList(fixed_args, params[[i]], keep.null = TRUE)
    args$series <- series
    args$seed <- as.integer((as.double(seed) + i) %% .Machine$integer.max)
    args$n_testing <- n_testing
    args$backtest_parallel <- FALSE
    args$return_train_probs <- FALSE
    args$verbose <- verbose
    do.call(craft_fit, args)
  })
  rows <- lapply(seq_len(n_models), function(i) {
    scored <- CRAFT:::.craft_autotune_score_fit(fits[[i]])
    cbind(data.frame(model_id = i, ok = TRUE, score = scored$score,
                     error = NA_character_), scored$overall,
          CRAFT:::.craft_autotune_param_row(params[[i]]))
  })
  results <- do.call(rbind, rows)
  selected <- utils::head(order(results$score), top_k)
  list(results = results, top_models = lapply(selected, function(i) {
    list(model_id = i, fit = fits[[i]], score = results$score[i])
  }))
}
