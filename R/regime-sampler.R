regime_sampler <- function(
        trafo_df,
        cluster_df,
        prob_vectors,
        weights = NULL,
        size = 1000L,
        return_type = c("stack", "index"),
        seed = NULL,
        eps = 1e-12,
        ...
) {
    dots <- list(...)
    if (!is.null(dots[["return"]])) return_type <- dots[["return"]]
    return_type <- match.arg(return_type)
    if (!is.null(seed)) set.seed(seed)
    
    trafo_df <- as.data.frame(trafo_df, stringsAsFactors = FALSE)
    cluster_df <- as.data.frame(cluster_df, stringsAsFactors = FALSE)
    
    if (nrow(trafo_df) != nrow(cluster_df)) {
        stop("`trafo_df` and `cluster_df` must have the same number of rows.")
    }
    
    if (!nrow(trafo_df)) stop("Cannot sample from an empty `trafo_df`.")
    
    target_names <- names(cluster_df)
    if (is.null(target_names) || any(!nzchar(target_names))) {
        target_names <- paste0("target_", seq_len(ncol(cluster_df)))
        names(cluster_df) <- target_names
    }
    
    if (is.null(weights)) weights <- rep(1, length(target_names))
    weights <- as.numeric(weights)
    if (length(weights) != length(target_names)) weights <- rep(1, length(target_names))
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) weights[] <- 1
    weights <- weights / sum(weights)
    names(weights) <- target_names
    
    log_score <- rep(0, nrow(cluster_df))
    
    for (target in target_names) {
        pv <- prob_vectors[[target]]
        if (is.null(pv)) next
        
        if (is.data.frame(pv) || is.matrix(pv)) {
            pv <- pv[1, , drop = TRUE]
        }
        
        pv <- as.numeric(pv)
        nm <- colnames(prob_vectors[[target]])
        if (is.null(nm)) nm <- names(prob_vectors[[target]])
        if (is.null(nm) && !is.null(colnames(as.matrix(prob_vectors[[target]])))) nm <- colnames(as.matrix(prob_vectors[[target]]))
        if (is.null(nm) || length(nm) != length(pv)) nm <- as.character(sort(unique(cluster_df[[target]])))
        
        names(pv) <- nm
        pv <- .mc_normalize_prob(pv)
        
        states <- as.character(cluster_df[[target]])
        p_row <- pv[states]
        p_row[!is.finite(p_row)] <- eps
        log_score <- log_score + weights[target] * log(pmax(p_row, eps))
    }
    
    prob <- .mc_softmax(log_score)
    idx <- sample.int(nrow(trafo_df), size = as.integer(size), replace = TRUE, prob = prob)
    
    if (return_type == "index") return(idx)
    out <- trafo_df[idx, , drop = FALSE]
    rownames(out) <- NULL
    out
}
