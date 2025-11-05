#' Internal function to compute Anti-Robinson or GAR/ RGAR score
#'
#' This function is not exported. It serves as the core computational engine
#' for evaluating Anti-Robinson (AR), Generalized Anti-Robinson (GAR), and
#' Relative GAR (RGAR) scores.
#'
#' @param mat_sorted A numeric, symmetric sorted distance matrix.
#' @param w An integer window size (if NULL, evaluates all triplets globally).
#' @param normalize Logical. If TRUE, returns a proportion; otherwise returns the raw count.
#'
#' @return A numeric value indicating the number (or proportion) of Anti-Robinson violations.
#' @keywords internal
# ---- Internal core function ----
.compute_gar_core <- function(mat_sorted, w = NULL, normalize = FALSE) {
  n <- nrow(mat_sorted)
  score <- 0
  count <- 0

  for (i in 1:n) {
    # 左側 i-w <= j < k < i
    j_start <- max(1, if (is.null(w)) 1 else i - w)
    j_end <- i - 2
    if (j_start <= j_end) {
      for (j in j_start:j_end) {
        k_start <- j + 1
        k_end <- i - 1
        if (k_start <= k_end) {
          for (k in k_start:k_end) {
            # 保證 j, k 在合法範圍
            if (j <= n && k <= n)
              score <- score + as.integer(mat_sorted[i, j] < mat_sorted[i, k])
            count <- count + 1
          }
        }
      }
    }

    # 右側 i < j < k <= i + w
    j_start <- i + 1
    j_end <- min(n - 1, if (is.null(w)) n - 1 else i + w)
    if (j_start <= j_end) {
      for (j in j_start:j_end) {
        k_start <- j + 1
        k_end <- min(n, if (is.null(w)) n else i + w)
        if (k_start <= k_end) {
          for (k in k_start:k_end) {
            if (j <= n && k <= n)
              score <- score + as.integer(mat_sorted[i, j] > mat_sorted[i, k])
            count <- count + 1
          }
        }
      }
    }
  }

  if (normalize) {
    return(if (count > 0) score / count else NA_real_)
  } else {
    return(score)
  }
}


#' Compute the Anti-Robinson (AR) score
#'
#' Calculates the total number of Anti-Robinson violations over all triplets
#' in the matrix using the specified ordering. This is equivalent to GAR with a full window.
#'
#' @param mat_sorted A numeric, symmetric sorted distance matrix.
#'
#' @return The AR score (the total number of structural violations).
#'
#' Please refer to \code{\link{GAP}} for complete usage examples.
#' @export
# ---- AR = GAR with w = NULL, no normalization ----
AR <- function(mat_sorted) {
  .compute_gar_core(mat_sorted, w = NULL, normalize = FALSE)
}

#' Compute the Generalized Anti-Robinson (GAR) score
#'
#' Calculates the number of Anti-Robinson violations within a specified window `w`,
#' allowing evaluation of local structural consistency in the reordered matrix.
#'
#' @param mat_sorted A numeric, symmetric sorted distance matrix.
#' @param w Window size (integer). If NULL, uses global comparisons (equivalent to AR).
#'
#' @return The GAR score (the total number of violations).
#'
#' Please refer to \code{\link{GAP}} for complete usage examples.
#' @export
# ---- GAR with specified w (default NULL = global) ----
GAR <- function(mat_sorted, w = NULL) {
  .compute_gar_core(mat_sorted, w = w, normalize = FALSE)
}

#' Compute the Relative Generalized Anti-Robinson (RGAR) score
#'
#' This function returns the relative GAR score, representing the proportion of
#' Anti-Robinson violations over the total number of evaluated triplets.
#'
#' @param mat_sorted A numeric, symmetric sorted distance matrix.
#' @param w Window size (integer). If NULL, uses global comparisons (equivalent to AR normalized).
#'
#' @return The RGAR score (between 0 and 1).
#'
#' Please refer to \code{\link{GAP}} for complete usage examples.
#' @export
# ---- RGAR = normalized GAR ----
RGAR <- function(mat_sorted, w = NULL) {
  .compute_gar_core(mat_sorted, w = w, normalize = TRUE)
}
