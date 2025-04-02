#' Compute Shannon Entropy on row-normalized data
#'
#' @param x A data.frame (with numeric columns)
#'   or a SummarizedExperiment (with an assay of numeric data).
#' @param assay_name (SummarizedExperiment only) The name of the assay
#'   to transform and compute Entropy on. If NULL, uses the first assay.
#' @param new_assay_name If you prefer to store Q-values in a
#'   *new* assay, provide a name. By default `Entropy`
#'
#' @import SummarizedExperiment
#'
#'
#' @return
#'   - If `x` is a data.frame: returns the same data.frame in which
#'     numeric columns have been replaced by their row-wise proportions,
#'     and an `Entropy` column is appended.
#'   - If `x` is a SummarizedExperiment: returns the same
#'     SummarizedExperiment in with a new assay (Default name is `Entropy`)
#'     and `rowData(x)$Entropy` is added.
#'
#' @examples
#' df <- data.frame(A = c(10, 0, 30),
#'                  B = c(5,  5, 0),
#'                  C = c(0,  0, 10))
#' # Original function row-normalizes the numeric columns, then calculates Entropy
#' df2 <- entropy(df)
#' df2
#'
#' # SummarizedExperiment example
#'  library(SummarizedExperiment)
#'  mat <- matrix(c(10,0,30, 5,5,0, 0,0,10), nrow=3, byrow=TRUE)
#'  se <- SummarizedExperiment(assays = list(counts=mat))
#'  se2 <- entropy(se)
#'  assay(se2, "counts")  # row-normalized
#'  rowData(se2)$Entropy  # computed entropy
#'
#'@export

entropy <- function(x,
                    assay_name = NULL,
                    new_assay_name = 'Entropy') {
  if (inherits(x, "SummarizedExperiment")) {
    #----------------------#
    # SummarizedExperiment
    #----------------------#
    if (is.null(assay_name)) {
      assay_name <- assayNames(x)[1]
      if (is.na(assay_name)) {
        stop("No assay found in SummarizedExperiment.")
      }
    }
    mat <- assay(x, assay_name)
    if (is.null(mat)) {
      stop("No assay named '", assay_name, "' found in the SummarizedExperiment.")
    }

    # 1. Compute row sums (use 1 if sum <= 0)
    row_sums <- rowSums(mat, na.rm = TRUE)
    row_sums <- ifelse(row_sums > 0, row_sums, 1)

    # 2. Row-normalize the assay
    mat_norm <- mat / row_sums

    # 3. Compute Shannon entropy:
    #    We'll do it vectorized:
    log_mat <- matrix(0, nrow = nrow(mat_norm), ncol = ncol(mat_norm))
    non_zero <- mat_norm > 0
    log_mat[non_zero] <- log2(mat_norm[non_zero])

    # row-wise sum of p_ij*log2(p_ij)
    ent_vals <- -rowSums(mat_norm * log_mat)

    # 4. Store the normalized assay back
    assay(x, new_assay_name) <- mat_norm

    # 5. Save Entropy in rowData
    rowData(x)$Entropy <- ent_vals

    return(x)

  } else if (is.data.frame(x)) {
    #--------------#
    # data.frame
    #--------------#
    # 1. Identify numeric columns
    numeric_cols <- vapply(x, is.numeric, logical(1))
    # If no numeric columns, nothing to do
    if (!any(numeric_cols)) {
      warning("No numeric columns found; returning data unchanged.")
      return(x)
    }

    # 2. Row sums
    mat <- as.matrix(x[, numeric_cols, drop = FALSE])
    row_sums <- rowSums(mat, na.rm = TRUE)
    row_sums <- ifelse(row_sums > 0, row_sums, 1)

    # 3. Row-normalize
    mat_norm <- mat / row_sums

    # 4. Compute Shannon entropy
    log_mat <- matrix(0, nrow = nrow(mat_norm), ncol = ncol(mat_norm))
    non_zero <- mat_norm > 0
    log_mat[non_zero] <- log2(mat_norm[non_zero])
    ent_vals <- -rowSums(mat_norm * log_mat)

    # 5. Replace numeric columns in x with normalized data
    x[, numeric_cols] <- mat_norm

    # 6. Add 'Entropy' column
    x[["Entropy"]] <- ent_vals

    return(x)
  } else {
    stop("Input must be a data.frame or SummarizedExperiment.")
  }
}
