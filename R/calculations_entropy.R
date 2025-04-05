#' Compute Shannon Entropy on row-normalized data
#'
#' @param x A data.frame (with numeric columns)
#'   or a SummarizedExperiment (with an assay of numeric data).
#' @param assay_name (SummarizedExperiment only) The name of the assay
#'   to transform and compute Entropy on. If NULL, uses the first assay.
#' @param new_assay_name If you prefer to store Q-values in a
#'   *new* assay, provide a name. By default `Entropy`
#'
#' @import SummarizedExperiment airway
#'
#'
#' @return
#'  \itemize{
#'   \item If \code{x} is a data.frame: returns the same data.frame in which
#'     numeric columns have been replaced by their row-wise proportions,
#'     and an \code{Entropy} column is appended.
#'   \item If \code{x} is a SummarizedExperiment: returns the same
#'     SummarizedExperiment in with a new assay (Default name is \code{Entropy})
#'     and \code{rowData(x)$Entropy} is added.
#'     }
#'
#' @examples
#' library(SummarizedExperiment)
#' library(airway)
#' data('airway')
#'
#' se = airway
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#' df = assay(se) |> as.data.frame()
#' df = entropy(df)
#'
#' ## The function adds a new column called Entropy and transform all the counts accordingly
#' head(df)
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#'
#' ## The function adds a new assay called 'Entropy' with the transformed counts.
#' ## This name can be modified with the `new_assay_name` parameter
#' ## In the rowData dataframe a new column called Entropy is added.
#' se2 = entropy(se, new_assay_name = 'Entropy')
#' se2
#'
#' ## In case the experiment has multiple assays, the function allows you to choose which assay to use.
#' new_matrix =  matrix(data = sample(x = seq(1, 100000),
#'                                    size = nrow(se) * ncol(se),
#'                                    replace = TRUE),
#'                      nrow = nrow(se),
#'                      ncol = ncol(se))
#' rownames(new_matrix) = rownames(se)
#' colnames(new_matrix) = colnames(se)
#'
#' ## Creating a new assay called new counts
#' assay(se, 'new_counts') = new_matrix
#'
#'
#' ## Saving the entropy values as Entropy_newmatrix using the assay 'new counts'
#' se2 = entropy(se,
#'               new_assay_name = 'Entropy_newmatrix',
#'               assay_name = 'new_counts')
#'
#' se2
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
