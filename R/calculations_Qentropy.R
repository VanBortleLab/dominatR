#' Compute Q-Entropy using existing row-normalized data + Entropy
#'
#' @description
#' #' Transform entropy scores into categorical entropy scores
#' \eqn{Q_{ij} = \mathrm{Entropy}_i - \log_2(x_{ij})}, or \code{Inf} if
#' \eqn{x_{ij} == 0}.
#'
#'  @details
#' For each row \eqn{i} and column \eqn{j}, \eqn{Q_{ij}} is defined as
#' \eqn{\mathrm{Entropy}_i - \log_2\bigl(x_{ij}\bigr)} if \eqn{x_{ij}} is
#' positive, or \code{Inf} otherwise.
#'
#' @param x A data.frame (already processed by `entropy()`)
#'   or a SummarizedExperiment (already processed by `entropy()`).
#' @param assay_name (SummarizedExperiment only) The name of the assay
#'   whose row-normalized data will be replaced by Q-values. If NULL,
#'   uses the first assay.
#' @param new_assay_name If you prefer to store Q-values in a
#'   *new* assay, provide a name. By default `Qentropy`
#'
#' @return
#'  \itemize{
#'   \item If \code{x} is a data.frame: returns the same data.frame with numeric
#'     columns replaced by \eqn{Q_{ij}} values and \code{Entropy} column removed.
#'   \item If \code{x} is a SummarizedExperiment: returns the same object with
#'     the specified assay replaced by \eqn{Q_{ij}} values (or a new assay
#'     if \code{new_assay_name} is set) and \code{rowData(x)$Entropy} removed.
#'     }
#'
#' @import SummarizedExperiment
#'
#' @examples
#'
#'library(SummarizedExperiment)
#'library(airway)
#'data('airway')
#'
#'se = airway
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'df = assay(se) |> as.data.frame()
#'
#'## Entropy needs to be calculated first
#'df = entropy(df)
#'
#'## Then you can apply the Qentropy function
#'df = Qentropy(df)
#'
#'head(df)
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#'
#'## Calculate Entropy first
#'se2 = entropy(se, new_assay_name = 'Entropy')
#'
#'## Transform entropy into Qentropy. new_assay_name specify a new assay
#'## where data is going to be stored. Assay_name must have Entropy transformed values
#'## By default function will look for an assay_name `Entropy` and assign a new assay to `Qentropy`
#'se2 = Qentropy(se2, new_assay_name = 'Qentropy', assay_name = 'Entropy')
#'
#'se2
#'

#'@export


Qentropy <- function(x,
                      assay_name = 'Entropy',
                      new_assay_name = 'Qentropy') {
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

    # Must have rowData(x)$Entropy
    if (!("Entropy" %in% colnames(rowData(x)))) {
      stop("rowData(x)$Entropy not found. Please run `entropy()` first.")
    }
    e_vals <- rowData(x)$Entropy

    # Q_{ij} = Inf if mat[i,j] == 0, else e_vals[i] - log2(mat[i,j])
    q_mat <- matrix(Inf, nrow = nrow(mat), ncol = ncol(mat))

    non_zero <- mat > 0
    # index of non-zero entries
    idx <- which(non_zero, arr.ind = TRUE)
    # row index
    i_idx <- idx[,1]
    # col index
    j_idx <- idx[,2]

    # e_vals[i] - log2(mat[i,j])
    q_mat[non_zero] <- e_vals[i_idx] - log2(mat[non_zero])

    rownames(q_mat) <- rownames(mat)
    colnames(q_mat) <- colnames(mat)

    # Overwrite or store in a new assay
    if (is.null(new_assay_name)) {
      # Overwrite the same assay
      assay(x, assay_name) <- q_mat
    } else {
      assay(x, new_assay_name) <- q_mat
    }

    return(x)

  } else if (is.data.frame(x)) {
    #--------------#
    # data.frame
    #--------------#
    # Must have 'Entropy'
    if (!("Entropy" %in% colnames(x))) {
      stop("'Entropy' column not found in data.frame. Did you run `entropy()` first?")
    }

    # Numeric columns (excluding 'Entropy')
    numeric_cols <- vapply(x, is.numeric, logical(1))
    numeric_cols["Entropy"] <- FALSE

    if (!any(numeric_cols)) {
      # If we have no numeric cols other than Entropy, just drop Entropy
      x[["Entropy"]] <- NULL
      return(x)
    }

    # The data in these numeric columns should be the row-normalized values
    mat <- as.matrix(x[, numeric_cols, drop = FALSE])
    e_vals <- x[["Entropy"]]

    # Build Q matrix
    q_mat <- matrix(Inf, nrow = nrow(mat), ncol = ncol(mat))
    non_zero <- mat > 0
    idx <- which(non_zero, arr.ind = TRUE)
    i_idx <- idx[,1]
    j_idx <- idx[,2]

    q_mat[non_zero] <- e_vals[i_idx] - log2(mat[non_zero])

    # Overwrite numeric columns
    x[, numeric_cols] <- q_mat

    x[["Entropy"]] <- NULL

    return(x)
  } else {
    stop("Input must be a data.frame or SummarizedExperiment.")
  }
}
