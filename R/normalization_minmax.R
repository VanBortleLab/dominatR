#' Min-Max Normalization
#'
#' @description
#' Scales each column of a matrix (or SummarizedExperiment assay) so that
#' the minimum value in that column is mapped to \code{new_min} and the
#' maximum value is mapped to \code{new_max}
#'
#'
#' @param x A numeric \code{matrix}, \code{data.frame}, or
#'   \code{SummarizedExperiment}.
#' @param new_min The lower bound of the new range (default 0).
#' @param new_max The upper bound of the new range (default 1).
#' @param assay_name If \code{x} is a SummarizedExperiment, name of the assay
#'   to normalize. Defaults to the first assay if none is specified.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new
#'   assay to store the normalized data. If \code{NULL}, overwrites the
#'   assay specified by \code{assay_name}.
#'
#' @return
#'   \itemize{
#'     \item If \code{x} is a data.frame or matrix, returns a matrix of
#'       column-wise scaled values (same dimensions as \code{x}).
#'     \item If \code{x} is a SummarizedExperiment, returns the same
#'       SummarizedExperiment object with the chosen or new assay replaced
#'       by the scaled values.
#'   }
#' @examples
#' library(SummarizedExperiment)
#' library(airway)
#' data('airway')
#'
#' se = airway
#'
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'
#' df = assay(se)
#'
#' df1 = minmax_normalization(df)
#'
#' apply(df1, 2, range)
#'
#' ## Using a new range
#' df1 = minmax_normalization(df, new_min = 5, new_max = 10)
#'
#' apply(df1, 2, range)
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#'
#' # If now new_assay_name is provided, then overwrites existing assay
#' se2 = minmax_normalization(se)
#'
#' apply(assay(se2), 2, range)
#'
#'
#' # If new new_assay_name, normalization stored in a new object
#' se2 = minmax_normalization(se, new_assay_name = 'minmax_counts')
#'
#' apply(assay(se2, 'minmax_counts'), 2, range)
#'
#' # A specific assay can also be selected
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
#' se2 = minmax_normalization(se,
#'                            new_assay_name = 'minmax_counts_new',
#'                            assay_name = 'new_counts')
#'
#' apply(assay(se2, 'minmax_counts_new'), 2, range)
#'
#' ## Using a different range
#' se2 = minmax_normalization(se,
#'                            new_assay_name = 'minmax_counts_new',
#'                            assay_name = 'new_counts',
#'                            new_min = 10,
#'                            new_max = 20)
#'
#' apply(assay(se2, 'minmax_counts_new'), 2, range)
#'
#' @export
minmax_normalization = function(x,
                                new_min = 0,
                                new_max = 1,
                                assay_name = NULL,
                                new_assay_name = NULL)
{
  if (inherits(x, "SummarizedExperiment")) {
    #-----------------------------#
    # SummarizedExperiment branch
    #-----------------------------#
    if (is.null(assay_name)) {
      # Default: use the first assay
      all_assays <- assayNames(x)
      if (length(all_assays) < 1) {
        stop("No assays found in the SummarizedExperiment.")
      }
      assay_name <- all_assays[[1]]
    }

    mat <- assay(x, assay_name)
    if (is.null(mat)) {
      stop("No assay named '", assay_name, "' found in the SummarizedExperiment.")
    }

    if (!is.numeric(mat)) {
      stop("Selected assay is not numeric.")
    }


  mins = apply(mat, 2, min, na.rm = TRUE)
  maxs = apply(mat, 2, max, na.rm = TRUE)

  mins = matrix(mins, nrow = nrow(mat), ncol = ncol(mat), byrow = T)
  maxs = matrix(maxs, nrow = nrow(mat), ncol = ncol(mat), byrow = T)

  diffs = maxs - mins

  new_rang = new_max - new_min

  mat = (mat - mins) * new_rang / diffs + new_min


  if (is.null(new_assay_name)) {
    # Overwrite existing assay
    assay(x, assay_name) <- mat
  } else {
    # Create a new assay
    assay(x, new_assay_name) <- mat
  }

  return(x)
  } else if (is.data.frame(x) || is.matrix(x)) {
  #-------------------#
  # data.frame/matrix
  #-------------------#
  # Convert data.frame to matrix if needed
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  # Ensure numeric
  if (!is.numeric(x)) {
    stop("Input data is not numeric.")
  }

    mat = x

    mins = apply(mat, 2, min, na.rm = TRUE)
    maxs = apply(mat, 2, max, na.rm = TRUE)

    mins = matrix(mins, nrow = nrow(mat), ncol = ncol(mat), byrow = T)
    maxs = matrix(maxs, nrow = nrow(mat), ncol = ncol(mat), byrow = T)

    diffs = maxs - mins

    new_rang = new_max - new_min

    mat = (mat - mins) * new_rang / diffs + new_min

    return(mat)

  } else {
    stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
  }
}
