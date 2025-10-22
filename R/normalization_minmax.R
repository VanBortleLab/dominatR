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
#'
#' @importFrom SummarizedExperiment assay assayNames SummarizedExperiment
#' @importFrom SummarizedExperiment assay<- rowData<-
#'
#' @examples
#' library(SummarizedExperiment)
#' library(airway)
#' data('airway')
#'
#' se <- airway
#'
#' # Only use a random subset of 1000 rows
#' set.seed(123)
#' idx <- sample(seq_len(nrow(se)), size = min(1000, nrow(se)))
#' se <- se[idx, ]
#'
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'
#' df <- assay(se)
#'
#' df1 <- minmax_normalization(df)
#'
#' apply(df1, 2, range)
#'
#' ## Using a new range
#' df1 <- minmax_normalization(df, new_min = 5, new_max = 10)
#'
#' apply(df1, 2, range)
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#'
#' # If now new_assay_name is provided, then overwrites existing assay
#' se2 <- minmax_normalization(se)
#'
#' apply(assay(se2), 2, range)
#'
#'
#' # If new new_assay_name, normalization stored in a new object
#' se2 <- minmax_normalization(se, new_assay_name = 'minmax_counts')
#'
#' apply(assay(se2, 'minmax_counts'), 2, range)
#'
#' # A specific assay can also be selected
#' new_matrix <-  matrix(data = sample(x = seq(1, 100000),
#'                                    size = nrow(se) * ncol(se),
#'                                    replace = TRUE),
#'                      nrow = nrow(se),
#'                      ncol = ncol(se))
#' rownames(new_matrix) <- rownames(se)
#' colnames(new_matrix) <- colnames(se)
#'
#' ## Creating a new assay called new counts
#' assay(se, 'new_counts') <- new_matrix
#'
#' se2 <- minmax_normalization(se,
#'                            new_assay_name = 'minmax_counts_new',
#'                            assay_name = 'new_counts')
#'
#' apply(assay(se2, 'minmax_counts_new'), 2, range)
#'
#' ## Using a different range
#' se2 <- minmax_normalization(se,
#'                            new_assay_name = 'minmax_counts_new',
#'                            assay_name = 'new_counts',
#'                            new_min = 10,
#'                            new_max = 20)
#'
#' apply(assay(se2, 'minmax_counts_new'), 2, range)
#'
#' @export
minmax_normalization <- function(x,
                                new_min = 0,
                                new_max = 1,
                                assay_name = NULL,
                                new_assay_name = NULL)
{
    if (inherits(x, "SummarizedExperiment")) {
    #-----------------------------#
    # SummarizedExperiment branch
    #-----------------------------#
    m <- .get_matrix(se = x, a_name = assay_name)

    mat <- m$mat
    assay_name <- m$assay_name

    mins <- apply(mat, 2, min, na.rm = TRUE)
    maxs <- apply(mat, 2, max, na.rm = TRUE)

    mins <- matrix(mins, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
    maxs <- matrix(maxs, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)

    diffs <- maxs - mins

    new_rang <- new_max - new_min

    mat <- (mat - mins) * new_rang / diffs + new_min


    x <- .assing_assay(se = x, matrix = mat, a_name = assay_name,
                new_a_name = new_assay_name)

    return(x)
    } else if (is.data.frame(x) || is.matrix(x)) {
    #-------------------#
    # data.frame/matrix
    #-------------------#
    # Convert data.frame to matrix if needed
    mat <- .get_matrix_df(x)


    mins <- apply(mat, 2, min, na.rm = TRUE)
    maxs <- apply(mat, 2, max, na.rm = TRUE)

    mins <- matrix(mins, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
    maxs <- matrix(maxs, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)

    diffs <- maxs - mins

    new_rang <- new_max - new_min

    mat <- (mat - mins) * new_rang / diffs + new_min

    return(mat)

    } else {
        stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
    }
}
