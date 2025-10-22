#' Quantile Normalization
#'
#' @description
#' Normalizes read counts by the quantile normalization method:
#'   \enumerate{
#'     \item Each sample (column) is sorted, and values at each rank are
#'     averaged across columns
#'     \item Each sample's values are replaced with the average of their
#'     respective rank
#'     \item If \code{log_trans = TRUE}, applies \code{log2(QN + 1)}
#'     transformation
#'   }
#'
#' @details
#' If \code{x} is a \code{SummarizedExperiment}, the function will extract the
#' assay using \code{assay_name}, apply quantile normalization, and return a
#' new or updated assay. If \code{x} is a matrix or data.frame, normalization is
#' applied directly to the input matrix.
#'
#' @param x A numeric \code{matrix} or \code{data.frame} of gene counts,
#'   or a \code{SummarizedExperiment} containing such counts.
#'   \describe{
#'     \item{If a \code{SummarizedExperiment},}{the function applies
#'       normalization to the specified assay (via \code{assay_name}).}
#'     \item{If a \code{data.frame}/\code{matrix},}{the normalization is
#'     applied directly.}
#'   }
#' @param log_trans Logical. If \code{TRUE}, apply \code{log2(... + 1)}
#'   transform to the quantile-normalized values.
#' @param assay_name If \code{x} is a SummarizedExperiment, name of the assay
#'   to normalize. Defaults to the first assay if not specified.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new
#'   assay in which to store the quantile-normalized (or log2-transformed)
#'   values. If \code{NULL}, overwrites the original assay.
#'
#' @return A numeric \strong{matrix} of quantile-normalized (or log2-normalized)
#'   values if \code{x} is a data.frame or matrix. If \code{x} is a
#'   SummarizedExperiment, returns the modified SummarizedExperiment with the
#'   normalized data placed in the existing or new assay.
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
#' df <- assay(se)
#'
#' ## Without log transformation
#' df_qn <- quantile_normalization(df, log_trans = FALSE)
#' df_qn[1:5, 1:5]
#'
#' ## With log transformation
#' df_qn_log <- quantile_normalization(df, log_trans = TRUE)
#' df_qn_log[1:5, 1:5]
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#'
#' ## Overwrite existing assay
#' se2 <- quantile_normalization(se)
#' assay(se2)[1:5, 1:5]
#'
#' ## Store result in new assay
#' se3 <- quantile_normalization(se, new_assay_name = "quant_norm")
#' assay(se3, "quant_norm")[1:5, 1:5]
#'
#' ## Use specific input assay (simulate new one)
#' new_matrix <- matrix(
#'   data = sample(x = seq(1, 100000), size = nrow(se) * ncol(se),
#'   replace = TRUE),
#'   nrow = nrow(se),
#'   ncol = ncol(se)
#' )
#' rownames(new_matrix) <- rownames(se)
#' colnames(new_matrix) <- colnames(se)
#'
#' ## Create a new assay in the SummarizedExperiment
#' assay(se, "new_counts") <- new_matrix
#'
#' ## Normalize the new assay and store it under a new name
#' se4 <- quantile_normalization(se, assay_name = "new_counts",
#' new_assay_name = "quant_new")
#' assay(se4, "quant_new")[1:5, 1:5]
#' @export
#'
quantile_normalization <- function(x,
                                    log_trans       = FALSE,
                                    assay_name      = NULL,
                                    new_assay_name  = NULL) {

    #---------------------------#
    # SummarizedExperiment path
    #---------------------------#
    if (inherits(x, "SummarizedExperiment")) {

        m <- .get_matrix(se = x, a_name = assay_name)

        mat <- m$mat
        assay_name <- m$assay_name

    ## Calculations
    # 1) Rank matrix
    ranks <- apply(mat, 2, rank, ties.method = "min")

    # 2) Sorted matrix
    sorted <- apply(mat, 2, sort)

    # 3) Row means
    means <- rowMeans(sorted)

    # 4) Map means to original ranks
    norm <- apply(ranks, 2, function(r) means[r])

    # 5) Preserve row/col names
    rownames(norm) <- rownames(mat)
    colnames(norm) <- colnames(mat)

    # 6) Optional log2
    if (log_trans) {
        norm <- log2(norm + 1)
    }

    # 7) Store result in new or existing assay

    x <- .assing_assay(se = x, matrix = norm, a_name = assay_name,
                new_a_name = new_assay_name)

    return(x)

    #---------------------------#
    # data.frame / matrix path
    #---------------------------#
    } else if (is.data.frame(x) || is.matrix(x)) {

    x <- .get_matrix_df(x)

    ## Calculations
    # 1) Rank matrix
    ranks <- apply(x, 2, rank, ties.method = "min")

    # 2) Sorted matrix
    sorted <- apply(x, 2, sort)

    # 3) Row means
    means <- rowMeans(sorted)

    # 4) Map means to original ranks
    norm <- apply(ranks, 2, function(r) means[r])

    # 5) Preserve row/col names
    rownames(norm) <- rownames(x)
    colnames(norm) <- colnames(x)

    # 6) Optional log2
    if (log_trans) {
        norm <- log2(norm + 1)
    }

    return(norm)

    } else {
        stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
    }
}
