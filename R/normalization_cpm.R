#' Counts Per Million normalization
#'
#' @description
#' Normalizes a count matrix (or a SummarizedExperiment assay) by the
#' counts-per-million (CPM) method. Specifically:
#'   \enumerate{
#'     \item If \code{log_trans = TRUE}, a \code{log2(x + 1)} transform is
#'     applied afterward.
#'   }
#'
#' @param x A \code{matrix}, \code{data.frame}, or a \code{SummarizedExperiment}
#'   object.
#' @param log_trans Logical. If \code{TRUE}, apply \code{log2(... + 1)}
#'   transform to the CPM-normalized values.
#' @param assay_name If \code{x} is a \code{SummarizedExperiment}, name of the
#'   assay to normalize (defaults to the first assay). Ignored otherwise.
#' @param new_assay_name If \code{x} is a \code{SummarizedExperiment}, name of
#'   a new assay where results should be stored (defaults to \code{NULL},
#'   meaning the existing assay is overwritten).
#'
#' @return
#'   \itemize{
#'     \item If \code{x} is a \code{matrix} or \code{data.frame}, returns a
#'       \strong{matrix} of CPM-normalized (and optionally
#'       \code{log2}-transformed) counts.
#'     \item If \code{x} is a \code{SummarizedExperiment}, returns the same
#'       \code{SummarizedExperiment} object with the specified assay replaced
#'       or a new assay created containing the CPM-normalized data.
#'   }
#'
#' @importFrom SummarizedExperiment assay assayNames SummarizedExperiment
#' @importFrom SummarizedExperiment assay<- rowData<-
#'
#' @examples
#'
#'library(SummarizedExperiment)
#'library(airway)
#'data('airway')
#'
#'se = airway
#'
#' # Only use a random subset of 1000 rows
#' set.seed(123)
#' idx <- sample(seq_len(nrow(se)), size = min(1000, nrow(se)))
#' se <- se[idx, ]
#'
#'# -------------------------------
#'# 1) Using a data.frame
#'# -------------------------------
#'
#'df = assay(se)
#'
#'## Without log transformation
#'df1 = cpm_normalization(df, log_trans = FALSE)
#'
#'df1[1:5,1:5]
#'
#'## With log transformation
#'df1 = cpm_normalization(df, log_trans = TRUE)
#'
#'df1[1:5,1:5]
#'
#'# -------------------------------
#'# 2) Using a SummarizedExperiment
#'# -------------------------------
#'
#'# If now new_assay_name is provided, then overwrites existing assay
#'se2 = cpm_normalization(se, log_trans = FALSE)
#'
#'se2
#'head(assay(se2))
#'
#'# If new new_assay_name, normalization stored in a new object
#'se2 = cpm_normalization(se, log_trans = FALSE, new_assay_name = 'cpm_counts')
#'
#'se2
#'head(assay(se2, 'cpm_counts'))
#'
#'# A specific assay can also be selected
#'new_matrix =  matrix(data = sample(x = seq(1, 100000),
#'                                   size = nrow(se) * ncol(se),
#'                                   replace = TRUE),
#'                     nrow = nrow(se),
#'                     ncol = ncol(se))
#'rownames(new_matrix) = rownames(se)
#'colnames(new_matrix) = colnames(se)
#'
#'## Creating a new assay called new counts
#'assay(se, 'new_counts') = new_matrix
#'
#'se2 = cpm_normalization(se, new_assay_name = 'cpm_counts_new', assay_name =
#' 'new_counts')
#'
#'se2
#'head(assay(se2, 'cpm_counts_new'))
#' @export

cpm_normalization <- function(x,
                                log_trans       = FALSE,
                                assay_name      = NULL,
                                new_assay_name  = NULL) {
    if (inherits(x, "SummarizedExperiment")) {
        #-----------------------------#
        # SummarizedExperiment branch
        #-----------------------------#
        m <- .get_matrix(se = x, a_name = assay_name)

        mat <- m$mat
        assay_name <- m$assay_name

        col_millions <- colSums(mat, na.rm = TRUE) / 1e6

        col_millions[col_millions <= 0] <- 1

        mat_cpm <- sweep(mat, 2, col_millions, FUN="/")

        if (log_trans) {
            mat_cpm <- log2(mat_cpm + 1)
        }

        x <- .assing_assay(se = x, matrix = mat_cpm, a_name = assay_name,
                new_a_name = new_assay_name)

        return(x)

    } else if (is.data.frame(x) || is.matrix(x)) {
    #-------------------#
    # data.frame/matrix
    #-------------------#
    # Convert data.frame to matrix if needed

    x <- .get_matrix_df(x)

    col_millions <- colSums(x, na.rm = TRUE) / 1e6
    col_millions[col_millions <= 0] <- 1

    x_cpm <- sweep(x, 2, col_millions, FUN="/")

    if (log_trans) {
        x_cpm <- log2(x_cpm + 1)
    }

    return(x_cpm)

    } else {
        stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
    }
}

