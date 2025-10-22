#' RPKM Normalization
#'
#' @description
#' Normalizes read counts by the RPKM (Reads Per Kilobase per Million
#' mapped reads) method:
#' \enumerate{
#'   \item Normalize counts by library size (column sums), scaled to millions.
#'   \item Divide each gene's value by its length in kilobases.
#'   \item If \code{log_trans = TRUE}, applies \code{log2(RPKM + 1)}.
#' }
#'
#' @details
#' If \code{x} is a \code{SummarizedExperiment}, the function looks for a
#' numeric column named \code{"gene_length"} in \code{rowData(x)}. That column
#' must have length equal to the number of rows in the assay being normalized.
#'
#' @param x A numeric \code{matrix} or \code{data.frame} of gene counts,
#'   or a \code{SummarizedExperiment} containing such counts.
#'   \describe{
#'     \item{If a \code{SummarizedExperiment},}{the function retrieves
#'       \code{gene_length} from \code{rowData(x)$gene_length}.}
#'     \item{If a \code{data.frame}/\code{matrix},}{the user must provide the
#'       \code{gene_length} argument.}
#'   }
#' @param gene_length A numeric vector of gene lengths (one per row), used
#'   only if \code{x} is a data.frame or matrix. Must match the number of rows
#'   in \code{x}. Ignored if \code{x} is a SummarizedExperiment.
#' @param log_trans Logical. If \code{TRUE}, apply \code{log2(... + 1)}
#'   transform to the RPKM-normalized values.
#' @param assay_name If \code{x} is a SummarizedExperiment, name of the assay to
#'   normalize. Defaults to the first assay if not specified.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new
#'   assay in which to store the RPKM (or log2-RPKM). If \code{NULL}, overwrites
#'   the assay specified in \code{assay_name}.
#'
#' @return A numeric \strong{matrix} of RPKM or log2(RPKM + 1) values if
#'   \code{x} is a data.frame or matrix. If \code{x} is a SummarizedExperiment,
#'   returns the modified SummarizedExperiment with the RPKM data placed in the
#'   existing or new assay.
#'
#' @importFrom SummarizedExperiment assay assayNames SummarizedExperiment
#' @importFrom SummarizedExperiment assay<- rowData<-
#'
#' @examples
#' library(SummarizedExperiment)
#' library(airway)
#' data('airway')
#'
#' se = airway
#'
#' # Only use a random subset of 1000 rows
#' set.seed(123)
#' idx <- sample(seq_len(nrow(se)), size = min(1000, nrow(se)))
#' se <- se[idx, ]
#'
#' ### Adding a column in rowData regarding the gene_length
#' rowData(se)$gene_length = rowData(se)$gene_seq_end -
#' rowData(se)$gene_seq_start
#'
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'
#' gene_length = rowData(se)$gene_length
#' df = assay(se)
#'
#' ## Without log transformation
#' df = rpkm_normalization(df, gene_length = gene_length)
#' df[1:5, 1:5]
#'
#' ## With log transformation
#' df = rpkm_normalization(df, gene_length = gene_length, log_trans = TRUE)
#' df[1:5, 1:5]
#'
#' # -------------------------------
#' # 2) Using a SummarizedExperiment
#' # -------------------------------
#'
#' # If no new_assay_name is provided, then overwrites existing assay
#' se2 = rpkm_normalization(se, log_trans = FALSE)
#' head(assay(se2))
#'
#' # If new_assay_name is given, normalization stored in a new assay
#' se2 = rpkm_normalization(se, log_trans = FALSE, new_assay_name =
#' 'rpkm_counts')
#' head(assay(se2, 'rpkm_counts'))
#'
#' # Creating a new assay to test specific input
#' new_matrix = matrix(data = sample(x = seq(1, 100000),
#'                                   size = nrow(se) * ncol(se),
#'                                   replace = TRUE),
#'                     nrow = nrow(se),
#'                     ncol = ncol(se))
#' rownames(new_matrix) = rownames(se)
#' colnames(new_matrix) = colnames(se)
#'
#' assay(se, 'new_counts') = new_matrix
#' se2 = rpkm_normalization(se, new_assay_name = 'rpkm_counts_new',
#' assay_name = 'new_counts')
#' head(assay(se2, 'rpkm_counts_new'))
#'
#' @export
rpkm_normalization <- function(x,
                                gene_length     = NULL,
                                log_trans       = FALSE,
                                assay_name      = NULL,
                                new_assay_name  = NULL) {

    #---------------------------
    # SummarizedExperiment path
    #---------------------------
    if (inherits(x, "SummarizedExperiment")) {
        m <- .get_matrix(se = x, a_name = assay_name)

        mat <- m$mat
        assay_name <- m$assay_name

        # Retrieve gene length from rowData
        gene_len_vec <- .get_gene_length_se(se = x, col = 'gene_length')

        ## Calculations
        # 1) Library size (col sums per million)
        lib_size <- colSums(mat, na.rm = TRUE) / 1e6
        lib_size[lib_size <= 0] <- 1

        # 2) Normalize to per million
        rpk <- sweep(mat, 2, lib_size, "/")

        # 3) Normalize per kb
        gene_kb <- gene_len_vec / 1000
        rpkm_mat <- sweep(rpk, 1, gene_kb, "/")

        # 4) Optional log2(... +1)
        if (log_trans) {
            rpkm_mat <- log2(rpkm_mat + 1)
        }

        # 5) Store result in new or existing assay

        x <- .assing_assay(se = x, matrix = rpkm_mat, a_name = assay_name,
                new_a_name = new_assay_name)

        return(x)

    #------------------------#
    # data.frame / matrix path
    #------------------------#
    } else if (is.data.frame(x) || is.matrix(x)) {
        x <- .get_matrix_df(x)

        gene_length <- .eval_gene_length_df(gene_length, df = x)

    ## Calculations
    # 1) Library size
    lib_size <- colSums(x, na.rm = TRUE) / 1e6
    lib_size[lib_size <= 0] <- 1

    # 2) Normalize to per million
    rpk <- sweep(x, 2, lib_size, "/")

    # 3) Normalize per kb
    gene_kb <- gene_length / 1000
    rpkm_mat <- sweep(rpk, 1, gene_kb, "/")

    # 4) Optional log2(... +1)
    if (log_trans) {
        rpkm_mat <- log2(rpkm_mat + 1)
    }

    # Return the RPKM (or log2-RPKM) matrix
    return(rpkm_mat)

    } else {
        stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
    }
}
