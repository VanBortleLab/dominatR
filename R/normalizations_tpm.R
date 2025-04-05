#' TPM Normalization
#'
#' @description
#' Normalizes read counts by the TPM (Transcripts Per Million) method:
#'   \enumerate{
#'     \item If \code{log_trans = TRUE}, applies \code{log2(TPM + 1)}.
#'   }
#'
#' @details
#' If \code{x} is a \code{SummarizedExperiment}, this function looks for a
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
#'   only if \code{x} is a data.frame or matrix. Must match the number of rows in
#'   \code{x}. Ignored if \code{x} is a SummarizedExperiment.
#' @param log_trans Logical. If \code{TRUE}, apply \code{log2(... + 1)} transform
#'   to the TPM-normalized values.
#' @param assay_name If \code{x} is a SummarizedExperiment, name of the assay to
#'   normalize. Defaults to the first assay if not specified.
#' @param new_assay_name If \code{x} is a SummarizedExperiment, name of a new
#'   assay in which to store the TPM (or log2-TPM). If \code{NULL}, overwrites
#'   the assay specified in \code{assay_name}.
#'
#' @return A numeric \strong{matrix} of TPM or log2(TPM + 1) values if
#'   \code{x} is a data.frame or matrix. If \code{x} is a SummarizedExperiment,
#'   returns the modified SummarizedExperiment with the TPM data placed in the
#'   existing or new assay.
#'
#' @examples
#'library(SummarizedExperiment)
#'library(airway)
#'data('airway')
#'
#'se = airway
#'
#'### Adding a column in rowData regarding the gene_length
#'rowData(se)$gene_length = rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
#'
#'# -------------------------------
#'# 1) Using a data.frame
#'# -------------------------------
#'
#'gene_length = rowData(se)$gene_length
#'
#'df = assay(se)
#'
#'## Without log transformation
#'df = tpm_normalization(df, gene_length = gene_length)
#'
#'df[1:5, 1:5]
#'
#'## With log transformation
#'df = tpm_normalization(df, gene_length = gene_length, log_trans = TRUE)
#'
#'df[1:5, 1:5]
#'
#'# -------------------------------
#'# 2) Using a SummarizedExperiment
#'# -------------------------------
#'
#'# If now new_assay_name is provided, then overwrites existing assay
#'se2 = tpm_normalization(se, log_trans = FALSE)
#'
#'head(assay(se2))
#'
#'# If new new_assay_name, normalization stored in a new object
#'se2 = tpm_normalization(se, log_trans = FALSE, new_assay_name = 'tpm_counts')
#'
#'head(assay(se2, 'tpm_counts'))
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
#'se2 = tpm_normalization(se, new_assay_name = 'tpm_counts_new', assay_name = 'new_counts')
#'
#'se2
#'
#'head(assay(se2, 'tpm_counts_new'))
#' @export

tpm_normalization <- function(x,
                               gene_length     = NULL,
                               log_trans       = FALSE,
                               assay_name      = NULL,
                               new_assay_name  = NULL) {

  #---------------------------
  # SummarizedExperiment path
  #---------------------------
  if (inherits(x, "SummarizedExperiment")) {
    if (is.null(assay_name)) {
      all_assays <- SummarizedExperiment::assayNames(x)
      if (length(all_assays) < 1) {
        stop("No assays found in the SummarizedExperiment.")
      }
      assay_name <- all_assays[[1]]
    }
    mat <- SummarizedExperiment::assay(x, assay_name)
    if (is.null(mat)) {
      stop("No assay named '", assay_name, "' found in the SummarizedExperiment.")
    }
    if (!is.numeric(mat)) {
      stop("Selected assay is not numeric. Please provide numeric data for TPM normalization.")
    }

    # Retrieve gene lengths from rowData
    rd <- rowData(x)
    if (!("gene_length" %in% colnames(rd))) {
      stop("No 'gene_length' column found in rowData(x). ",
           "Please add rowData(x)$gene_length or use data.frame/matrix mode.")
    }
    gene_len_vec <- rd[["gene_length"]]
    if (!is.numeric(gene_len_vec)) {
      stop("'gene_length' in rowData(x) must be numeric.")
    }

    ## Calculations
    # 1) RPK
    rpk <- mat

    rpk <- sweep(mat, 1, gene_len_vec, "/") * 1000

    # 2) Column sums
    col_scaling <- colSums(rpk, na.rm=TRUE)
    # avoid zero or negative
    col_scaling[col_scaling <= 0] <- 1

    # 3) TPM
    tpm_mat <- sweep(rpk, 2, col_scaling, '/') * 1e6

    # 4) Optional log2(... +1)
    if (log_trans) {
      tpm_mat <- log2(tpm_mat + 1)
    }

    # 5) Store result in new or existing assay
    if (is.null(new_assay_name)) {
      assay(x, assay_name) <- tpm_mat
    } else {
      assay(x, new_assay_name) <- tpm_mat
    }

    return(x)

    #------------------------#
    # data.frame / matrix path
    #------------------------#
  } else if (is.data.frame(x) || is.matrix(x)) {
    # Convert data.frame to matrix if needed
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    if (!is.numeric(x)) {
      stop("Input 'x' must contain numeric data for TPM normalization.")
    }
    # Check user-provided gene_length
    if (is.null(gene_length)) {
      stop("You must provide 'gene_length' for data.frame/matrix input.")
    }
    if (!is.numeric(gene_length)) {
      stop("Argument 'gene_length' must be numeric.")
    }
    if (length(gene_length) != nrow(x)) {
      stop("Length of 'gene_length' (", length(gene_length),
           ") must match the number of rows (", nrow(x), ") in 'x'.")
    }


    ## Calculations
    # 1) RPK
    rpk <- x

    rpk <- sweep(rpk, 1, gene_length, "/") * 1000

    # 2) Column sums
    col_scaling <- colSums(rpk, na.rm=TRUE)
    # avoid zero or negative
    col_scaling[col_scaling <= 0] <- 1

    # 3) TPM
    tpm_mat <- sweep(rpk, 2, col_scaling, '/') * 1e6

    # 4) Optional log2(... +1)
    if (log_trans) {
      tpm_mat <- log2(tpm_mat + 1)
    }

    # Return the TPM (or log2-TPM) matrix
    return(tpm_mat)

  } else {
    stop("Input must be a matrix/data.frame or a SummarizedExperiment.")
  }
}
