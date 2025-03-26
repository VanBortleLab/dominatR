
#' Transcripts per Million
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param gene_length  A vector with the length for the genes in the dataframe in basepair format
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#'
#' @description Normalizes a count matrix by the TPM method
#' @return A Matrix with normalized counts.
#' @export
tpm_normalization = function(df, gene_length, log_trans = FALSE){

  ### Divide counts by Gene Length
  df1 <- 1000 * df / gene_length

  ### Sequencing Depth
  counts <- colSums(df1)/1e6

  ### Divide previous by Sequencing depth
  df1 <- mapply('/', df1, counts)

  if(log_trans == FALSE){
    return(df1)
  } else {
    df1 = log2(df1 + 1)
    return(df1)
  }

}
