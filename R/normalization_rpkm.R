#' Reads Per Kilobase Million
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param gene_length  A vector with the length for the genes in the dataframe in basepair format
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#'
#' @description Normalizes a count matrix by the RPKM method
#' @return A Matrix with normalized counts.
#' @export
rpkm_normalization = function(df, gene_length, log_trans = FALSE) {
  sum = colSums(df) / 1e6

  df = mapply('/', df, sum)

  df = df * 1000 / gene_length

  if (log_trans == FALSE) {
    return(df)
  } else {
    df = log2(df + 1)
    return(df)
  }

}
