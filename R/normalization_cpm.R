

#' Counts Per Million
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#'
#' @description Normalizes a count matrix by the counts per million method
#' @return A Matrix with normalized counts.
#' @export
#'
cpm_normalization = function(df, log_trans = FALSE) {
  sum = colSums(df) / 1e6

  df = mapply('/', df, sum)

  if (log_trans == FALSE) {
    return(df)
  } else {
    df = log2(df + 1)
    return(df)
  }
}
