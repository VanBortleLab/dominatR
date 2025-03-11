#' Quantile Normalization
#'
#' @param df A dataframe or matrix, must contain only numerical values
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#' @description Normalizes a count matrix by the quantile normalization method
#' @return A Matrix with quantile normalized values
#' @export
#'
quant_normalization <- function(df, log_trans = F){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)

  if(log_trans == F){
    return(df_final)
  } else {
    df_final = log2(df_final + 1)
    return(df_final)
  }
}
