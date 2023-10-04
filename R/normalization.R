#' Quantile Normalization
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#' @description Normalizes a count matrix by the quantile normalization method
#' @return A Matrix with quantile normalized values
#'
#'
quant.normalization <- function(df, log_trans = F){
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


#' Counts Per Million
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#'
#' @description Normalizes a count matrix by the counts per million method
#' @return A Matrix with normalized counts.
#'
#'
cpm = function(df, log_trans = F){
  sum = colSums(df)/1e6


  df = mapply('/', df, sum)

  if(log_trans == F){
  return(df)
  } else {
  df = log2(df + 1)
  return(df)
  }
}


#' Transcripts per Million
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param gene_length  A vector with the length for the genes in the dataframe in basepair format
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#'
#' @description Normalizes a count matrix by the TPM method
#' @return A Matrix with normalized counts.
#'

tpm = function(df, gene_length, log_trans = F){

  ### Divide counts by Gene Length
  df1 <- 1000 * df / gene_length

  ### Sequencing Depth
  counts <- colSums(df1)/1e6

  ### Divide previous by Sequencing depth
  df1 <- mapply('/', df1, counts)

  if(log_trans == F){
    return(df1)
  } else {
    df1 = log2(df1 + 1)
    return(df1)
  }

}

#' Reads Per Kilobase Million
#'
#' @param df A dataframe or matrix, must contain only numerical valeus
#' @param gene_length  A vector with the length for the genes in the dataframe in basepair format
#' @param log_trans Logical attribute, determines if normalized data should be log2 transformed. Note: An extra count is added.
#'
#' @description Normalizes a count matrix by the RPKM method
#' @return A Matrix with normalized counts.
#'
rpkm = function(df, gene_length, log_trans = F){

  sum = colSums(df)/1e6

  df = mapply('/', df, sum)

  df = df * 1000 / gene_length

  if(log_trans == F){
    return(df)
  } else {
    df = log2(df + 1)
    return(df)
  }

}
