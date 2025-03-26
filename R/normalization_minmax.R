#' MinMax Normalization
#'
#' @param x data that needed to be normalized
#' @param min the minimum value of the final normalized data
#' @param max the maximum value of the final normalized data
#' @return a vector with normalized data
#' @export
#'
#'
minmax_normalization=function(x,min,max)
{
  MIN=min(x, na.rm = TRUE)
  MAX=max(x, na.rm = TRUE)
  y=as.double((x-MIN)*(max-min)/(MAX-MIN)+min)
  return(y)
}

