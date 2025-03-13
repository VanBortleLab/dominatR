
#' Center of Mass
#' @description Returns coordinates of center of mass from the data provided.
#'
#' @param data a table that provides mass of different objects
#' @param x a vector that provides x coordinates of  objects
#' @param y a vector that provides y coordinates of  objects
#' @return
#' Returns a table with 2 columns(i.e. comx and comy) containing x coordinate and y coordinate of the center of mass
#' @export
#'
#'
centmass=function(data,x=c(0,1,0.5) ,y=c(0,0,sqrt(3)/2)){

  sum = rowSums(data)
  data = data.frame(comx=rowSums(as.matrix(data)%*%diag(x)),
                    comy=rowSums(as.matrix(data)%*%diag(y)))

  data = data/sum

return(data)}
