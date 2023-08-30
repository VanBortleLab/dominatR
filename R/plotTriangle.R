#' Vrtx - creates regular polygon
#' @description Returns coordinates of a regular polygon with n vertices. It plots the polygon with the flag "plot=T"
#'
#' @param n an integer that indicates the number of vertex of the polygon
#' @param plot a boolean T/F that indicates wheteher to plot the polygon or not
#'
#' @return t
#' Returns the coordinates of the polygon
#' @export
#'
vrtx=function(n,plot=F)
{
  #n=8
  t=data.frame(x=sapply(0:(n-1),  function(x) sin(pi*(360/n)*x/180)),
               y=sapply(0:(n-1),  function(x) cos(pi*(360/n)*x/180)))
  if(plot==T){ plot(c(1,-1),c(1,-1),col="transparent")
    polygon(t)}
  return(t)
}

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
getCOM=function(data,x=c(0,1,0.5) ,y=c(0,0,sqrt(3)/2))
{return(data.frame(comx=rowSums(as.matrix(data)%*%diag(x)),comy=rowSums(as.matrix(data)%*%diag(y)))/rowSums(data))}






#' Plotting Triangle with center of mass
#'
#' @param data a table with three columns
#' @param entropyrange a vector with 2 elements in the range [0,1] i.e. c(minimum entropy , maximum entropy)
#' only the entries that has their entropy within this range are highlighted
#' @param magnituderange a vector with 2 elements in the range [0,1] i.e. c(minimum magnitude , maximum magnitude)
#' only the entries that has their magnitude within this range are highlighted
#' @param col a vector of 3 elements that contains color names
#' @param output_table a boolean value of T or F indicating whether to return the output table or not
#' @param plotAll a boolean value
#' If T, all the points that are out of entropyrange or magnitude range are plotted in whitesmoke color, else they are not plotted at all.
#' @param cex a numeric value indicating the size of the points
#' @param pch a integer value indicating different shapes of the points
#'
#' @return
#'  Returns a table with the coordinates of center of mass entropy(in the range of 0 and 1) and magnitude(in the range of 0 and 1)
#' @export
#' @import dplyr
#'
#'
plotTriangle=function(data,entropyrange=c(0,1),magnituderange=c(0,1),col=c("red","green","blue"),output_table=T,plotAll=T,cex=1,pch=16)
{
  # require(dplyr)
  n=3
  #data=log2(data)
  colnames(data)=c("a","b","c")
  final=cbind(data,getCOM(data,vrtx(n)$x,vrtx(n)$y))

  #Finding domination score
  final=final %>% mutate(entropy=MinMaxNorm2(ent(data),0,1), magnitude=MinMaxNorm2(a^2+b^2+c^2,0,1))
  final=final %>% mutate(max=pmax(a,b,c)) %>% mutate(color=ifelse(a==max,col[1],ifelse(b==max,col[2],col[3])))
  final$color[!(final$entropy>entropyrange[1] & final$entropy<entropyrange[2] & final$magnitude>magnituderange[1] & final$magnitude<magnituderange[2])]="whitesmoke"

  #plotting
  {
    plot(c(1,-1),c(1,-1),pch=16,xlab = "",ylab="",frame.plot = F,ann = F,axes = F,col="transparent")
    if(plotAll==T){
      badpoints=final %>% filter(color=="whitesmoke")
      points(badpoints$comx,badpoints$comy,col=badpoints$color,pch=pch,cex=cex)}

    goodpoints=final %>% filter(color!="whitesmoke")
    points(goodpoints$comx,goodpoints$comy,col=goodpoints$color,pch=pch,cex=cex)
    polygon(vrtx(n)$x,vrtx(n)$y,pch=16)
  }##Luka##

  if(output_table==T) return(final)
}


#' MinMax Normalization
#'
#' @param x data that needed to be normalized
#' @param min the minimum value of the final normalized data
#' @param max the maximum value of the final normalized data
#' @return a vector with normalized data
#' @export
#'
#'
MinMaxNorm2=function(x,min,max)
{
  MIN=min(x)
  MAX=max(x)
  y=as.double((x-MIN)*(max-min)/(MAX-MIN)+min)
  return(y)
}



# finding entropy
ent=function(data)
{
  data=data/rowSums(data)
  eata=log2(data)
  fata=eata %>% mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))
  return(-rowSums(data*fata))
}



#' Plotting a dummy triangle
#'
#' @export
#' @import dplyr
#'
plotTriangle_dummy=function()
{
  par(mar=c(1,0.5,1,0.5))
  #require(dplyr)
  repeat{
    ok=as.integer(readline(prompt=paste0("Enter the number of random points you want to plot or enter 0 to end = ","\n")))
    if(ok==0) {dev.off();break}
    if(ok>50000){
      ok2=readline(prompt=paste0('Plotting more than 50 thousands points may hang your system. Enter "YES" if you wanna proceed = ',"\n"))
      if(ok2!="YES") next}

    data=data.frame(a=runif(ok,0,1),b=runif(ok,0,1),c=runif(ok,0,1))
    data=data %>% mutate(color=rgb(a,b,c,maxColorValue = 1))
    data=cbind(data,getCOM(data[,c(1,2,3)]))
    data=data %>% arrange(a^2+b^2+c^2)

    #dev.off()

    plot(c(0,1,0.5),c(0,0,sqrt(3)/2),pch=16,xlab = "",ylab="",frame.plot = F,ann = F,axes = F,col="transparent")
    points(data$comx,data$comy,col=data$color,pch=16,cex=0.5)
    polygon(c(0,1,0.5),c(0,0,sqrt(3)/2),pch=16)
    text(0.2,0.85,"Plotting",cex=0.5)
    text(0.2,0.8,ok,cex=0.7)
    text(0.2,0.75,"random points",cex=0.5)
  }
}
