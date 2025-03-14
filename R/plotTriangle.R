#' Plotting a dummy triangle
#'
#' @export
#' @import dplyr
#'
plot_Triangle_dummy=function()
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
