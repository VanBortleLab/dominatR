#' Dominance Plot in Two Dimensions - Plot Rope
#'
#' Creates a rope plot based on the entropy scores for each row in a dataset
#'
#' @param data A data frame with two columns, one for each variable of interest
#' @param push_text The distance from the label to the plot
#' @param rope_width Width of the rope, default is 1
#' @param rope_color Color of the rope
#' @param rope_border A logical parameter, if TRUE, the rope will have a border
#' @param col A list of two colors used for the observations that pass the threshold
#' @param col_bg The color used for the observations that do not pass the threshold
#' @param pch A list of two shapes used for the observations that pass the threshold
#' @param pch_bg The shape of the points that do not pass the threshold
#' @param cex Point size
#' @param entropyrange Range to filter observations based on their entropy score
#' @param maxvaluerange Range to filter observations based on their value
#' @param plotAll A logical parameter, if TRUE, all the points will be plot, if FALSE only points above threshold
#' @param label A logical parameter. if TRUE, the names of each variable will be written on each side of the rope
#' @param title Title for the plot
#' @param output_table A logical parameter, if TRUE the dataframe used for making the plot will be returned
#'
#' @returns
#' A rope plot and a dataframe with the information used to build the plot
#'
#' @export
#' @import dplyr
#' @seealso [entropy()], [Qentropy()], [centmass()]
#'
#'
plot_rope=function(data,
                  push_text=1,
                  rope_width=1,
                  rope_color="#CCCCCCCC",
                  rope_border=T,
                  col=c("red","blue"),
                  col_bg="whitesmoke",
                  pch = c(21,21),
                  pch_bg=19,
                  cex=1,
                  entropyrange=c(0,Inf),
                  maxvaluerange=c(0,Inf),
                  plotAll=T,
                  label=T,
                  title="title",
                  output_table=T
)
{

  rope_width = rope_width * 0.25
  # data=as.data.frame(abs(matrix(rnorm(nu*2,100,sd=100),nrow=nu,ncol=2)))
  n = ncol(data)
  px = c(-1, 1)
  py = c(0, 0)

  colu = colnames(data)
  colnames(data) = c("a", "b")

  final = cbind(data, centmass(data, px, py))
  final$comy = rope_width * sample(c(-500:500), nrow(final), replace = T) /
    500
  final = final %>%
    dplyr::mutate(color = ifelse(a > b, col[1], col[2]))

  ##Filtering on the basis of maxvalue and entropy
  final = final %>% dplyr::mutate(maxvalue = pmax(a, b), entropy = entropy(data)$Entropy)
  final$color[!(final$entropy >= entropyrange[1] &
                  final$entropy <= entropyrange[2])] = col_bg
  final$color[!(final$max >= maxvaluerange[1] &
                  final$max <= maxvaluerange[2])] = col_bg

  #plotting background
  # par(pty="s", mar = 0.5*c(5, 5, 5, 5))
  plot(
    0,
    0,
    xlim = px * push_text * 1.5,
    pch = 16,
    xlab = "",
    ylab = "",
    frame.plot = F,
    ann = T,
    axes = F,
    col = "transparent",
    main = title
  )
  polygon(
    x = c(-1, -1, 1, 1),
    y = rope_width * c(-1, 1, 1, -1),
    col = rope_color,
    border = rope_border
  )
  if (plotAll == T) {
    badpoints = final %>% dplyr::filter(color == col_bg)
    points(
      badpoints$comx,
      badpoints$comy,
      pch = pch_bg,
      cex = cex,
      col = unlist(ifelse(
        pch_bg %in% c(21:25), "black", list(badpoints$color)
      )),
      bg = unlist(ifelse(
        pch_bg %in% c(21:25), list(badpoints$color), "black"
      ))
    )
  }

  goodpoints = final %>% dplyr::filter(color != col_bg)

  #plotting highlighted points
  type1 =  goodpoints %>% filter(color == col[1])
  points(
    type1$comx,
    type1$comy,
    pch = pch[1],
    cex = cex,
    col = unlist(ifelse(
      pch[1] %in% c(21:25), "black", list(type1$color)
    )),
    bg = unlist(ifelse(
      pch[1] %in% c(21:25), list(type1$color), "black"
    ))
  )
  type2 =  goodpoints %>% filter(color == col[2])
  points(
    type2$comx,
    type2$comy,
    pch = pch[2],
    cex = cex,
    col = unlist(ifelse(
      pch[2] %in% c(21:25), "black", list(type2$color)
    )),
    bg = unlist(ifelse(
      pch[2] %in% c(21:25), list(type2$color), "black"
    ))
  )



  if (label == T)
    text(colu, x = px * push_text * 1.2, y = py)
  if (output_table == T)
    return(final)

}
