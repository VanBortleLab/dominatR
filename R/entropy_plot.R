#' Entropy - Shannon Information
#'
#' @param test A Data frame, must contain numerical columns
#' @description Calculates the Entropy value for each row across numerical columns.
#'
#'
#' @return
#' Returns an extra column to the dataframe with a value related to the Entropy value.
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr
#' @references
#' https://en.wikipedia.org/wiki/Entropy_(information_theory)
#' @examples
#' data(rnapol_chip)
#' test =  entropy(rnapol_chip)
#'
#'
#'
entropy = function(test){
  test = test |>  mutate(suma = rowSums(across(where(is.numeric))))
  test <- test %>%
    mutate(across(where(is.numeric), ~ . / suma))

  test = test |> select(-suma)

  test <- test %>%
    rowwise() %>%
    mutate(Entropy = -sum(c_across(where(is.numeric)) * ifelse(c_across(where(is.numeric)) == 0,0, log2(c_across(where(is.numeric))))))

  return(test)
}


#' Categorical Entropy
#'
#' @param test A Data frame, must contain numerical columns and Entropy should be calculated already
#'
#' @description Calculates the Categorical Entropy value related to each column. Entropy calculation should be performed before using this function
#'
#'
#'
#' @return
#' Returns the data frame with modified values for numerical columns.
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr
#' @examples
#' data(rnapol_chip)
#' test =  entropy(rnapol_chip)
#' test = Qentropy(test)
#'
#' @seealso [func(entropy)]
Qentropy = function(test){
  test <- test %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~ Entropy - ifelse(. == 0, 0, log2(.))))

  test = test |> select(-Entropy)

  return(test)
}



#' Entropy Plotter
#'
#' @param n Total number of variables to account for Entropy calculation. It is the total number of columns
#' @param data A numerical dataframe. Each column represents the variable used to calculate the entropy and build the plot
#' @param rect A logical parameter. If TRUE, points will be displayed in a single rect line
#' @param back_alpha A numerical parameter, controls the intensity of the background for the plot
#' @param label It determines if the variables are labeled around the circle or displayed as a legend. Two options are available, 'Curve' or 'Legend'
#' @param labels A vector that contains the numerical location of the variables that want to be labeled in the plot. The default is to present labels for each variable, but it can be modified if a vector with the numerical location of the variable of interest is provided
#'
#' @return
#' Returns a list of objects. It contains the domination plot, a data frame with the calculated entropy, and a dataframe with the calculated categorical entropy
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr ggforce geomtextpath
#'
#' @seealso [func(entropy), func(Qentropy)]
#'
#'
#'
circle = function(n, data, rect = F, back_alpha = 0.05, label = c('curve', 'legend'), labels = c(1:ncol(data))) {
  area = back_alpha
  a = ifelse(n >15, 80, 70)
  b = ifelse(n > 15, 95, 110)
  size = ifelse(n >15, 2,3)
  rect1 = rect
  deg = 2 * pi / n
  deg_sp = (2 * pi / n) / 2
  labels1 = colnames(data)[labels]
  location = data.frame(
    col = colnames(data),
    deg = pi / 2 - deg * (1:n),
    start = a * pi / 180 - deg_sp - deg * (1:n - 1),
    end = b * pi / 180 - deg_sp - deg * (1:n),
    curv_start = pi / 2 - deg_sp - deg * (1:n - 1) ,
    curv_end = pi / 2 - deg_sp - deg * (1:n)
  )

  data3 = data
  location$labels = NA
  location$labels[location$col %in%labels1] =   location$col[location$col %in%labels1]

  location = location |> mutate(x = 100.5*(cos(curv_start)), xend = 100.5*(cos(curv_end)), y=100.5*(sin(curv_start)), yend =100.5*(sin(curv_end)))
  arc = NULL

  for (i in 1:n) {
    p = ggplot() + geom_arc(aes(
      x0 = 0,
      y0 = 0,
      r = 100,
      start = deg_sp + deg * (i - 1),
      end = deg_sp + deg * (i)
    ))
    poly = rbind(c(0, 0),
                 data.frame(x = ggplot_build(p)$data[[1]]$x, y = ggplot_build(p)$data[[1]]$y),
                 c(0, 0))
    poly$type = colnames(data)[i]

    arc = rbind(arc, poly)
  }

  a = 100 / (n - 1)
  data1 = data.frame(Entropy = n:1, y = c(0, 1:(n - 1) * a))
  lm = lm(y ~ Entropy, data = data1)

  data = entropy(data)
  data_1 = Qentropy(data)
  data_1$col = apply(data_1, 1, function(row)
    colnames(data_1)[which.min(row)])
  data = data |> select(Entropy) |> mutate(Entropy = 2 ^ Entropy)
  data$rad = predict(lm, newdata = data)
  data = data |> bind_cols(data_1[, ncol(data_1)])
  data = data |> left_join(location, join_by(col))
  data = data |> mutate(rand_deg = sample(seq(from=start, to=end, length.out = 10), 1))

  data = data |> select(Entropy, col, rad, deg, rand_deg) |> mutate(
    x = rad * cos(ifelse(rect1 == TRUE, deg, rand_deg)),
    y = rad * sin(ifelse(rect1 == TRUE, deg, rand_deg)),
    alpha = 1 - Entropy / n
  )

  data1 = data.frame(rad = a * (1:(n - 1)))

  if(label == 'legend'){
    circle = ggplot() +  theme_minimal() +
      geom_polygon(
        data = arc,
        aes(x, y, fill = type),
        color = 'gray90',
        show.legend = F,
        alpha = area
      ) +
      geom_circle(data = data1,
                  aes(x0 = 0, y0 = 0, r = rad),
                  col = 'gray80')  +
      theme(
        aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.title = element_blank()
      ) + geom_jitter(
        data = data,
        aes(x, y, col = col, alpha = alpha),
        pch = 19,
        size = ifelse(n >15, 0.8, 1.5),
        width = 3,
        height = 3
      ) +
      scale_alpha(guide = 'none') +
      coord_equal(xlim = c(-105,105), ylim = c(-105,105))
  }
  else{
    circle = ggplot() +  theme_minimal() +
      geom_labelcurve(data = location, aes(x = x, xend = xend, y = y, yend = yend, label = labels), linecolour = 'white', curvature = -0.5, size = size, na.rm = T) +
      geom_polygon(
        data = arc,
        aes(x, y, fill = type),
        color = 'gray90',
        show.legend = F,
        alpha = area
      ) +
      geom_circle(data = data1,
                  aes(x0 = 0, y0 = 0, r = rad),
                  col = 'gray80')  +
      theme(
        aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.title = element_blank()
      ) + geom_jitter(
        data = data,
        aes(x, y, col = col, alpha = alpha),
        pch = 19,
        size = ifelse(n >15, 0.8, 1.5),
        width = 3,
        height = 3, show.legend = F
      ) +
      scale_alpha(guide = 'none') +
      coord_equal(xlim = c(-105,105), ylim = c(-105,105))
  }

  return(circle) }

