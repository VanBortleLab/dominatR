#' Entropy - Shannon Information
#'
#' @param test A Data frame, must contain numerical columns
#' @description Calculates the Entropy value for each row across numerical columns.
#'
#'
#' @return
#' Returns an extra column to the dataframe with a value related to the Entropy value. :)
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr
#' @references
#' https://en.wikipedia.org/wiki/Entropy_(information_theory)
#' @examples
#' data(rnapol_score)
#' test =  entropy(rnapol_score)
#'
#'
#'
entropy = function(test){
  test = test |>  mutate(suma = rowSums(across(where(is.numeric))))
  test = test |> mutate(suma = ifelse(suma > 0 , suma, 1))
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
#' data(rnapol_score)
#' test =  entropy(rnapol_score)
#' test = Qentropy(test)
#'
#' @seealso [entropy()]
Qentropy = function(test){
  test <- test %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~ ifelse(. == 0, Inf, Entropy - log2(.))))

  test = test |> select(-Entropy)

  return(test)
}



#' Entropy Plotter
#'
#' Creates a circular plot based on the Entropy for each row in a dataset.
#'
#' @param n Total number of variables to account for Entropy calculation. It is the total number of columns
#' @param data A numerical dataframe. Each column represents the variable used to calculate the entropy and build the plot
#' @param str A logical parameter. If TRUE, points will be displayed in a single rect line
#' @param back_alpha A numerical parameter, controls the intensity of the background for the plot
#' @param label It determines if the variables are labeled around the circle or displayed as a legend. Two options are available, 'Curve' or 'Legend'
#' @param variables A vector that contains the names of the columns that should appear in the plot when label == 'Curve'.
#'@param title Title of the plot.
#'@param threshold A numerical value that determines the minimum value for rows to be filtered. Any row with a rowSums value less than this threshold is removed from the subsequent analysis. Default is 0
#'@param col_variable The name of a variable with categorical data that wants to be displayed. The function in general adds color to the observations based on their column names, if each observation has a categorical value linked to it, by providing the name of that variable, data points will be colored accordingly.
#'@param point_size Size for the plot points, default is 3.
#'@param text_size Size for the text. 'It works when label is around the circle.
#'@param line_col Color for the plot outline.
#'@param out_line Color for the outtermost circle
#' @return
#' Returns a list of objects. It contains the domination plot, an a dataframe with the plotting coordinates with entropy and dominant variable for each point
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr ggforce geomtextpath
#'
#' @seealso [entropy()], [Qentropy()]
#'
#'
#'
plot_circle = function(n, data, str = F,
                       back_alpha = 0.05, label = c('curve', 'legend'), variables = colnames(data),
                       title = NULL, threshold = 0, col_variable = NULL,
                       point_size = 3, text_size = 3,  line_col = 'gray90', out_line = 'black'){
  if(length(col_variable > 0)){
    colnames(data)[which(colnames(data) == col_variable)] = 'Factor'
    colnames(data) = substr(colnames(data), 1, 10)
    variables = substr(variables,1,10)

    ## removing rows where sum = 0
    data = data |> relocate(Factor)
    data = data |> dplyr::filter(rowSums(select_if(data, is.numeric)) > threshold)

    area = back_alpha
    a = ifelse(n >15, 80, 70)
    b = ifelse(n > 15, 100, 110)
    sizex = point_size
    textsize = text_size
    rect1 = str
    deg = 2 * pi / n
    deg_sp = (2 * pi / n) / 2
    labels1 = variables
    location = data.frame(
      col = colnames(select_if(data, is.numeric)),
      deg = pi / 2 - deg * (1:n),
      start = a * pi / 180 - deg_sp - deg * (1:n - 1),
      end = b * pi / 180 - deg_sp - deg * (1:n),
      curv_start = pi / 2 - deg_sp - deg * (1:n - 1) ,
      curv_end = pi / 2 - deg_sp - deg * (1:n)
    )

    location$labels = NA
    location$labels[location$col %in%labels1] =   location$col[location$col %in%labels1]

    rad_label = ifelse(n == 3, 120, 100.5)

    location = location |> mutate(x = rad_label*(cos(curv_start)), xend = rad_label*(cos(curv_end)), y=rad_label*(sin(curv_start)), yend =rad_label*(sin(curv_end)))
    arc = NULL


    numeric_cols <- colnames(data)[sapply(data, is.numeric)]

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
      poly$type = numeric_cols[i]

      arc = rbind(arc, poly)
    }

    arc$type = factor(arc$type, levels = numeric_cols)

    a = 100 / (n - 1)
    data1 = data.frame(Entropy = n:1, y = c(0, 1:(n - 1) * a))
    lm = lm(y ~ Entropy, data = data1)

    data = entropy(data)
    data_1 = Qentropy(data)

    numeric_cols1 = sapply(data_1, is.numeric)

    data_1$col = suppressWarnings(
      apply(data_1[, numeric_cols1], 1, function(row) {

        min_cols = colnames(data_1)[numeric_cols1][row == min(row)]

        if(length(min_cols) > 1) {
          min_cols = sample(min_cols, 1)
        }
        return(min_cols)
      })
    )

    data = data |> select(Entropy, Factor) |> mutate(Entropy = 2 ^ Entropy)
    data$rad = predict(lm, newdata = data)
    data = data |> bind_cols(data_1[, ncol(data_1)])
    data = data |> left_join(location, join_by(col))
    data = data |> mutate(rand_deg = sample(seq(from=start, to=end, length.out = 10), 1))

    data = data |> select(Entropy, Factor, col, rad, deg, rand_deg) |> mutate(
      x = rad * cos(ifelse(rect1 == TRUE, deg, rand_deg)),
      y = rad * sin(ifelse(rect1 == TRUE, deg, rand_deg)),
      alpha = 1 - Entropy / n
    )

    data$col = factor(data$col, levels = numeric_cols)

    data1 = data.frame(rad = a * (1:(n - 1)))

    if(label == 'legend'){
      circle = ggplot2::ggplot() +  theme_minimal() +
        geom_point(data = data,
                   aes(x, y, fill = col),
                   pch = 21,
                   col = 'white',
                   size = ifelse(n >15, 1.2, 2))+
        geom_polygon(
          data = arc,
          aes(x, y), color = 'white',fill = 'white') +
        geom_polygon(
          data = arc,
          aes(x, y, fill = type),
          color = line_col,
          show.legend = F,
          alpha = area
        ) +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad),
                    col = line_col)  +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad[nrow(data1)]),
                    col = out_line)  +
        theme(
          aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor  = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.background = element_blank()
        ) + geom_jitter(
          data = data,
          aes(x, y, fill2 = Factor, alpha = alpha, col2 = Factor ),
          pch = 21,
          size = sizex,
          width = 3,
          height = 3
        ) |> rename_geom_aes(new_aes = c('fill' = 'fill2', 'col' = 'col2'))  +
        scale_alpha(guide = 'none') +
        #guides(fill = 'none') +
        coord_equal(xlim = c(-105,105), ylim = c(-105,105)) + ggtitle(title) +
        scale_fill_manual(aesthetics = 'fill2', values = scales::hue_pal()(n)) +
        scale_color_manual(aesthetics = 'col2', values = rep('black', n))
    }
    else{
      circle = ggplot() +  theme_minimal() +
        geom_labelcurve(data = location, aes(x = x, xend = xend, y = y, yend = yend, label = labels), linecolour = 'white', curvature = -0.5, size = textsize, na.rm = T, show.legend = F) +
        geom_polygon(
          data = arc,
          aes(x, y, fill = type),
          color = line_col,
          show.legend = F,
          alpha = area
        ) +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad),
                    col = line_col)  +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad[nrow(data1)]),
                    col = out_line) +
        theme(
          aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor  = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()
        ) + geom_jitter(
          data = data,
          aes(x, y, fill2 = Factor, alpha = alpha, col2 = Factor),
          pch = 21,
          size = sizex,
          width = 3,
          height = 3, show.legend = T
        ) |> rename_geom_aes(new_aes = c('fill' = 'fill2', 'col' = 'col2')) +
        scale_alpha(guide = 'none') +
        guides(fill = 'none') +
        coord_equal(xlim = c(-105,105), ylim = c(-105,105)) + ggtitle(title) +
        scale_fill_manual(aesthetics = 'fill2', values = scales::hue_pal()(n)) +
        scale_color_manual(aesthetics = 'col2', values = rep('black', n))
    }
    data = data |> select(Entropy, Factor, col, x, y) |> mutate(Entropy = log2(Entropy))

  } else {
    colnames(data) = substr(colnames(data), 1, 10)
    variables = substr(variables,1,10)

    ## removing rows where sum = 0
    data = data |> dplyr::filter(rowSums(select_if(data, is.numeric)) > threshold)

    area = back_alpha
    a = ifelse(n >15, 80, 70)
    b = ifelse(n > 15, 100, 110)
    sizex = point_size
    textsize = text_size
    rect1 = str
    deg = 2 * pi / n
    deg_sp = (2 * pi / n) / 2
    labels1 = variables
    location = data.frame(
      col = colnames(select_if(data, is.numeric)),
      deg = pi / 2 - deg * (1:n),
      start = a * pi / 180 - deg_sp - deg * (1:n - 1),
      end = b * pi / 180 - deg_sp - deg * (1:n),
      curv_start = pi / 2 - deg_sp - deg * (1:n - 1) ,
      curv_end = pi / 2 - deg_sp - deg * (1:n)
    )

    location$labels = NA
    location$labels[location$col %in%labels1] =   location$col[location$col %in%labels1]

    rad_label = ifelse(n == 3, 120, 100.5)

    location = location |> mutate(x = rad_label*(cos(curv_start)), xend = rad_label*(cos(curv_end)), y=rad_label*(sin(curv_start)), yend =rad_label*(sin(curv_end)))
    arc = NULL

    numeric_cols <- colnames(data)[sapply(data, is.numeric)]

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
      poly$type = numeric_cols[i]

      arc = rbind(arc, poly)
    }

    arc$type = factor(arc$type, levels = numeric_cols)

    a = 100 / (n - 1)
    data1 = data.frame(Entropy = n:1, y = c(0, 1:(n - 1) * a))
    lm = lm(y ~ Entropy, data = data1)

    data = entropy(data)
    data_1 = Qentropy(data)
    numeric_cols1 = sapply(data_1, is.numeric)

    data_1$col = suppressWarnings(
      apply(data_1[, numeric_cols1], 1, function(row) {

        min_cols = colnames(data_1)[numeric_cols1][row == min(row)]

        if(length(min_cols) > 1) {
          min_cols = sample(min_cols, 1)
        }
        return(min_cols)
      })
    )
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
    data$col = factor(data$col, levels = numeric_cols)

    data1 = data.frame(rad = a * (1:(n - 1)))

    if(label == 'legend'){
      circle = ggplot2::ggplot() +  theme_minimal() +
        geom_polygon(
          data = arc,
          aes(x, y, fill = type),
          color = line_col,
          show.legend = F,
          alpha = area
        ) +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad),
                    col = line_col)  +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad[nrow(data1)]),
                    col = out_line) +
        theme(
          aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor  = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()
        ) + geom_jitter(
          data = data,
          aes(x, y, fill2 = col, alpha = alpha, col2 = col),
          pch = 21,
          size = sizex,
          width = 3,
          height = 3
        ) |> rename_geom_aes(new_aes = c('fill' = 'fill2', 'col' = 'col2')) +
        scale_alpha(guide = 'none') +
        coord_equal(xlim = c(-105,105), ylim = c(-105,105)) + ggtitle(title) +
        scale_fill_manual(aesthetics = 'fill2', values = scales::hue_pal()(n)) +
        scale_color_manual(aesthetics = 'col2', values = rep('black', n))
    }
    else{
      circle = ggplot() +  theme_minimal() +
        geom_labelcurve(data = location, aes(x = x, xend = xend, y = y, yend = yend, label = labels), linecolour = 'white', curvature = -0.5, size = textsize, na.rm = T) +
        geom_polygon(
          data = arc,
          aes(x, y, fill = type),
          color = line_col,
          show.legend = F,
          alpha = area
        ) +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad),
                    col = line_col)  +
        geom_circle(data = data1,
                    aes(x0 = 0, y0 = 0, r = rad[nrow(data1)]),
                    col = out_line) +
        theme(
          aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor  = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()
        ) + geom_jitter(
          data = data,
          aes(x, y, fill2 = col, alpha = alpha, col2 = col),
          pch = 21,
          size = sizex,
          width = 1,
          height = 1, show.legend = F
        ) |> rename_geom_aes(new_aes = c('fill' = 'fill2', 'col' = 'col2')) +
        scale_alpha(guide = 'none') +
        coord_equal(xlim = c(-105,105), ylim = c(-105,105)) + ggtitle(title) +
        scale_fill_manual(aesthetics = 'fill2', values = scales::hue_pal()(n)) +
        scale_color_manual(aesthetics = 'col2', values = rep('black', n))
    }
    data = data |> select(Entropy, col, x, y) |> mutate(Entropy = log2(Entropy))}
  return(list(circle, data, data_1 ))
  }

#' Dominance Distribution Plotter
#'
#' This function depends on the plot.circle function when the parameter col_variable is used. It plots the dominance frequency distribution of each of the categories present in col_variable.
#'
#' @param n Number of Variables, must be the same as the ones used in the function
#' @param circle An object type list obtained from the plot.circle function
#'
#' @return
#' Returns a list of objects, a frequency distribution plot for each category and the respective dataframe
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr ggforce geomtextpath
#'
#'
#'
plot_circle_frequency = function(n, circle){
  data_ent = circle[[2]]
  breaks_log2 = log2(1:n)
  data_ent$bin = cut(data_ent$Entropy, breaks = breaks_log2,
                     labels = paste0(floor(2^(round(breaks_log2[-length(breaks_log2)], 2)))),
                     include.lowest = T)

  data_ent = data_ent |> group_by(bin, Factor) |> dplyr::select(bin, Factor) |> summarise(n = n(), .groups = 'drop')
  data_ent = data_ent |> group_by(Factor) |> mutate(proportion = n/sum(n))
  data_ent$Factor = factor(data_ent$Factor)

  plot_stat = data_ent |> ggplot() + geom_line(aes(x = bin, y = proportion, col = Factor, group = Factor), linewidth = 1) +
    theme_minimal() + theme(
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.text.x = element_text(
      size = 10 ,
      angle = 45,
      vjust = 0.5
    ),
    legend.position = 'none',
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face =
                                  'bold'),
    legend.key.size = unit(15, 'points'),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    )
  ) +
    ylab('Proportion of genes') +
    xlab('Dominance') + facet_wrap(~Factor, ncol = 1)

  return(list(plot_stat, data_ent))
}
