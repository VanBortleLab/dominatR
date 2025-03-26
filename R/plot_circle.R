#' Dominance Plot in any Dimension - Plot Circle
#'
#' Creates a circular plot based on the entropy scores for each row in a dataset
#'
#' @param n Total number of variables to account for entropy calculation.
#' @param data A numerical dataframe. Each column represent the variable used to calculate the entropy and build the plot
#' @param str  A logical parameter. If TRUE points will be displayed in a single rect line
#' @param back_alpha A numerical parameter, controls the intensity of the background for the plot
#' @param label It determies if the variables are labeled around the circle or displayed as a legend. Two options are available 'curve' or 'legend'
#' @param variables A vector that contains the name of the columns that should appear in the plot when label == 'curve'
#' @param title Title of the plot
#' @param threshold A numerical value that determines the minmum values for rows to be filtered. Any row with a rowSums value less than this threshold is removed from the subsequent analysis. Default is 0
#' @param col_variable The name of a variable with categorical data that wants to be displayed. The function in general adds color to the observations based on their column names, if each observation has a category linked to it, by providing the name of that variable, data points will be colored accordingly
#' @param point_size Size for the plot points, default is 3
#' @param text_size Size for the text. It works when label is 'curve'
#' @param line_col Color for the plot outline.
#' @param out_line Color for the outtermost circle
#' @param background_colors Colors for the background of the plot
#' @param point_fill_colors Colors for the points in the plot
#' @param point_line_colors Colors for the outline of the point in the plot
#'
#' @returns
#' Returns a list of objects. It contains the circular dominance plot and a dataframe with the plotting coordinates, entropy and the dominant variable for each row
#' @export
#' @import dplyr forcats lubridate purrr readr stringr tibble tidyr ggforce geomtextpath ggnewscale ggplot2
#'
#' @seealso [entropy()], [Qentropy()]
#'
#'
plot_circle = function(n,
                       data,
                       str = FALSE,
                       back_alpha = 0.05,
                       label = c('curve', 'legend'),
                       variables = colnames(data),
                       title = NULL,
                       threshold = 0,
                       col_variable = NULL,
                       point_size = 3,
                       text_size = 3,
                       line_col = 'gray90',
                       out_line = 'black',
                       background_colors = NULL,
                       point_fill_colors = NULL,
                       point_line_colors = NULL){

  if(!is.null(col_variable) && length(col_variable) > 0){

    colnames(data)[which(colnames(data) == col_variable)] = 'Factor'
    ## removing rows where sum = 0
    data = data |> relocate(Factor)
  }

  ## Shortening column names
  colnames(data) = substr(colnames(data), 1, 10)
  variables = substr(variables,1,10)

  ## Subsetting for observations below the threshold
  data = data |>
    dplyr::filter(rowSums(select_if(data, is.numeric)) > threshold)

  ## Key Plot setting
  area = back_alpha
  a = ifelse(n >15, 80, 70)
  b = ifelse(n > 15, 100, 110)
  deg = 2 * pi / n
  deg_sp = (2 * pi / n) / 2
  labels1 = variables

  ## Plot parameters
  sizex = point_size
  textsize = text_size
  rect1 = str

  ## Location data
  location = data.frame(
    col = colnames(select_if(data, is.numeric)),
    deg = pi / 2 - deg * (1:n),
    start = a * pi / 180 - deg_sp - deg * (1:n - 1),
    end = b * pi / 180 - deg_sp - deg * (1:n),
    curv_start = pi / 2 - deg_sp - deg * (1:n - 1) ,
    curv_end = pi / 2 - deg_sp - deg * (1:n)
  )

  ## For those specific variables we want to label
  location$labels = NA
  location$labels[location$col %in%labels1] =   location$col[location$col %in%labels1]

  rad_label = ifelse(n == 3, 120, 100.5)

  location = location |>
    mutate(
      x = rad_label * (cos(curv_start)),
      xend = rad_label * (cos(curv_end)),
      y = rad_label * (sin(curv_start)),
      yend = rad_label * (sin(curv_end))
    )


  ## Generate arc polygons with the information
  numeric_cols <- colnames(data)[sapply(data, is.numeric)]

  arc = NULL
  for (i in 1:n) {
    p = ggplot() +
      geom_arc(aes(
        x0 = 0,
        y0 = 0,
        r = 100,
        start = deg_sp + deg * (i - 1),
        end = deg_sp + deg * (i)
      ))
    poly = rbind(c(0, 0),
                 data.frame(x = ggplot_build(p)$data[[1]]$x,
                            y = ggplot_build(p)$data[[1]]$y),
                 c(0, 0))
    poly$type = numeric_cols[i]

    arc = rbind(arc, poly)
  }

  arc$type = factor(arc$type, levels = numeric_cols)

  ## Calculating Entropy and fitting a small linear model for prediction.
  a = 100 / (n - 1)
  data1 = data.frame(Entropy = n:1,
                     y = c(0, 1:(n - 1) * a))

  lm = lm(y ~ Entropy, data = data1)

  data = entropy(data)

  save_e = data$Entropy

  data_1 = Qentropy(data)

  numeric_cols1 = sapply(data_1, is.numeric)

  ## Determine column with minimum categorical entropy
  data_1$col = suppressWarnings(
    apply(data_1[, numeric_cols1], 1, function(row) {

      min_cols = colnames(data_1)[numeric_cols1][row == min(row)]

      if(length(min_cols) > 1) {
        min_cols = sample(min_cols, 1)
      }
      return(min_cols)
    })
  )


  ## Mapping entropy values to radial coordinates


  if(!is.null(col_variable) && length(col_variable) > 0){
    data = data |>
      select(Entropy, Factor) |>
      mutate(Entropy = 2 ^ Entropy)


    data$rad = predict(lm, newdata = data)
    data = data |> bind_cols(data_1[, ncol(data_1)])
    data = data |> left_join(location, join_by(col))
    data = data |> mutate(rand_deg = sample(seq(from=start, to=end, length.out = 10), 1))

    ## Computing final positions for each observation
    data = data |> select(Entropy, Factor, col, rad, deg, rand_deg) |> mutate(
      x = rad * cos(ifelse(rect1 == TRUE, deg, rand_deg)),
      y = rad * sin(ifelse(rect1 == TRUE, deg, rand_deg)),
      alpha = 1 - Entropy / n
    )

  } else {
    data = data |>
      select(Entropy) |>
      mutate(Entropy = 2 ^ Entropy)


    data$rad = predict(lm, newdata = data)
    data = data |> bind_cols(data_1[, ncol(data_1)])
    data = data |> left_join(location, join_by(col))
    data = data |> mutate(rand_deg = sample(seq(from=start, to=end, length.out = 10), 1))

    ## Computing final positions for each observation
    data = data |> select(Entropy, col, rad, deg, rand_deg) |> mutate(
      x = rad * cos(ifelse(rect1 == TRUE, deg, rand_deg)),
      y = rad * sin(ifelse(rect1 == TRUE, deg, rand_deg)),
      alpha = 1 - Entropy / n
    )

  }

  ## factor levels
  data$col = factor(data$col, levels = numeric_cols)

  data1 = data.frame(rad = a * (1:(n - 1)))

  ## Base plot

  polygon_colors <- if (!is.null(background_colors) && length(background_colors) > 0) {
    background_colors
  } else {
    scales::hue_pal()(n)
  }

  theme_set(theme_minimal())
  themx = theme_update(
    aspect.ratio = 1,
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.background = element_blank())

  base_plot = ggplot2::ggplot() +
    geom_polygon(
      data = arc,
      aes(x, y, fill = type),
      color = line_col,
      show.legend = FALSE,
      alpha = area
    ) +
    geom_circle(data = data1,
                aes(x0 = 0, y0 = 0, r = rad),
                col = line_col)  +
    geom_circle(aes(x0 = 0,
                    y0 =0,
                    r = data1$rad[nrow(data1)]),
                col = out_line, inherit.aes = FALSE) +
    coord_equal(xlim = c(-105,105), ylim = c(-105,105)) +
    scale_fill_manual(values = polygon_colors, na.value = 'whitesmoke' ) +
    ggtitle(title)


  if(!is.null(point_fill_colors) && length(point_fill_colors) > 0) {
    point_fill_colors
  } else {
    scales::hue_pal()(n)
  }


  if(!is.null(col_variable) && length(col_variable) > 0 && label == 'legend'){



    n2 = length(unique(data$Factor))

    point_fill = if (!is.null(point_fill_colors) && length(point_fill_colors) > 0) {
      point_fill_colors
    } else {
      scales::hue_pal()(n2)
    }

    point_line  = if (!is.null(point_line_colors) && length(point_line_colors) > 0) {
      point_line_colors
    } else {
      scales::hue_pal()(n2)
    }


    circle_plot = base_plot +
      new_scale('color') +
      new_scale('fill') +
      geom_jitter(data = data,
                  aes(x, y,
                      fill = Factor,
                      alpha = alpha,
                      col = Factor),
                  pch = 21,
                  size = sizex,
                  width = 3,
                  height = 3) +
      scale_fill_manual(values = point_fill, na.value = 'whitesmoke') +
      scale_color_manual(values = point_line, na.value = 'whitesmoke', guide = 'none') +
      scale_alpha(guide = 'none')

  }else if(!is.null(col_variable) && length(col_variable) > 0 && label == 'curve'){

    n2 = length(unique(data$Factor))

    point_fill = if (!is.null(point_fill_colors) && length(point_fill_colors) > 0) {
      point_fill_colors
    } else {
      scales::hue_pal()(n2)
    }

    point_line  = if (!is.null(point_line_colors) && length(point_line_colors) > 0) {
      point_line_colors
    } else {
      scales::hue_pal()(n2)
    }


    circle_plot = base_plot +
      new_scale('color') +
      new_scale('fill') +
      geom_jitter(data = data,
                  aes(x, y,
                      fill = Factor,
                      alpha = alpha,
                      col = Factor),
                  pch = 21,
                  size = sizex,
                  width = 3,
                  height = 3,
                  show.legend = FALSE) +
      geom_textcurve(data = location,
                     aes(x = x,
                         xend = xend,
                         y = y,
                         yend = yend,
                         label = labels),
                     linecolour = NA,
                     curvature = -0.4,
                     size = textsize,
                     na.rm = TRUE,
                     show.legend = FALSE) +
      scale_fill_manual(values = point_fill,
                        na.value = 'whitesmoke') +
      scale_color_manual(values = point_line,
                         na.value = 'whitesmoke', guides = 'none') +
      scale_alpha(guide = 'none')



  }else if(is.null(col_variable) && label == 'legend'){



    point_fill = if (!is.null(point_fill_colors) && length(point_fill_colors) > 0) {
      point_fill_colors
    } else {
      scales::hue_pal()(n)
    }

    point_line  = if (!is.null(point_line_colors) && length(point_line_colors) > 0) {
      point_line_colors
    } else {
      scales::hue_pal()(n)
    }


    circle_plot = base_plot +
      new_scale('color') +
      new_scale('fill') +
      geom_jitter(data = data,
                  aes(x, y,
                      fill = col,
                      alpha = alpha,
                      col = col),
                  pch = 21,
                  size = sizex,
                  width = 3,
                  height = 3) +
      scale_fill_manual(values = point_fill, na.value = 'whitesmoke') +
      scale_color_manual(values = point_line, na.value = 'whitesmoke', guide = 'none') +
      scale_alpha(guide = 'none')



  }else if(is.null(col_variable) && label == 'curve'){



    point_fill = if (!is.null(point_fill_colors) && length(point_fill_colors) > 0) {
      point_fill_colors
    } else {
      scales::hue_pal()(n)
    }

    point_line  = if (!is.null(point_line_colors) && length(point_line_colors) > 0) {
      point_line_colors
    } else {
      scales::hue_pal()(n)
    }

    circle_plot = base_plot +
      new_scale('color') +
      new_scale('fill') +
      geom_jitter(data = data,
                  aes(x, y,
                      fill = col,
                      alpha = alpha,
                      col = col),
                  pch = 21,
                  size = sizex,
                  width = 3,
                  height = 3,
                  show.legend = FALSE) +
      geom_textcurve(data = location,
                     aes(x = x,
                         xend = xend,
                         y = y,
                         yend = yend,
                         label = labels),
                     linecolour = NA,
                     curvature = -0.4,
                     size = textsize,
                     na.rm = TRUE,
                     show.legend = FALSE) +
      scale_fill_manual(values = point_fill,
                        na.value = 'whitesmoke') +
      scale_color_manual(values = point_line,
                         na.value = 'whitesmoke', guide = 'none') +
      scale_alpha(guide = 'none')
  } else {

    print('There is something wrong with the data. Returning NULL')
    return(NULL)
  }



  data$Entropy = save_e



  return(list(circle_plot, data))

}
