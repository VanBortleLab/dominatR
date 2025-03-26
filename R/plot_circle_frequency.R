#' Dominance Distribution Plotter
#'
#' This function depends on the plot_circle function when the parameter col_variable is used. It plots the dominance frequency distribution of each of the categories present in col_variable.
#'
#' @param n Number of Variables, must be the same as the ones used in the function
#' @param circle An object type list obtained from the plot_circle function
#' @param single A logical parameter. If TRUE a single plot will be produced otherwise a plot will be faceted
#' @param legend A logical parameter. If TRUE legends will appear
#' @param numb_columns A numerical value. Number of columns for plotting when single is FALSE
#' @param filter_class Name or names of classes of interest for display
#' @param point_size Size of the point
#'
#' @returns
#' Returns a list of objects, a frequency distribution plot for each category and the respective dataframe
#'
#' @export
#' @seealso [entropy()], [Qentropy()], [plot_circle()]
#'
plot_circle_frequency = function(n,
                                 circle,
                                 single = FALSE,
                                 legend = TRUE,
                                 numb_columns = 1,
                                 filter_class = NULL,
                                 point_size = 2){
  data_ent = circle[[2]]

  breaks_log2 = log2(seq(0,n) + 0.5)


  data_ent$bin = cut(as.numeric(data_ent$Entropy),
                     breaks =breaks_log2,
                     labels =  2^(breaks_log2[-length(breaks_log2)]) |> ceiling(),
                     include.lowest = FALSE, right = TRUE)

  data_ent = data_ent |>
    group_by(bin, Factor) |>
    dplyr::select(bin, Factor) |>
    dplyr::summarise(n = dplyr::n(), .groups = 'drop')

  data_ent = data_ent |>
    group_by(Factor) |>
    mutate(proportion = n/sum(n))

  data_ent$Factor = factor(data_ent$Factor)


  if(!is.null(filter_class)){
    data_ent = data_ent |>
      dplyr::filter(Factor %in% filter_class)
  }

  plot_stat = data_ent |>
    ggplot() +
    geom_line(aes(x = bin,
                  y = proportion,
                  col = Factor,
                  group = Factor),
              linewidth = 1,
              show.legend = legend) +
    geom_point(aes(x = bin,
                   y = proportion,
                   fill = Factor,
                   group = Factor),
               show.legend = FALSE,
               pch = 21,
               size = point_size) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(colour = 'gray95', linetype = 'dotted'),
      panel.spacing.y = unit(0.4, "cm"),
      panel.spacing.x = unit(0.2, 'cm'),
      strip.placement = "outside",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face =
                                    'bold'),
      legend.key.size = unit(15, 'points'),
      axis.title.y = element_text( size = 20, color = 'black', face = 'bold', vjust = 1),
      axis.title.x = element_text(size = 20, color = 'black', face = 'bold'),
      axis.text.x = element_text(size = 15,color = 'black',hjust = 0.5),
      axis.text.y = element_text(size = 15, color = 'black'),
      axis.ticks = element_line(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 1
      )
    ) +
    xlab('Dominance') +
    scale_y_continuous(limits = c(0,1), n.breaks =5) +
    scale_x_discrete(limits = factor(seq(1, n)) )

  if(single == FALSE){
    plot_stat = plot_stat + facet_wrap(~Factor, ncol = numb_columns)
  }
  return(list(plot_stat, data_ent))
}
