#' Dominance Distribution Plotter
#'
#' This function depends on the plot_circle function when the parameter col_variable is used. It plots the dominance frequency distribution of each of the categories present in col_variable.
#'
#' @param n Number of Variables, must be the same as the ones used in the function
#' @param circle An object type list obtained from the plot_circle function
#' @param single A logical parameter. If TRUE a single plot will be produced otherwise a plot will be faceted
#' @param legend A logical parameter. If TRUE legends will appear
#' @param numb_columns A numerical value. Number of columns for plotting when single is FALSE
#'
#' @returns
#' Returns a list of objects, a frequency distribution plot for each category and the respective dataframe
#'
#' @export
#' @seealso [entropy()], [Qentropy()], [plot_circle()]
#'
plot_circle_frequency = function(n, circle, single = F, legend = T, numb_columns = 1){
  data_ent = circle[[2]]

  breaks_log2 = log2(0:n + 0.5)


  data_ent$bin = cut(as.numeric(data_ent$Entropy),
                     breaks =breaks_log2,
                     labels =  2^(breaks_log2[-length(breaks_log2)]) |> ceiling(),
                     include.lowest = F, right = T)

  data_ent = data_ent |> group_by(bin, Factor) |> dplyr::select(bin, Factor) |> summarise(n = n(), .groups = 'drop')
  data_ent = data_ent |> group_by(Factor) |> mutate(proportion = n/sum(n))
  data_ent$Factor = factor(data_ent$Factor)

  plot_stat = data_ent |>
    ggplot() +
    geom_line(aes(x = bin,
                  y = proportion,
                  col = Factor,
                  group = Factor), linewidth = 1,
              show.legend = legend) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor  = element_blank(),
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
    xlab('Dominance')

  if(single == F){
    plot_stat = plot_stat + facet_wrap(~Factor, ncol = numb_columns)
  }
  return(list(plot_stat, data_ent))
}
