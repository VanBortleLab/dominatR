
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
