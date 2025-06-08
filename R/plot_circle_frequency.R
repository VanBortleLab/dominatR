#' Dominance–Entropy Frequency Plot
#'
#' @description
#' Visualises how often each categorical level ( `Factor` ) is dominant at a
#' given entropy score.  The function expects the \strong{second} element of the
#' list returned by \code{\link{plot_circle}()}.
#'
#' @param n              Integer. Number of numeric variables used in
#'   \code{plot_circle()}.
#' @param circle         The list returned by \code{plot_circle()}.
#' @param single         Logical.  If \code{TRUE} draw one combined panel;
#'   otherwise facet by \code{Factor}.
#' @param legend         Logical.  Show a legend for the plot
#' @param numb_columns   Faceting columns when \code{single = FALSE}.
#' @param filter_class   Character vector of levels to keep; \code{NULL} keeps all.
#' @param point_size     Numeric.  Size of jitter points.
#'
#' @return
#' A list with
#' \itemize{
#'   \item \code{plot_stat} — a \link[ggplot2]{ggplot} object.
#'   \item \code{data}      — the aggregated frequency table.
#' }
#' @export
#'
#' @seealso \code{\link{plot_circle}}
#'
#' @examples
#'library(SummarizedExperiment)
#'library(airway)
#'data('airway')
#'se = airway
#'
#'## Normalize the data first using tpm_normalization
#'rowData(se)$gene_length = rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
#'
#'se = tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
#'
#'## Creating a plot_circle list using the 'gene_biotype' column as factor
#'plot_test <- plot_circle(
#'  x = se,
#'  n = 8,
#'  column_variable_factor = 'gene_biotype',
#'  entropyrange     = c(0,Inf),
#'  magnituderange   = c(0, Inf),
#'  label  = 'legend',
#'  output_table = TRUE,
#'  assay_name = 'tpm_norm'
#')
#'
#'## Using the plot_test object created above
#'## Default
#'plot <- plot_circle_frequency(n = 8,
#'                      circle = plot_test,
#'                      single = TRUE,
#'                      legend = TRUE,
#'                      numb_columns = 1,
#'                      filter_class = NULL,
#'                      point_size = 2)
#'
#'plot[[1]]
#'
#'## Facetting by factor is possible, adjusting the number of columns
#'plot <- plot_circle_frequency(n = 8,
#'                      circle = plot_test,
#'                      single = FALSE,
#'                      legend = TRUE,
#'                      numb_columns = 3,
#'                      filter_class = NULL,
#'                      point_size = 2)
#'
#'plot[[1]]
#'
#'## Subsetting by a specific class present in Factor
#'plot_circle_frequency(n = 8,
#'                      circle = plot_test,
#'                      single = FALSE,
#'                      legend = TRUE,
#'                      numb_columns = 1,
#'                      filter_class = c('protein_coding', 'snoRNA', 'miRNA'),
#'                      point_size = 2)
#'
#'plot[[1]]
#'
plot_circle_frequency <- function(n,
                                  circle,
                                  single        = FALSE,
                                  legend        = TRUE,
                                  numb_columns  = 1,
                                  filter_class  = NULL,
                                  point_size    = 2)
{
  ## ---------------- 1. sanity checks ----------------------------------------
  stopifnot(length(circle) >= 2,
            is.list(circle),
            "data.frame" %in% class(circle[[2]]))
  dat <- circle[[2]]

  if (!"Factor" %in% colnames(dat)){
    stop("`plot_circle_frequency()` requires the `plot_circle()` call ",
         "to include a categorical column called 'Factor'.")
  }
  if (!is.numeric(dat$Entropy)){
    stop("Second element of `circle` does not contain numeric 'Entropy'.")
  }

  ## ---------------- 2. bin Entropy in log2 bands ----------------------------
  breaks   <- log2(seq(0, n) + 0.5)
  labels   <- ceiling(2^(breaks[-length(breaks)]))                         # 1,2,4,8,…
  dat$bin  <- factor(
    cut(dat$Entropy,
        breaks      = breaks,
        labels      = labels,
        include.lowest = FALSE, right = TRUE),
    levels = labels
  )

  ## ---------------- 3. fast frequency & proportion table --------------------
  tab <- as.data.frame.table(table(dat$bin, dat$Factor),
                             responseName = "n",
                             stringsAsFactors = FALSE)
  names(tab) <- c("bin", "Factor", "n")
  tab$bin    <- factor(tab$bin, levels = labels)

  # proportions within each Factor
  totals     <- tapply(tab$n, tab$Factor, sum)
  tab$proportion <- tab$n / totals[tab$Factor]


  ## optional filtering
  if (!is.null(filter_class)){
    tab <- tab[tab$Factor %in% filter_class, , drop = FALSE]
  }

  ## ---------------- 4. build ggplot -----------------------------------------
  p <- ggplot(tab,
              aes(x = bin, y = proportion,
                  group = Factor, colour = Factor, fill = Factor)) +
    geom_line(linewidth = 1, show.legend = legend) +
    geom_point(data = tab |> filter(proportion > 0),
               aes(x = bin, y = proportion),
               shape = 21,  size = point_size, show.legend = FALSE, col = 'black') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
    scale_x_discrete(limits = factor(seq(1, n))) +
    labs(x = "Dominance", y = "Proportion") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(colour = "grey92", linetype = "dotted"),
      panel.border     = element_rect(colour = "black", fill = NA),
      legend.title     = element_blank()
    )

  if (!single){
    p <- p + facet_wrap(~ Factor, ncol = numb_columns)
  }

  return(list(plot_stat = p, data = tab))
}
