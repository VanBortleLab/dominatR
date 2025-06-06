#' Circular Dominance Plot (N ≥ 3 variables)
#'
#' @description
#' Produces a radial dominance plot in which each observation is located by:
#' \itemize{
#'   \item \strong{Angle (θ)} – the variable with the greatest value
#'         (ties broken at random).
#'   \item \strong{Radius (r)} – a monotone mapping of the row‐wise Shannon
#'         entropy: points with low entropy (one variable dominates)
#'         are near the edge; points with high entropy lie toward the centre.
#' }
#' The circle is partitioned into \eqn{n} coloured slices; optional factor
#' information can colour/jitter points independently.  Labels for each
#' slice may be drawn as curved text on the circle or shown in a legend.
#'
#' @param x A numeric \code{data.frame}, \code{matrix}, or
#'   a \code{SummarizedExperiment}.
#' @param n Integer (\eqn{\ge 3}). How many numeric variables to visualise.
#'   Must match \code{length(column_variable_factor)} when supplied.
#' @param column_variable_factor Character. Name of a column (or rowData
#'   column in a SummarizedExperiment) holding a categorical variable whose
#'   levels will colour the points.  If \code{NULL} (default) points are
#'   coloured by their dominant variable.
#' @param variables_highlight Character vector naming which variables
#'   should receive curved text labels when \code{label = "curve"}.
#'   Defaults to all variables.
#' @param entropyrange,magnituderange Numeric length-2 vectors.
#'   Rows falling outside either interval are excluded from the plot/data.
#' @param background_alpha_polygon Alpha level (0–1) for the coloured
#'   background slices.
#' @param background_polygon Character vector of slice fill colours;
#'   defaults to \code{scales::hue_pal()(n)}.  \code{background_na_polygon}
#'   sets the colour for missing values.
#' @param point_size Numeric; plotted point size.
#' @param point_fill_colors,point_line_colors Optional colour vectors for
#'   point fill / outline.
#' @param background_na_polygon,point_fill_na_colors,point_line_na_colors  Sets the colour for missing values.
#' @param line_col Colour for the inner grid / slice borders.
#' @param out_line Colour for the outermost circle.
#' @param straight_points Logical. If TRUE points are plotted in a straight line.
#' @param label Either \code{"legend"} (default) to list variables in a
#'   legend or \code{"curve"} to print them around the rim.
#' @param text_label_curve_size Numeric font size for curved labels.
#' @param assay_name (SummarizedExperiment only) Which assay to use.
#'   Defaults to the first assay.
#' @param output_table Logical.  Also return the underlying data frame?
#'
#' @return
#' If \code{output_table = TRUE} a list with:
#' \itemize{
#'   \item \code{circle_plot} — a \link[ggplot2]{ggplot} object;
#'   \item \code{data}        — the augmented data frame containing
#'         entropy, radius, (x,y) coordinates, dominant variable and
#'         optional factor.
#' }
#' Otherwise only the \code{ggplot} object is returned.
#'
#' @details
#' \subsection{Radius mapping}{
#' A linear map is used
#' \deqn{ r \;=\; 100 \,\frac{n - 2^{H}}{n-1} }
#' where \eqn{H} is the Shannon entropy of the row after log base 2, so
#' \eqn{H \in [0,\log_2 n]}.
#' }
#'
#' @import SummarizedExperiment tidyverse forcats ggforce ggnewscale lubridate purrr readr stringr tibble tidyr dplyr ggplot2 utils
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' library(airway)
#' library(tidyverse)
#' data('airway')
#' se = airway
#'
#' ## Normalize the data first using tpm_normalization
#' rowData(se)$gene_length = rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
#'
#' se = tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
#'
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'
#' df <- assay(se, 'tpm_norm') |> as.data.frame()
#'
#' ## For simplicity let's rename the columns
#' colnames(df) <- paste('Column_', 1:8, sep ='')
#'
#' # Default
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   entropyrange     = c(0, 3),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend', output_table = FALSE
#' )
#'
#' # Filtering by entropy, 8 variables, max entropy value is log2(8)
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   entropyrange     = c(2, 3),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend', output_table = FALSE
#' )
#'
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend', output_table = FALSE
#' )
#'
#' # Aesthetics modification
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'curve',
#'   output_table = FALSE
#' )
#'
#' # It is possible to highlight only a specific variable
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   background_alpha_polygon = 0.2,
#'   background_na_polygon = 'transparent',
#'   background_polygon = c('Column_1'  = 'indianred',
#'                          'Column_3' = 'lightblue',
#'                          'Column_5' = 'lightgreen'),
#'   point_fill_colors = c('Column_1'  = 'darkred',
#'                         'Column_3' = 'darkblue',
#'                         'Column_5' = 'darkgreen'),
#'   point_line_colors = c('Column_1'  = 'black',
#'                         'Column_3' = 'black',
#'                         'Column_5' = 'black')
#' )
#'
#' # Let's create a factor column in our df
#' df$factor <- sample(c('A', 'B', 'C', 'D'), size = nrow(df), replace = TRUE)
#'
#' # It is possible to visualize things by this specific factor column using
#' # column_variable_factor
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   column_variable_factor = 'factor',
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   background_alpha_polygon = 0.2,
#'   background_na_polygon = 'transparent',
#'   background_polygon = c('Column_1'  = 'indianred',
#'                          'Column_3' = 'lightblue',
#'                          'Column_5' = 'lightgreen')
#' )
#'
#' # Colors can be modified
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   column_variable_factor = 'factor',
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'curve',
#'   output_table = FALSE,
#'   background_alpha_polygon = 0.02,
#'   background_na_polygon = 'transparent',
#'   point_fill_colors = c('A' = 'black',
#'                         'B' = 'gray',
#'                         'C' = 'white',
#'                         'D' = 'orange'),
#'   point_line_colors = c('A' = 'black',
#'                         'B' = 'gray',
#'                         'C' = 'white',
#'                         'D' = 'orange')
#' )
#'
#' # Size of the points can be modified too
#' plot_circle(
#'   x = df,
#'   n = 8,
#'   point_size =  2,
#'   column_variable_factor = 'factor',
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'curve',
#'   output_table = FALSE,
#'   background_alpha_polygon = 0.02,
#'   background_na_polygon = 'transparent',
#'   point_fill_colors = c('A' = 'black',
#'                         'B' = 'gray',
#'                         'C' = 'white',
#'                         'D' = 'orange'),
#'   point_line_colors = c('A' = 'black',
#'                         'B' = 'gray',
#'                         'C' = 'white',
#'                         'D' = 'orange')
#' )
#'
#' # Retrieving a dataframe with the results used for plotting, set output_table = TRUE
#' plot <- plot_circle(
#'   x = df,
#'   n = 8,
#'   point_size =  2,
#'   column_variable_factor = 'factor',
#'   entropyrange     = c(0, 2),
#'   magnituderange   = c(0, Inf),
#'   label  = 'curve',
#'   output_table = TRUE,
#'   background_alpha_polygon = 0.02,
#'   background_na_polygon = 'transparent',
#'   point_fill_colors = c('A' = 'black',
#'                         'B' = 'gray',
#'                         'C' = 'white',
#'                         'D' = 'orange'),
#'   point_line_colors = c('A' = 'black',
#'                         'B' = 'gray',
#'                         'C' = 'white',
#'                         'D' = 'orange')
#' )
#'
#'
#' # The first object is the plot
#' plot[[1]]
#'
#' # The second the dataframe with information for each row, including
#' # Entropy and the variable that dominates that particular observation.
#'
#'
#' head(plot[[2]])
#'
#'
#'
#' # -------------------------------
#' # 1) Using a SummarizedExperiment
#' # -------------------------------
#' # Changing column names
#' colnames(se) <- paste('Column_', 1:8, sep ='')
#'
#' # Default
#' plot_circle(
#'   x = se,
#'   n = 8,
#'   entropyrange     = c(0, 3),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   assay_name = 'tpm_norm'
#' )
#'
#' # Filtering High Entropy genes
#' plot_circle(
#'   x = se,
#'   n = 8,
#'   entropyrange     = c(0, 1.5),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   assay_name = 'tpm_norm'
#' )
#'
#' # Filtering Low Entropy genes
#' plot_circle(
#'   x = se,
#'   n = 8,
#'   entropyrange     = c(2, 3),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   assay_name = 'tpm_norm'
#' )
#'
#'
#' # Using a character column from rowData
#'
#' plot_circle(
#'   x = se,
#'   n = 8,
#'   column_variable_factor = 'gene_biotype',
#'   entropyrange     = c(2,3),
#'   magnituderange   = c(0, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   assay_name = 'tpm_norm'
#' )
#'
#' plot_circle(
#'   x = se,
#'   n = 8,
#'   column_variable_factor = 'gene_biotype',
#'   point_size = 3,
#'   entropyrange     = c(0,1.5),
#'   magnituderange   = c(2, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   assay_name = 'tpm_norm',
#' )
#'
#' # Highlighting only a class of interest
#'
#' plot_circle(
#'   x = se,
#'   n = 8,
#'   column_variable_factor = 'gene_biotype',
#'   point_size = 3,
#'   entropyrange     = c(0,1.5),
#'   magnituderange   = c(2, Inf),
#'   label  = 'legend',
#'   output_table = FALSE,
#'   assay_name = 'tpm_norm',
#'   point_fill_colors = c('miRNA' = 'orange'),
#'   point_line_colors = c('miRNA' = 'orange')
#' )
#'
#'
#' # Retrieving a dataframe with the results used for plotting, set output_table = TRUE
#'
#' plot <- plot_circle(
#'   x = se,
#'   n = 8,
#'   column_variable_factor = 'gene_biotype',
#'   point_size = 3,
#'   entropyrange     = c(0,1.5),
#'   magnituderange   = c(2, Inf),
#'   label  = 'legend',
#'   output_table = TRUE,
#'   assay_name = 'tpm_norm',
#'   point_fill_colors = c('miRNA' = 'orange'),
#'   point_line_colors = c('miRNA' = 'orange')
#' )
#'
#' # It returns a list.
#' # The first object is the plot
#' plot[[1]]
#'
#' # The second the dataframe with information for each row, including
#' # Entropy and the variable that dominates that particular observation.
#' head(plot[[2]])
#'


plot_circle <- function(x,
                        n,
                        column_variable_factor = NULL,
                        variables_highlight = NULL,
                        entropyrange     = c(0, Inf),
                        magnituderange   = c(0, Inf),
                        background_alpha_polygon = 0.05,
                        background_polygon = NULL,
                        background_na_polygon = 'whitesmoke',
                        point_size = 1,
                        point_fill_colors = NULL,
                        point_fill_na_colors = 'whitesmoke',
                        point_line_colors = NULL,
                        point_line_na_colors = 'whitesmoke',
                        straight_points = TRUE,
                        line_col = 'gray90',
                        out_line = 'black',
                        label  = c("legend", "curve"),
                        text_label_curve_size = 3,
                        assay_name       = NULL,
                        output_table     = TRUE)
{
  ## ---- 1. Pulling the data and the column factor if needed.
  if (inherits(x, "SummarizedExperiment")) {
    if (is.null(assay_name))
      assay_name <- SummarizedExperiment::assayNames(x)[1]
    mat <- SummarizedExperiment::assay(x, assay_name)
    if (!is.numeric(mat)) stop("Chosen assay is not numeric.")
    mat <- as.data.frame(mat)

    if(!is.null(column_variable_factor)){
      mat <- cbind(mat, rowData(x)[which(colnames(rowData(se)) == column_variable_factor)])
      colnames(mat)[ncol(mat)] = 'Factor'

    }
  } else {
    mat <- as.data.frame(x)

    if(!is.null(column_variable_factor)){
      colnames(mat)[which(colnames(mat) == column_variable_factor)] <- 'Factor'
      mat <- mat |> relocate(Factor)

    }
  }


  ## Plot parameters
  area <- background_alpha_polygon
  sizex <- point_size
  textsize <- text_label_curve_size
  rect1 <- straight_points

  if(is.null(variables_highlight)){
    variables_highlight <- colnames(mat)
  } else {
    variables_highlight
  }

  colnames(mat) <- substr(colnames(mat), 1, 10)
  variables_highlight <- substr(variables_highlight,1,10)

  ## ---- 2. Key Plotting Settings
  ## Key Plot setting
  ## Whole circle has an perimeter of 2pi.
  ## Constant Parameters
  inner_r <- ifelse(n > 15, 80, 70)   # deg inner arc radius
  outer_r <- ifelse(n > 15, 100, 110) # deg outter arc radius
  deg <- 2 * pi / n                   # width of each slice
  deg_sp <- deg / 2                   # half a slice
  rad_label <- ifelse(n == 3, 120, 100.5) # Label position
  labels1 <- variables_highlight


  ## Mid_angle, start, end. Radians
  deg_f <- pi/2 - deg * seq_len(n)
  start <- inner_r * pi/180 - deg_sp - deg * seq(0, n-1)
  end <- outer_r * pi/180 - deg_sp - deg * seq_len(n)
  curv_start <- pi/2 - deg_sp - deg * seq(0, n - 1)
  curv_end <-  pi/2 - deg_sp - deg * seq_len(n)
  x <- rad_label * cos(curv_start)
  xend <- rad_label * cos(curv_end)
  y <- rad_label * sin(curv_start)
  yend <- rad_label * sin(curv_end)

  location <- data.frame(
    col = colnames(select_if(mat, is.numeric)),
    deg = deg_f,
    start = start,
    end = end,
    x = x,
    xend = xend,
    y = y,
    yend = yend
  )


  ## For those specific variables we want to label
  location$labels <- NA
  location$labels[location$col %in% labels1] <- location$col[location$col %in%
                                                               labels1]

  ## ---- 3. Ploting coordinates for the each variable slice.

  numeric_cols <- colnames(mat)[vapply(mat, is.numeric, logical(1))]

  arc <- NULL
  for (i in seq(1, n)) {
    p <- ggplot() +
      geom_arc(aes(
        x0 = 0,
        y0 = 0,
        r = 100,
        start = deg_sp + deg * (i - 1),
        end = deg_sp + deg * (i)
      ))
    poly <- rbind(c(0, 0),
                  data.frame(x = ggplot_build(p)$data[[1]]$x, y = ggplot_build(p)$data[[1]]$y),
                  c(0, 0))
    poly$type <- numeric_cols[i]

    arc <- rbind(arc, poly)
  }

  arc$type <- factor(arc$type, levels = numeric_cols)




  ## ---- 4. Entropy calculation and filtering
  mat_num <- mat[,numeric_cols]
  range_values <- do.call(pmax, as.data.frame(mat_num))
  mat_entropy <- entropy(mat)
  entropy_values <- mat_entropy$Entropy

  a <- 100 / (n - 1)
  data1 <- data.frame(Entropy = seq(n, 1), y = c(0, seq(1, (n - 1)) * a))
  lm <- lm(y ~ Entropy, data = data1)
  rad <- predict(lm, newdata = mat_entropy |> mutate(Entropy = 2 ^ Entropy) )

  dominant_col <- max.col(mat_entropy[,numeric_cols], ties.method = 'random')
  dominant_col <- numeric_cols[dominant_col]

  mat <- entropy(mat)

  mat <- mat |> mutate(Entropy = entropy_values,
                       Range = range_values,
                       col = dominant_col,
                       rad = rad)

  # Filtering based on our range of interest
  mat <- mat |> filter(Entropy >= entropyrange[1] &
                         Entropy <= entropyrange[2] &
                         range_values >= magnituderange[1] &
                         range_values <= magnituderange[2])

  rn  <- row.names(mat)

  ## ---- 5. Creating the final dataframe
  # Joining location
  mat <- mat |> select(!all_of(numeric_cols)) |>
    left_join(location, join_by(col))

  mat <- mat |>
    rowwise() |>
    mutate(rand_deg = sample(seq(from = start, to = end, length.out = 10), 1),
           x = rad * cos(ifelse(rect1 == TRUE, deg, rand_deg)),
           y = rad * sin(ifelse(rect1 == TRUE, deg, rand_deg)),
           alpha = 1 - Entropy / n
    )


  mat <- mat |> select(!any_of(c('Range', 'start', 'end', 'xend', 'yend'))) |>
    as.data.frame()

  mat$col <- factor(mat$col, levels = numeric_cols)


  row.names(mat) <- rn

  ## ---- 6. Plotting
  a <- 100 / (n - 1)
  data1 <- data.frame(rad = a * (seq(1, (n - 1))))

  ### Colors for the background of each slice
  background_polygon <- if (!is.null(background_polygon) &&
                            length(background_polygon) > 0) {
    background_polygon
  } else {
    scales::hue_pal()(n)
  }

  theme_set(theme_minimal())
  themx <- theme_update(
    aspect.ratio = 1,
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.background = element_blank()
  )

  base_plot <- ggplot2::ggplot() +
    geom_polygon(
      data = arc,
      aes(x, y, fill = type),
      color = line_col,
      show.legend = FALSE,
      alpha = area
    ) +
    geom_circle(data = data1, aes(x0 = 0, y0 = 0, r = rad), col = line_col)  +
    geom_circle(aes(x0 = 0, y0 = 0, r = data1$rad[nrow(data1)]),
                col = out_line,
                inherit.aes = FALSE) +
    coord_equal(xlim = c(-105, 105), ylim = c(-105, 105)) +
    scale_fill_manual(name = 'Variable',
                      values = background_polygon,
                      na.value = background_na_polygon)


  ## Mapping decisions and colors
  has_factor <- !is.null(column_variable_factor)
  map_var      <- if (has_factor) quo(Factor) else quo(col)

  level <- if(has_factor) unique(mat$Factor) else unique(mat$col)
  n_level <- length(level)

  fill_pal <- if (!is.null(point_fill_colors) &&
                  length(point_fill_colors) > 0){
    point_fill_colors
  } else {
    scales::hue_pal()(n_level)}


  line_pal <- if (!is.null(point_line_colors) &&
                  length(point_line_colors) > 0){
    point_line_colors
  }else{
    fill_pal }


  ## Core plot
  circle_plot <- base_plot +
    new_scale_fill() +
    new_scale_color() +
    geom_jitter(
      data = mat,
      aes(
        x,
        y,
        fill = !!map_var,
        col = !!map_var,
        alpha = alpha
      ),
      pch = 21,
      size = sizex,
      width = 3,
      height = 3
    ) +
    scale_fill_manual(values = fill_pal,
                      na.value = point_fill_na_colors
    ) +
    scale_color_manual(values = line_pal,
                       na.value = point_line_na_colors,
                       guide = 'none'
    ) +
    scale_alpha(guide = 'none')

  ## Ad label as curve if needed
  if(label == 'curve') {
    circle_plot  <- circle_plot +
      geomtextpath::geom_textcurve(
        data = location,
        aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend,
          label = labels
        ),
        linecolour = NA,
        curvature = -0.4,
        size = textsize,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      guides(fill = ifelse(!is.null(column_variable_factor), 'legend', 'none'))
  }




  ## Final
  if (output_table){ list(circle_plot, mat)} else {circle_plot}
}
