#' Dominance Classification Plotter.
#'
#' Creates an abacus plot that aims to classify different observations in
#' determined number of groups across variables using an entropy-based method
#'
#' @param data A numerical dataframe. Each column represent the variable used
#' to calculate the entropy and build the plot
#' @param n Total number of variables to account for entropy calculation.
#' @param x_variable A reference column, this is the variable that will be
#' classified. A column with unique gene names or unique index for genes.
#' @param y_variables Columns of interest for visualization purposes.
#' @param title Title for the plot
#' @param percentiles Number of groups intended for classification. If 4 is
#' selected you will obtain the percentiles 0.25, 0.5, 0.75 and 1.00.
#' @param entropy_threshold A numerical value that determines the entropy
#' threshold for visualization.
#' @param point_fill_colors Filling colors for the points
#' @param point_line_colors Outline colors for the points
#' @param point_size Point Size
#' @param single A logical parameter. If TRUE a single plot will be produced.
#' @param numb_columns A numerical value. Number of columns for plotting when
#' single is FALSE
#' @param features_of_interest A vector with the names of genes or index of
#' genes of interest that should be plotted. This names should be present in
#' the x_variable column and have the same format.
#'
#' @returns
#' Returns a list of objects, an abacus classification plot and the respective
#' dataframe
#'
#' @export
#'
#'
plot_abacus <- function(data,
                        n,
                        x_variable,
                        y_variables,
                        title = NULL,
                        percentiles,
                        entropy_threshold = NULL,
                        point_fill_colors = NULL,
                        point_line_colors = NULL,
                        point_size = 2,
                        single = FALSE,
                        numb_columns = 1 ,
                        features_of_interest = NULL) {
    X_axis <- NULL
    Variable <- NULL

    point_fill <- if (!is.null(point_fill_colors) &&
                    length(point_fill_colors) > 0) {
        point_fill_colors
    } else {
        scales::hue_pal()(percentiles)
    }

        point_line  <- if (!is.null(point_line_colors) &&
                    length(point_line_colors) > 0) {
        point_line_colors
    } else {
        scales::hue_pal()(percentiles)
    }


        entropy_threshold <- if (!is.null(entropy_threshold) &&
                                length(entropy_threshold) > 0) {
        entropy_threshold
    } else {
        c(0, log2(n))
    }


    ## Calculating threshold for categorical entropy
    qentropies <- matrix(seq(1 / percentiles, 1, by = 1 / percentiles)) * n

    range <- log2(qentropies) - log2(1 / qentropies)

    range <- format(round(range, 2), 2)

    ### Calculating entropy
    data <- entropy(data)

    ### Looking for variable of interest for the x axis
    colnames(data)[which(colnames(data) == x_variable)] <- 'X_axis'

    ### subsetting for those
    entropy_ind <- data |>
        dplyr::select(X_axis, Entropy)

    entropy_ind <- entropy_ind |> dplyr::filter(Entropy >= entropy_threshold[1]
                                            & Entropy <= entropy_threshold[2])


    ### Calculating categorical entropy
    data <- Qentropy(data)

    ### Subsetting for y variables of interest
    y_vars <- which(colnames(data) %in% y_variables)

    ## Pivoting
    data <- data |>
        select(X_axis, all_of(y_vars)) |>
        pivot_longer(!X_axis, names_to = 'Variable', values_to = 'Qentropy')


    data$Qentropy[is.infinite(data$Qentropy)] <- NA

    #data$Qentropy = max_qentropy - data$Qentropy

    ## Transforming the data into groups
    data$bin <- cut(
        data$Qentropy,
        breaks = c(0, range),
        labels = format(round(
        seq(1 / percentiles, 1, by = 1 / percentiles), 2
    ), 2),
        include.lowest = TRUE,
        right = FALSE
    )

    data$Variable <- factor(data$Variable, levels = y_variables)
    data$bin <- factor(data$bin, levels = format(round(
        seq(1 / percentiles, 1, by = 1 / percentiles), 2
    ), 2))

    data <- data |> dplyr::filter(X_axis %in% entropy_ind$X_axis)
    data <- stats::na.omit(data)

    data <- data[order(data$Qentropy), ]

    if (!is.null(features_of_interest)) {
        data <- data |> dplyr::filter(X_axis %in% features_of_interest)
    }

    plot_heat <- data |> ggplot() +
        geom_point(aes(
        x = X_axis,
        y = Variable,
        fill = bin,
        col = bin
        ),
        pch = 21,
        size = point_size) +
        theme_minimal() +
        ggtitle(title) +
        theme(
            axis.text.x  = element_blank(),
            panel.border = element_rect(
            colour = 'black',
            fill = NA,
            linewidth = 2
        ),
        panel.grid.major.y = element_line(colour = 'black'),
        panel.grid.major.x = element_blank()
    ) +
        scale_size(range = c(0, 3)) +
        scale_x_discrete(expand = c(0.05, 0.05)) +
        scale_y_discrete(expand = c(0.1, 0.1), limits = y_variables) +
        scale_fill_manual(values = point_fill, na.value = NA) +
        scale_color_manual(values = point_line, na.value = NA) +
        guides(fill = guide_legend(title = 'Percentile'),
            color = guide_legend(title = 'Percentile')) +
        xlab('Genes') + ylab('')

    if (single == FALSE) {
        plot_heat <- plot_heat +
            facet_wrap( ~ bin,
                    ncol = numb_columns,
                    scales = 'free_x',
                    strip.position = 'right')
    }

    return(list(plot_heat, data))

}
