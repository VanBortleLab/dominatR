#' Rope (binary) dominance plot
#'
#' @description
#' Creates a rope-like visualization comparing two numeric columns
#' (e.g., "a" vs. "b"), with optional color filtering based on
#' maximum value range and entropy range.
#'
#' The plot is useful for visualising “winner-takes-all” behaviour in two-way
#' comparisons, e.g. gene expression in *A* and *B* conditions.
#'
#' @param x A \code{data.frame} or \code{matrix} with numeric columns, or
#'   a \code{SummarizedExperiment} containing such data in one of its assays.
#' @param column_name Character. The name of the two variables that will be
#' used for the analysis. By default it is \code{NULL}.
#' @param push_text Numeric. Expands or contracts text label positions along
#' the x-axis.
#' @param rope_width Numeric. Thickness of the "rope" drawn in the center.
#' @param rope_color Character. Color for the rope's fill.
#' @param rope_border Logical or a color. Whether or how to draw the rope
#' border.
#' @param col Character vector of length 2. Colors assigned when \code{a > b}
#'  or \code{b > a}, respectively.
#' @param col_bg Background color (used when a row is filtered out by entropy
#' or max value).
#' @param pch Integer or vector specifying point types for the main two
#' categories.
#' @param pch_bg Integer specifying the point type for the "gray" points
#' (if \code{plotAll=TRUE}).
#' @param cex Numeric. Expansion factor for point size.
#' @param entropyrange Numeric vector of length 2. Rows with \code{entropy}
#' outside this range become background color.
#' @param maxvaluerange Numeric vector of length 2. Rows with \code{max(a,b)}
#' outside this range become background color.
#' @param plotAll Logical. If \code{TRUE}, also draw "filtered" points in
#' \code{col_bg} color. If \code{FALSE}, only highlight active points.
#' @param label Logical. If \code{TRUE}, label the two columns near the rope
#' ends.
#' @param output_table Logical. If \code{TRUE}, return the processed data
#' frame with added columns.
#' @param assay_name (SummarizedExperiment only) Name of the assay containing
#'   the 2-column data. If not specified, the first assay is used.
#'
#' @return
#'   \itemize{
#'     \item If \code{output_table=TRUE}, returns a data frame with extra
#'       columns (\code{comx}, \code{comy}, \code{color}, \code{maxvalue},
#'       \code{entropy}) used in the plot.
#'     \item If \code{output_table=FALSE}, invisibly returns \code{NULL}.
#'   }
#'
#' @details
#' The function expects two numeric columns. If the experiment has more than
#' two columns, the name of the columns of interest can be specificed by using
#' the parameter \code{column_name}. If \code{x} is a
#' \code{SummarizedExperiment}, it extracts the indicated assay and extracts
#' the columns of interest
#'
#' It also uses:
#'   - \code{centmass()} for computing \code{comx}.
#'   - \code{entropy()} for computing Shannon entropy, stored in the
#'   \code{entropy} column.  Between two variables, entropy rangeS between 0
#'   and 1.
#'
#' The rope is drawn in the middle of the plot (the x-axis from -1 to 1, y = 0),
#' with thickness \code{rope_width}. Points are scattered in \code{comy}
#' direction for a bit of jitter within the rope.
#'
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' library(airway)
#' data('airway')
#' se <- airway
#'
#' # Only use a random subset of 1000 rows
#' set.seed(123)
#' idx <- sample(seq_len(nrow(se)), size = min(1000, nrow(se)))
#' se <- se[idx, ]
#'
#' ## Normalize the data first using tpm_normalization
#'
#' rowData(se)$gene_length = rowData(se)$gene_seq_end
#' - rowData(se)$gene_seq_start
#'
#' se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
#'
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'
#' df <- assay(se, 'tpm_norm')
#' df <- as.data.frame(df)
#'
#' # Choose two columns of interest, in this case 'SRR1039508'
#' # and SRR1039516'
#'
#' # Default Behaviour
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE)
#'
#' # Colors can be modified
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('darkgreen', 'darkred'))
#'
#' # Emphasis can be applied to highly dominant variables by controling
#' # entropy parameter,
#' # values outside of that range will be colored smokewhite.
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('darkgreen', 'darkred'),
#'           entropyrange = c(0,0.1))
#'
#' # Points in the center are a reflection of genes with expression levels = 0.
#' # This can be modified by adjusting the maxvalue range
#'
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('darkgreen', 'darkred'),
#'           entropyrange = c(0,0.1),
#'           maxvaluerange = c(2, Inf))
#'
#' # By controling entropy range, you can observe different types of genes.
#' # Values closer to 0 represent dominance and closer to 1 shareness.
#'
#' # Exploring genes with high normalized expression values across different
#' #' entropy ranges
#'
#'
#'
#' # Looking for genes with a Log2(TPM) score between 4 and 8
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('darkgreen', 'darkred'),
#'           entropyrange = c(0,0.1),
#'           maxvaluerange = c(4, 8))
#'
#'
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('darkgreen', 'darkred'),
#'           entropyrange = c(0.1,0.8),
#'           maxvaluerange = c(4, 8))
#'
#'
#' plot_rope(df,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('darkgreen', 'darkred'),
#'           entropyrange = c(0.8,1),
#'           maxvaluerange = c(4, 8))
#'
#' # -------------------------------
#' # 1) Using a SummarizedExperiment
#' # -------------------------------
#'
#' plot_rope(se,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('lightgreen', 'indianred'),
#'           entropyrange = c(0,0.1),
#'           maxvaluerange = c(4, 8))
#'
#'
#' plot_rope(se,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col =c('lightgreen', 'indianred'),
#'           entropyrange = c(0.1,0.8),
#'           maxvaluerange = c(4, 8))
#'
#'
#' plot_rope(se,
#'           column_name = c("SRR1039508", "SRR1039516"),
#'           output_table = FALSE,
#'           col = c('lightgreen', 'indianred'),
#'           entropyrange = c(0.8,1),
#'           maxvaluerange = c(4, 8))
#'
#' ### Obtaining the DF output for the analysis
#'
#' object <- plot_rope(se,
#'                    column_name = c("SRR1039508", "SRR1039516"),
#'                    output_table = TRUE,
#'                    col = c('lightgreen', 'indianred'),
#'                    entropyrange = c(0.8,1),
#'                    maxvaluerange = c(4, 8))
#'
#' head(object)
#'


plot_rope <- function(x,
                    column_name = NULL,
                    push_text     = 1,
                    rope_width    = 1,
                    rope_color    = "#CCCCCCCC",
                    rope_border   = TRUE,
                    col           = c("red", "blue"),
                    col_bg        = "whitesmoke",
                    pch           = c(21, 21),
                    pch_bg        = 19,
                    cex           = 1,
                    entropyrange  = c(0, Inf),
                    maxvaluerange = c(0, Inf),
                    plotAll       = TRUE,
                    label         = TRUE,
                    output_table  = TRUE,
                    assay_name    = NULL)
{
    #-------------------------#
    # 1) Acquire / check data
    #-------------------------#
    if (inherits(x, "SummarizedExperiment")) {
        # SummarizedExperiment path
        if (is.null(assay_name)) {
            all_assays <- SummarizedExperiment::assayNames(x)
        if (length(all_assays) < 1)
            stop("No assays found in the SummarizedExperiment.")
        assay_name <- all_assays[1]
        }
        mat <- SummarizedExperiment::assay(x, assay_name)
        if (!is.matrix(mat) || !is.numeric(mat)) {
            stop("The selected assay must be a numeric matrix.")
        }
        if (length(column_name) == 2){
            mat <- mat[,column_name]
        }
        if (ncol(mat) != 2) {
            stop("plot_rope() requires exactly 2 columns of data; found ",
                ncol(mat))
        }
        data <- as.data.frame(mat)
        original_colnames <- colnames(mat)
    } else if (is.data.frame(x) || is.matrix(x)) {
        # Data.frame or matrix
        if (is.data.frame(x)) {
            mat <- as.data.frame(x)
        }
        if (!is.numeric(as.matrix(x))) {
            stop("Data is not numeric.")
        }
        if (length(column_name) == 2){
            mat <- mat[,column_name]
        }
        if (ncol(mat) != 2) {
            stop("plot_rope() requires exactly 2 columns of data; found ",
                ncol(mat))
        }
        data <- as.data.frame(mat)
        original_colnames <- colnames(data)
    } else {
        stop("Input must be a data.frame/matrix or SummarizedExperiment.")
    }

    rope_width <- rope_width * 0.25

    # Rename columns to "a","b" for convenience
    colnames(data) <- c("a","b")

    #-------------------------#
    # 2) Compute positions
    #-------------------------#
    # We'll "center of mass" on x-axis from -1..1 (like a seesaw),
    #  using px=c(-1,1), py=c(0,0)
    #  Then add random y-jitter inside rope.
    px <- c(-1, 1)
    py <- c(0, 0)

    cm_df <- centmass(data, x_coord=px, y_coord=py)
    cm_df$comy <- rope_width * sample(-500:500, size=nrow(cm_df),
                                    replace=TRUE) / 500

    final <- cbind(data, cm_df)


    final$color <- ifelse(final$a > final$b, col[1], col[2])

    #-------------------------#
    # 3) Filter: maxvalue, entropy
    #-------------------------#
    # maxvalue
    final$maxvalue <- pmax(final$a, final$b)

    ent_df <- entropy(data)
    final$entropy <- ent_df$Entropy

    # If outside entropyrange, color => col_bg
    outside_entropy <- !(final$entropy >= entropyrange[1] & final$entropy
                        <= entropyrange[2])
    final$color[outside_entropy] <- col_bg

    # If outside maxvaluerange, color => col_bg
    outside_maxval <- !(final$maxvalue >= maxvaluerange[1]
                        & final$maxvalue <= maxvaluerange[2])
    final$color[outside_maxval] <- col_bg

    #-------------------------#
    # 4) Base R Plot
    #-------------------------#
    # set up an empty plot
    plot(
        x     = 0,
        y     = 0,
        xlim  = px * push_text,  # e.g. from -1.5..1.5 if push_text=1
        pch   = 16,
        xlab  = "",
        ylab  = "",
        frame.plot = FALSE,
        ann   = TRUE,
        axes  = FALSE,
        col   = "transparent"
    )

    # draw rope polygon
    polygon(
        x = c(-1, -1, 1, 1),
        y = rope_width * c(-1, 1, 1, -1),
        col    = rope_color,
        border = 'transparent'
    )


    # 4a) Possibly plot background points
    if (plotAll) {
        # those with color = col_bg
        badpoints <- final[final$color == col_bg, ]
        points(
            badpoints$comx,
            badpoints$comy,
            pch = pch_bg,
            cex = cex,
            col = if (pch_bg %in% 21:25) "black" else badpoints$color,
            bg  = if (pch_bg %in% 21:25) badpoints$color else "black"
        )
    }

    # 4b) Plot highlight points
    goodpoints <- final[final$color != col_bg, ]
    type1 <- goodpoints[goodpoints$color == col[1], ]
    type2 <- goodpoints[goodpoints$color == col[2], ]

    # type1
    points(
        type1$comx,
        type1$comy,
        pch = pch[1],
        cex = cex,
        col = if (pch[1] %in% 21:25) "black" else type1$color,
        bg  = if (pch[1] %in% 21:25) type1$color else "black"
    )

    # type2
    points(
        type2$comx,
        type2$comy,
        pch = pch[2],
        cex = cex,
        col = if (pch[2] %in% 21:25) "black" else type2$color,
        bg  = if (pch[2] %in% 21:25) type2$color else "black"
    )

    # draw rope polygon
    polygon(
        x = c(-1, -1, 1, 1),
        y = rope_width * c(-1, 1, 1, -1),
        col    = 'transparent',
        border = rope_border
    )

    #-------------------------#
    # 5) Optional Labels
    #-------------------------#
    if (label) {
        # Label the original column names near the left & right of the rope
        text(original_colnames,
            x = px * push_text,  # e.g. -1.2.. +1.2
            y = py)
    }

    #-------------------------#
    # 6) Return
    #-------------------------#
    if (output_table) {
        return(final)
    } else {
        invisible(NULL)
    }
}
