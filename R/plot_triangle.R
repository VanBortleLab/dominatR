#' Triangle (ternary) dominance plot
#'
#' @description
#' Creates a triangular (ternary) scatter plot for **three** numeric variables
#' Each point is coloured by the variable with the largest value and can be
#' filtered by (i) Entropy score ranging from (0 to 1.585) and
#' (ii) overall score
#'
#' The plot is useful for visualising “winner-takes-all” behaviour in three-way
#' comparisons, e.g. gene expression in *A*, *B*, *C* conditions.
#'
#' @param x            A numeric \code{data.frame}/\code{matrix} **or**
#'   a \code{SummarizedExperiment}.
#' @param column_name  Character. Names (or indices) of the three columns to
#'   visualise. If \code{NULL}, the first three numeric columns are used.
#' @param entropyrange Numeric. Keep points whose entropy lies inside
#'   this interval. Default is \code{c(0,Inf)}
#' @param maxvaluerange Numeric. Keep points whose values lies inside
#' this interval. Default is \code{c(0,Inf)}
#' @param col          Character. Colors for each variable.
#' @param background_col Character. Color for the observations outside
#' \code{entropyrange} and \code{maxvaluerange}
#' @param output_table Logical. If \code{TRUE} returns the processed data frame.
#' @param plotAll      Logical. If \code{TRUE}, filtered points are shown in
#'   \code{background_col}; if \code{FALSE}, they are omitted.
#' @param cex,pch      Base-graphics point size / symbol.
#' @param assay_name   (SummarizedExperiment only) Which assay to use. Default:
#'   the first assay.
#'
#' @return
#' If \code{output_table = TRUE}, a \code{data.frame} with the original three
#' columns plus:
#' \itemize{
#'   \item \code{comx}, \code{comy}  — Cartesian coordinates in the triangle;
#'   \item \code{color}              — final plotting colour;
#'   \item \code{entropy}            — Entropy scores for each gene;
#'   \item \code{max_counts}         — Maximum score across variables
#' }
#'
#' @details
#' The function expects three numeric columns. If the experiment has more than
#' three columns, the name of the columns of interest can be specified by using
#'  the parameter \code{column_name}. If \code{x} is
#'  a \code{SummarizedExperiment}, it extracts the indicated assay and extracts
#'   the columns of interest
#'
#' It also uses:
#'   - \code{centmass()} for computing \code{comx} and \code{comy}.
#'   - \code{entropy()} for computing Shannon entropy, stored in the
#'     \code{entropy} column.  Between three variables, entropy rangeS between
#'     0 and 1.585.
#'
#' The ternary vertices are fixed at
#' \eqn{( \sin(0),  \cos(0) )},
#' \eqn{( \sin(2\pi/3),  \cos(2\pi/3) )} and
#' \eqn{( \sin(4\pi/3),  \cos(4\pi/3) )}.
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
#' ## Normalize the data first using tpm_normalization
#' rowData(se)$gene_length <- rowData(se)$gene_seq_end -
#'                            rowData(se)$gene_seq_start
#'
#' se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
#'
#' # -------------------------------
#' # 1) Using a data.frame
#' # -------------------------------
#'
#' df <- assay(se, 'tpm_norm') |> as.data.frame()
#'
#'
#' # Choose three columns of interest, in this case 'SRR1039508', 'SRR1039516'
#' and 'SRR1039512'
#'
#' # Default Behaviour
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE)
#'
#' # Colors can be modified
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'))
#'
#' # Emphasis can be applied to highly dominant variables by controling
#' entropy parameter,
#' # values outside of that range will be colored smokewhite.
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(0, 0.1))
#'
#' # Points in the center are a reflection of genes with expression levels = 0.
#' # This can be modified by adjusting the maxvalue range
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(0, 0.1),
#'               maxvaluerange = c(0.1, Inf))
#'
#' # By controling entropy range, you can observe different types of genes.
#' # Values closer to 0 represent dominance and closer to 1.6 shareness.
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(0, 0.4),
#'               maxvaluerange = c(0.1, Inf))
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(0.4, 1.3),
#'               maxvaluerange = c(0.1, Inf))
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(1.3, Inf),
#'               maxvaluerange = c(0.1, Inf))
#'
#' # Same analysis can be performed by filtering out genes with low expression
#' values
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(1.2, Inf),
#'               maxvaluerange = c(2, Inf))
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(1.2, Inf),
#'               maxvaluerange = c(5, Inf))
#'
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(1.2, Inf),
#'               maxvaluerange = c(10, Inf))
#'
#' # Background points can be removed
#' plot_triangle(x = df,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('indianred', 'lightgreen', 'lightblue'),
#'               entropyrange = c(1.2, Inf),
#'               maxvaluerange = c(2, Inf),
#'               plotAll = FALSE)
#' # -------------------------------
#' # 1) Using a SummarizedExperiment
#' # -------------------------------
#'
#'
#' plot_triangle(x = se,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('darkred', 'darkgreen', 'darkblue'),
#'               entropyrange = c(0, 0.4),
#'               maxvaluerange = c(0.1, Inf),
#'               assay_name = 'tpm_norm')
#'
#' plot_triangle(x = se,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('darkred', 'darkgreen', 'darkblue'),
#'               entropyrange = c(0.4, 1.3),
#'               maxvaluerange = c(0.1, Inf),
#'               assay_name = 'tpm_norm')
#'
#' plot_triangle(x = se,
#'               column_name = c("SRR1039508", "SRR1039516", 'SRR1039512'),
#'               output_table = FALSE,
#'               col = c('darkred', 'darkgreen', 'darkblue'),
#'               entropyrange = c(1.3, Inf),
#'               maxvaluerange = c(0.1, Inf),
#'               assay_name = 'tpm_norm')
#'
#'
#' ### Obtaining the DF output for the analysis
#'
#' object = plot_triangle(x = se,
#'                        column_name = c("SRR1039508", "SRR1039516",
#'                        'SRR1039512'),
#'                        output_table = TRUE,
#'                        col = c('darkred', 'darkgreen', 'darkblue'),
#'                        entropyrange = c(1.3, Inf),
#'                        maxvaluerange = c(0.1, Inf),
#'                        assay_name = 'tpm_norm')
#'
#' head(object)
#'
#' @export

plot_triangle <- function(x,
                    column_name     = NULL,
                    entropyrange    = c(0, Inf),
                    maxvaluerange  = c(0, Inf),
                    col             = c("darkred", "darkgreen", "darkblue"),
                    background_col =  'whitesmoke',
                    output_table    = TRUE,
                    plotAll         = TRUE,
                    cex             = 1,
                    pch             = 16,
                    assay_name      = NULL)
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

    ### Extracting the columns of interest. Default the first three columns
    if (is.null(column_name)) column_name <- colnames(mat)[seq_len(3)]

    if (length(column_name) != 3)
        stop("Exactly three columns are required.")
    mat <- mat[, column_name, drop = FALSE]
    if (ncol(mat) != 3 || !is.numeric(mat))
        stop("Chosen columns must be numeric.")

    data <- as.data.frame(mat)
    original_colnames <- colnames(mat)

    } else if (is.data.frame(x) || is.matrix(x)) {
        # Data.frame or matrix
        if (is.data.frame(x)) {
            mat <- as.data.frame(x)
        }
    # Data.frame or matrix
    if (is.matrix(x)) {
        mat <- as.matrix(x)
    }
    if (!is.numeric(as.matrix(x))) {
        stop("Data is not numeric.")
    }
    if (length(column_name) == 3){
        mat <- mat[,column_name]
    }
    if (ncol(mat) != 3) {
        stop("plot_triangle() requires exactly 2 columns of data; found ",
            ncol(mat))
    }
    data <- as.data.frame(mat)
    original_colnames <- colnames(data)
    } else {
        stop("Input must be a data.frame/matrix or SummarizedExperiment.")
    }

    colnames(data) <- c("a", "b", "c")

    #-------------------------#
    # 2) Compute positions
    #-------------------------#
    verts <- data.frame(x = sin(2*pi*seq(0,2)/3),
                        y = cos(2*pi*seq(0,2)/3))
    cm     <- centmass(data, x_coord = verts$x, y_coord = verts$y)


    #-------------------------#
    # 3) Entropy and Max calculations
    #-------------------------#
    temp <- data |> mutate(max_counts = pmax(a, b, c)) |>
                    select(max_counts)



    data <- entropy(data)
    data <- cbind(data, cm, temp)

    rm(temp)

    data <- data |> mutate(max = pmax(a,b,c)) |>
        relocate(max, max_counts, comx, comy)


    #-------------------------#
    # 4) Colors and Dominance Assignments
    #-------------------------#

    data <- data |>
        mutate(color = ifelse( a == max, col[1],
                            ifelse(b == max, col[2], col[3])))

    keep <- data$Entropy >= entropyrange[1] &
        data$Entropy <= entropyrange[2] &
        data$max_counts >= maxvaluerange[1] &
        data$max_counts <= maxvaluerange[2]

    data$color[!keep] <- 'whitesmoke'

    #-------------------------#
    # 5) Plotting
    #-------------------------#
    par(pty = "s")
    plot(
        c(1, -1),
        c(1, -1),
        type = "n",
        axes = FALSE,
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        ann = FALSE
    )



    if (plotAll) {
        idx <- data$color == "whitesmoke"
        graphics::points(data$comx[idx],
                        data$comy[idx],
                        col = "whitesmoke",
                        pch = pch,
                        cex = cex)
    }
    idx <- data$color != "whitesmoke"
    graphics::points(data$comx[idx],
                    data$comy[idx],
                    col = data$color[idx],
                    pch = pch,
                    cex = cex)

    graphics::polygon(verts$x,
                    verts$y,
                    pch = 16)

    data <- data |> select(-max)

    if (output_table) return(data)
}
