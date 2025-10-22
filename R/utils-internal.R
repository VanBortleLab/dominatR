
.get_matrix <- function(se, a_name = NULL){
        an <- if(is.null(a_name)) assayNames(se)[1] else a_name
        if(is.na(an) || is.null(an))
                stop('No assays found in the SummarizedExperiment')
        mat <- SummarizedExperiment::assay(se, an)
        if(is.null(mat))
                stop('No assay named "', an, '" found')
        if(!is.numeric(mat))
        stop('Selected assay is not numeric')
        return(list(mat = mat, assay_name = an))
}

.assing_assay <- function(se, matrix,  a_name = NULL, new_a_name = NULL){
        stopifnot(inherits(se, 'SummarizedExperiment'))
        fin_name <- if (is.null(new_a_name)) a_name else new_a_name
        assay(se, fin_name) <- matrix

        return(se)
}

.get_matrix_df <- function(x){
        if (is.matrix(x)) {
        if (!is.numeric(x)) stop("Input must be numeric.")
        return(x)
        } else if (is.data.frame(x)) {
        dat <- as.matrix(x)
        if (!is.numeric(dat)) stop("Input must be numeric.")
        return(dat)
        } else {
        stop("Input must be a matrix or data.frame.")
        }
}


.get_gene_length_se <- function(se, col = 'gene_legth'){

        rd <- rowData(se)

        if(!(col %in% colnames(rd))) {
        stop("No '", col, "' column found in rowData(x). ",
        "Please add rowData(x)$", col, " or use data.frame/matrix mode.")
        }

        gene_len_vect <- rd[[col]]

        if (!is.numeric(gene_len_vect)) {
        stop("'gene_length' in rowData(x) must be numeric.")
        }

        return(gene_len_vect)
}

.eval_gene_length_df <- function(gene_length, df){
        if (is.null(gene_length)) {
        stop("You must provide 'gene_length' for data.frame/matrix input.")
        }
        if (!is.numeric(gene_length)) {
        stop("Argument 'gene_length' must be numeric.")
        }
        if (length(gene_length) != nrow(df)) {
        stop("Length of 'gene_length' (", length(gene_length),
        ") must match the number of rows (", nrow(df), ") in 'x'.")
        }

        return(gene_length)
}



.get_column_df <- function(mat, column_name, type = c(2, 3)){
        if (length(column_name) == type){
        mat <- mat[,column_name]
        }
        if (ncol(mat) != type) {
        stop("plot_rope() and plot_triangle()
        require exactly 2 or 3 columns of data respectively; found ",
        ncol(mat))
        }
        data <- as.data.frame(mat)
        return(data)
}

.arc_slice_polygon <- function(block, deg_sp, deg, type_label){

        start1 <- deg_sp + deg * (block - 1)
        end1 <- deg_sp + deg * (block)

        p <- ggplot() +
        geom_arc(aes(
        x0 = 0,
        y0 = 0,
        r = 100,
        start = start1,
        end = end1
        ))


        built <- ggplot_build(p)$data[[1]]
        arc_xy <- built[, c("x", "y")]

        poly <- rbind(
        data.frame(x = 0, y = 0),
        data.frame(x = arc_xy$x, y = arc_xy$y),
        data.frame(x = 0, y = 0)
        )

        poly$type <- rep(type_label, nrow(poly))

        return(poly)

}




