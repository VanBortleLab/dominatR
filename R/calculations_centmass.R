
#' Compute the "center of mass" for rows of a data frame or SummarizedExperiment
#'
#' @description
#' For each row of the numeric data, \code{centmass()} computes a 2D center of mass
#' with coordinates (\code{comx}, \code{comy}). The \code{x_coord} and \code{y_coord}
#' vectors specify the location for each column's "mass."
#'
#' The original usage assumes a ternary coordinate system by default, but this
#' can be generalized to any scenario where columns represent discrete "masses"
#' at known (x,y) positions.
#'
#' By default, \code{x_coord = c(0, 1, 0.5)} and
#' \code{y_coord = c(0, 0, sqrt(3)/2)}, which correspond to the corners
#' of an equilateral triangle (often used in ternary plots).
#'
#' @param x A data.frame (with numeric columns) or a SummarizedExperiment.
#' @param x_coord Numeric vector of length equal to the number of columns
#'   in \code{x}, specifying the x-coordinates of each column's mass.
#' @param y_coord Numeric vector of length equal to the number of columns
#'   in \code{x}, specifying the y-coordinates of each column's mass.
#' @param assay_name If \code{x} is a SummarizedExperiment, the name of the
#'   assay to use. Defaults to the first assay if not specified.
#'
#' @return
#' \itemize{
#'   \item If \code{x} is a data.frame, returns a new \code{data.frame} with
#'     columns \code{comx} and \code{comy}.
#'   \item If \code{x} is a SummarizedExperiment, returns the same object but
#'     with two new columns \code{comx} and \code{comy} in \code{rowData(x)}.
#'}
#'
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#'library(SummarizedExperiment)
#'library(airway)
#'data('airway')
#'
#'se = airway
#' # Let's subset for the first 3 columns for this example
#'se = se[,1:3]
#'# -------------------------------
#'# 1) Using a data.frame
#'# -------------------------------
#'
#'
#'df = assay(se) |> as.data.frame()
#'
#'df = centmass(df)
#'head(df)
#'
#'# -------------------------------
#'# 2) Using a SummarizedExperiment
#'# -------------------------------
#'
#'se2 = centmass(se)
#'
#'## X and Y coordinates are stored in rowData(se2)
#'head(rowData(se2))
#'
#'
centmass <- function(
    x,
    x_coord = c(0, 1, 0.5),
    y_coord = c(0, 0, sqrt(3)/2),
    assay_name = NULL
) {
  #----------------------------------#
  # Handle SummarizedExperiment case #
  #----------------------------------#
  if (inherits(x, "SummarizedExperiment")) {
    if (is.null(assay_name)) {
      # use the first assay by default
      all_assays <- assayNames(x)
      if (length(all_assays) < 1) {
        stop("No assays found in the SummarizedExperiment.")
      }
      assay_name <- all_assays[[1]]
    }

    mat <- assay(x, assay_name)
    if (is.null(mat)) {
      stop("No assay named '", assay_name, "' found in the SummarizedExperiment.")
    }

    # check length of x_coord, y_coord
    if (ncol(mat) != length(x_coord) || ncol(mat) != length(y_coord)) {
      stop("Length of x_coord or y_coord does not match the number of columns in the assay.")
    }

    # row sums
    row_sums <- rowSums(mat, na.rm = TRUE)
    # avoid division by 0
    row_sums <- ifelse(row_sums > 0, row_sums, 1)

    # compute weighted sums in x and y directions
    # vectorize by repeating x_coord across rows
    comx_vals <- rowSums(mat * matrix(rep(x_coord, each=nrow(mat)),
                                      nrow=nrow(mat)), na.rm=TRUE)
    comy_vals <- rowSums(mat * matrix(rep(y_coord, each=nrow(mat)),
                                      nrow=nrow(mat)), na.rm=TRUE)

    # divide by row sum
    comx_vals <- comx_vals / row_sums
    comy_vals <- comy_vals / row_sums

    # store in rowData
    SummarizedExperiment::rowData(x)$comx <- comx_vals
    SummarizedExperiment::rowData(x)$comy <- comy_vals

    return(x)

    #-----------------------#
    # Handle data.frame case
    #-----------------------#
  } else if (is.data.frame(x)) {
    # Convert entire data.frame to matrix (assuming it's numeric columns
    # or user has made sure it is).
    # If your real use-case has non-numeric columns that you need to skip,
    # you can subset them similarly to the entropy approach.
    mat <- as.matrix(x)

    # Dimension checks
    if (ncol(mat) != length(x_coord) || ncol(mat) != length(y_coord)) {
      stop("Length of x_coord or y_coord does not match the number of columns in 'x'.")
    }

    row_sums <- rowSums(mat, na.rm = TRUE)
    row_sums <- ifelse(row_sums > 0, row_sums, 1)

    comx_vals <- rowSums(mat * matrix(rep(x_coord, each=nrow(mat)), nrow=nrow(mat)), na.rm=TRUE)
    comy_vals <- rowSums(mat * matrix(rep(y_coord, each=nrow(mat)), nrow=nrow(mat)), na.rm=TRUE)

    # divide by row sum
    comx_vals <- comx_vals / row_sums
    comy_vals <- comy_vals / row_sums

    # return a small data.frame with comx, comy
    return(data.frame(comx = comx_vals, comy = comy_vals))

  } else {
    stop("Input must be either a data.frame or a SummarizedExperiment.")
  }
}
