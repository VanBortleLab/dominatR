#' RNA Polymerase CHIP-Seq Data
#'
#' A subset of data obtained from CHIP-Seq Data for RNA Polymerase I, II and III
#'
#' @format A dataframe with 1635 rows and 4 variables:
#' \describe{
#' \item{UU}{Type of transcript}
#' \item{pol1}{Normalized peak call value for RNA Polymerase I Chip-Seq}
#' \item{pol2}{Normalized peak call value for RNA Polymerase II Chip-Seq}
#' \item{pol3}{Normalized peak call value for RNA Polymerase III Chip-Seq}}
#'
#' @source  {Created by The VanBortle lab at UIUC to serve as an example}
#' @eval data(rnapol_chip)
'rnapol_chip'


#' ATAC-Seq Significance for POL3 Active Genes
#'
#' A subset of tissues with a list of POL3 Genes that are active.
#'
#' @format A dataframe with 1415 rows and 19 variables:
#' \describe{
#' \item{Chr}{Coordinates for RNA of interest}
#' \item{Gene}{Chromosome where the RNA of interest is located}
#' \item{Index}{Unique RNA Sequence identifier - Can be retrieved from RNA Central}}
#'
#'
#' @source  {Created by The VanBortle lab at UIUC to serve as an example}
#' @eval data(ATAC_Seq_Tissue)
'ATAC_Seq_Tissue'
