#' RNA Polymerase raw counts dataframe
#'
#' A subset of data obtained from CHIP-Seq Data for RNA Polymerase I, II and III
#'
#' @docType data
#' @format A dataframe with 1061 rows and 7 variables:
#' \describe{
#' \item{Chr}{Chromosome where the gene is located}
#' \item{Start}{Gene start coordinate}
#' \item{Stop}{Gene stop coordinate}
#' \item{RNA_Type}{Type of transcript}
#' \item{pol1}{Raw counts for RNA Polymerase I Chip-Seq}
#' \item{pol2}{Raw counts  for RNA Polymerase II Chip-Seq}
#' \item{pol3}{Raw counts for RNA Polymerase III Chip-Seq}}
#'
#' @source  {Created by The VanBortle lab at UIUC to serve as an example}
#' @usage data(rnapol_counts)
#' @examples
#' data('rnapol_counts')
#' head(rnapol_counts)
#'
#'
'rnapol_counts'


#' RNA Polymerase score dataframe
#'
#' A subset of data obtained from CHIP-Seq Data for RNA Polymerase I, II and III. Values are represented as the -log10(p-value)
#'
#' @docType data
#' @format A dataframe with 1061 rows and 7 variables:
#' \describe{
#' \item{Chr}{Chromosome where the gene is located}
#' \item{Start}{Gene start coordinate}
#' \item{Stop}{Gene stop coordinate}
#' \item{RNA_Type}{Type of transcript}
#' \item{pol1}{Score for RNA Polymerase I Chip-Seq}
#' \item{pol2}{Score for RNA Polymerase II Chip-Seq}
#' \item{pol3}{Score for RNA Polymerase III Chip-Seq}}
#'
#' @source  {Created by The VanBortle lab at UIUC to serve as an example}
#' @usage data(rnapol_score)
#' @examples
#' data('rnapol_score')
#' head(rnapol_score)
'rnapol_score'




#' ATAC-Seq rawcounts for POL3 Genes Dataframe
#'
#' A list of tissues and the corresponding counts for RNA POL3 genes for each of them.
#'
#' @docType data
#' @format A dataframe with 9817 rows and 26 variables:
#' \describe{
#' \item{core_type}{Category for genes based on their expression across tissues}
#' \item{Chr}{Chromosome where the RNA of interest is located}
#' \item{Start}{Gene start coordinate}
#' \item{End}{Gene end coordinate}
#' \item{Gene}{Gene name}
#' \item{Index}{Unique RNA Sequence identifier - Can be retrieved from RNA Central}
#' \item{Type}{Type of RNA POL3 Transcript}
#' \item{Tissue}{Remaining columns contain the name of different assessed tissues}}
#'
#'
#' @source  {Created by The VanBortle lab at UIUC to serve as an example}
#' @usage data(atac_tissue_counts)
#' @examples
#' data('atac_tissue_counts')
#' head(atac_tissue_counts)
'atac_tissue_counts'

#' ATAC-Seq Score for POL3 Genes Dataframe
#'
#' A list of tissues and the corresponding counts for RNA POL3 genes for each of them. Values are represented as the -log10(p-value)
#'
#' @docType data
#' @format A dataframe with 9817 rows and 23 variables:
#' \describe{
#' \item{core_type}{Category for genes based on their expression across tissues}
#' \item{Chr}{Chromosome where the RNA of interest is located}
#' \item{Start}{Gene start coordinate}
#' \item{End}{Gene end coordinate}
#' \item{Gene}{Gene name}
#' \item{Index}{Unique RNA Sequence identifier - Can be retrieved from RNA Central}
#' \item{Type}{Type of RNA POL3 Transcript}
#' \item{Tissue}{Remaining columns contain the name of different assessed tissues}}
#'
#'
#' @source  {Created by The VanBortle lab at UIUC to serve as an example}
#' @usage data(atac_tissue_counts)
#' @examples
#' data('atac_tissue_score')
#' head(atac_tissue_score)
'atac_tissue_score'
