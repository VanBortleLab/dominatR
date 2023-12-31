% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{rnapol_score}
\alias{rnapol_score}
\title{RNA Polymerase score dataframe}
\format{
A dataframe with 1061 rows and 7 variables:
\describe{
\item{Chr}{Chromosome where the gene is located}
\item{Start}{Gene start coordinate}
\item{Stop}{Gene stop coordinate}
\item{RNA_Type}{Type of transcript}
\item{pol1}{Score for RNA Polymerase I Chip-Seq}
\item{pol2}{Score for RNA Polymerase II Chip-Seq}
\item{pol3}{Score for RNA Polymerase III Chip-Seq}}
}
\source{
{Created by The VanBortle lab at UIUC to serve as an example}
}
\usage{
rnapol_score
}
\description{
A subset of data obtained from CHIP-Seq Data for RNA Polymerase I, II and III. Values are represented as the -log10(p-value)
}
\keyword{datasets}
