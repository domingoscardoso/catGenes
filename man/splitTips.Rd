% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitTips.R
\name{splitTips}
\alias{splitTips}
\title{Parse tip labels into a dataframe with multiple columns for each content}
\usage{
splitTips(tiplabels = NULL)
}
\arguments{
\item{tiplabels}{A vector of tip labels from a phylogeny.}
}
\value{
A dataframe.
}
\description{
This function is designed to split tip labels into a dataframe
object. Each column in the dataframe corresponds to the genus, species,
infraspecific name, and associated voucher and GenBank information, whenever
available. If doubtful particles such as 'aff.' or 'cf.' are provided, they
will also be separated into a distinct column.
}
\examples{
\dontrun{
library(catGenes)

data(Harpalyce_bayes_tree)

df <- splitTips(tiplabels = Harpalyce_bayes_tree@phylo$tip.label)
}

}
\author{
Domingos Cardoso
}
