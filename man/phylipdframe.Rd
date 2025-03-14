% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylipdframe.R
\name{phylipdframe}
\alias{phylipdframe}
\title{Writes a PHYLIP-formatted DNA alignment from a dataframe-formatted DNA alignment}
\usage{
phylipdframe(x, file,
             dropmisseq = TRUE,
             endgaps.to.miss = TRUE)
}
\arguments{
\item{x}{The object to be written, any two-column-sized \code{data.frame} where
the first column contains the taxon names and the second column the DNA sequence.
Otherwise, the object may be a list-formatted NEXUS file as originally
imported with \code{\link{ape}}'s function \code{\link{read.nexus.data}}.}

\item{file}{Either a character string naming a file or a \code{\link{connection}}
open for writing.}

\item{dropmisseq}{Logical, if \code{FALSE} the function will not drop species
with empty DNA sequence. After running the concatenating function
\code{\link{catmultGenes}} using missdata = \code{TRUE}, and then using
\code{\link{dropSeq}} to remove duplicated accessions of the same species,
you might find useful to keep dropmisseq = \code{TRUE} so as to save each
individual DNA alignment by also removing species that fully miss the sequence
data.}

\item{endgaps.to.miss}{Logical, if \code{FALSE} the function will not replace
terminal GAPs into missing character (?).}
}
\description{
Writes \code{data.frame} formatted DNA alignment or \code{list}
formatted NEXUS file as originally imported with \code{\link{ape}}'s function
\code{\link{read.nexus.data}} into a PHYLIP-formatted file. It is useful for
writing each gene dataset from within the resulting list of compared gene datasets,
after running the concatenating functions \code{\link{catfullGenes}} and
\code{\link{catmultGenes}}. The function is also useful for saving into PHYLIP
format the original list-formatted NEXUS object as read by \code{\link{read.nexus.data}},
after making specific changes in such original individual alignment (e.g. corrections
of species names).
}
\examples{
\dontrun{
data(Gaya)
catdf <- catfullGenes(Gaya,
                      multiaccessions = FALSE,
                      shortaxlabel = TRUE,
                      missdata = FALSE,
                      outgroup = "Abutilon_costicalyx")

ITS <- catdf[[1]]
petLpsbE <- catdf[[2]]
rpl16 <- catdf[[3]]

phylipdframe(ITS, file = "filename.nex")
phylipdframe(petLpsbE, file = "filename.nex")
phylipdframe(rpl16, file = "filename.nex")
}

}
\seealso{
\code{\link{catfullGenes}}

\code{\link{catmultGenes}}

\code{\link{dropSeq}}
}
\author{
Domingos Cardoso
}
