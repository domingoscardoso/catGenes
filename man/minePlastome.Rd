% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minePlastome.R
\name{minePlastome}
\alias{minePlastome}
\title{Read and download targeted loci from plastome sequences in GenBank}
\usage{
minePlastome(genbank = NULL,
             taxon = NULL,
             voucher = NULL,
             CDS = TRUE,
             genes = NULL,
             verbose = TRUE,
             dir = "RESULTS_minePlastome")
}
\arguments{
\item{genbank}{A vector comprising the GenBank accession numbers specifically
corresponding to the plastome sequence targeted for locus mining.}

\item{taxon}{A vector containing the taxon name linked to the plastome sequence.
In the absence of this information, the function will default to the existing
nomenclature linked to the plastome, as originally provided in GenBank.}

\item{voucher}{A vector containing relevant voucher information linked to the
plastome sequence. If this information is supplied, the function will promptly
append it immediately following the taxon name of the downloaded targeted
sequence.}

\item{CDS}{a logical controlling whether the targeted loci are protein coding
genes, otherwise the function understands that entered gene names are e.g.
intron or intergenic spacer regions.}

\item{genes}{A vector of one or more gene names as annotated in GenBank.}

\item{verbose}{Logical, if \code{FALSE}, a message showing each step during
the GenBank search will not be printed in the console in full.}

\item{dir}{The path to the directory where the mined DNA sequences
in a fasta format file will be saved. The default is to create a directory
named \strong{RESULTS_minePlastome} and the sequences will be saved within a
subfolder named after the current date.}
}
\description{
A function built on the \href{https://bioconductor.org/packages/release/bioc/html/genbankr.html}{genbankr}
package, designed to establish a connection with the GenBank database. This
function reads plastome sequences using provided accession numbers, extracting
and formatting any specified targeted loci, and finally writing them in a
fasta file format.
}
\examples{
\dontrun{
library(catGenes)
library(dplyr)

data(GenBank_accessions)

GenBank_plastomes <- GenBank_accessions \%>\%
  filter(!is.na(Plastome)) \%>\%
  select(c("Species", "Voucher", "Plastome"))

minePlastome(genbank = GenBank_plastomes$Plastome,
             taxon = GenBank_plastomes$Species,
             voucher = GenBank_plastomes$Voucher,
             CDS = TRUE,
             genes = c("matK", "rbcL"),
             verbose = TRUE,
             dir = "RESULTS_minePlastome")
}

}
\author{
Domingos Cardoso
}
