#' The catGenes package is under development and is intended to help researchers
#' in phylogenetics and phylogenomics to build a fully concatenated or combined
#' (non-interleaved) dataset by automatically comparing individual DNA alignments.
#' Other features also involve mining DNA sequences or targeted loci from
#' plastomes and mitochondrial genomes from GenBank, as well as plotting
#' phylogenetic figures.
#'
#' The package's main functions \code{\link{catfullGenes}} and \code{\link{catmultGenes}}
#' compare at least two individual DNA alignments with differing (or equal) number
#' of sequences and creates equally-sized DNA alignments that are useful for
#' concatenated phylogenetic or phylogenomic analyses. By using the functions
#' \code{\link{writeNexus}} or \code{\link{writePhylip}} the user can write the
#' catfullGenes- or catmultGenes-derived list of DNA alignments in NEXUS or PHYLIP
#' format, respectively, which can serve for downstream model-based phylogenetic
#' analyses in MrBayes, BEAST, RAxML or PAUP.
#'
#' Another set of functions in this package allow the user to write a \code{data.frame}
#' formatted DNA alignment (two-column-sized table including the taxon names and
#' corresponding DNA sequence) as NEXUS, PHYLIP or FASTA formats: \code{\link{nexusdframe}},
#' \code{\link{phylipdframe}}, \code{\link{fastadframe}}. These are useful for
#' readily writing each gene dataset from within the resulting list of compared
#' gene datasets, after running the functions \code{\link{catfullGenes}},
#' \code{\link{catmultGenes}} or \code{\link{dropSeq}}.
#'
#' For the most recent version of the catGenes, you are directed to
#' package's page on github (\url{http://www.github.com/domingoscardoso/catGenes}).
#'
#' @name catGenes-package
#'
#' @aliases catGenes-package
#'
#' @title Tools for combining individual DNA alignments for multigene analysis in
#' phylogenetics and phylogenomics, and other phylogenetic features
#'
#' @author \strong{Domingos Cardoso}\cr
#' (email: \email{cardosobot@@gmail.com};
#' Website: \url{https://biologia.ufba.br/domingos-benicio-oliveira-silva-cardoso})
#'
#' @author \strong{Quezia Cavalcante}\cr
#' (email: \email{queziacs@@yahoo.com.br})
#'
#' @keywords package
#'
#' @details \tabular{ll}{
#' Package: \tab catGenes\cr
#' Type: \tab Package\cr
#' Version: \tab 0.10\cr
#' Date: \tab 2020-10-25\cr
#' }
#'
#' @references Cardoso, D. & Cavalcante, Q. (2024). catGenes: a new R package for
#' combining multiple DNA alignments for multigene analysis in phylogenetics and
#' phylogenomics.
#'
#' @import dplyr magrittr tidyr stringr R.utils rmarkdown ape shiny shinydashboard shinyjs tibble stats genbankr ggtree ggtext ggplot2 cowplot phangorn phytools flora glue
#'
"_PACKAGE"

NULL
