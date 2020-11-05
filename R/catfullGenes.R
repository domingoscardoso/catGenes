#' Compares gene datasets for combined phylogenetic analysis
#'
#' @author Domingos Cardoso and Quezia Cavalcante
#'
#' @description Compares a list of "n" gene datasets (individual DNA alignments)
#' and makes them with the same number of taxa, ready for combined, multigene
#' phylogenetic analysis. This function only works for DNA alignments where species
#' do not have duplicated sequences. For concatenation when species have multiple
#' accessions, please use \code{\link{catmultGenes}}.
#'
#' @usage
#' catfullGenes(\dots,
#'              shortaxlabel = TRUE,
#'              missdata = TRUE,
#'              outgroup = NULL)
#'
#' @param ... a list of nexus-formatted gene datasets as read by ape's \code{\link{read.nexus.data}}
#' or at least two individually ape-read objects of nexus-formatted gene datasets.
#'
#' @param shortaxlabel Logical, if \code{FALSE} the final individual gene dataset will maintain
#' the accession numbers associated with each species or sequence.
#'
#' @param missdata Logical, if \code{FALSE} the comparison will exclude any species
#' that lacks a complete sequence for one of the input gene dataset.
#'
#' @param outgroup Provide the outgroup taxa (either one taxon name or a vector of
#' multiple taxon names that are present in all individual gene dataset) if the
#' concatenation is intended to maintain incomplete taxa (taxa missing the sequence for a particular gene).
#'
#' @return A list of dataframes of the equally-sized gene dataset, where the first column "species"
#' include all taxon names and the second column "sequence" include the DNA sequence for the corresponding taxon.
#'
#' @seealso \code{\link{writeNexus}}
#' @seealso \code{\link{writePhylip}}
#' @seealso \code{\link{nexusdframe}}
#' @seealso \code{\link{phylipdframe}}
#' @seealso \code{\link{fastadframe}}
#'
#' @examples \dontrun{
#' data(Gaya)
#' catdf <- catfullGenes(Gaya,
#'                       shortaxlabel = TRUE,
#'                       missdata = FALSE,
#'                       outgroup = "Abutilon_costicalyx")
#'
#' outgrouptaxa <- c("Abutilon_costicalyx", "Abutilon_itatiaie")
#' catdf <- catfullGenes(Gaya,
#'                       shortaxlabel = TRUE,
#'                       missdata = FALSE,
#'                       outgroup = outgrouptaxa)
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom stats setNames na.omit
#'
#' @export

catfullGenes <- function(...,
                         shortaxlabel = TRUE,
                         missdata = TRUE,
                         outgroup = NULL) {

  # Loading all individual genes into a single named list
  datset <- .namedlist(...)

  if(length(datset) == 1){
    datset <- datset[[1]]
  }

  numberdatset <- length(datset)

  if (numberdatset == 1 | numberdatset == 0) {
    stop("You must provide at least TWO gene datasets in the following format:
         datset=gene1, gene2, gene3... or a list of genes in a single vector
         Find help also at DBOSLab-UFBA
         (Domingos Cardoso; cardosobot@gmail.com)")
  }


  # Adjusting species labels when they have cf or aff
  adjust_cf <- lapply(datset, function(x) grepl("_cf_", names(x)))
  adjust_aff <- lapply(datset, function(x) grepl("_aff_", names(x)))

  if(any(unlist(adjust_cf))){
    cat("Any gene dataset includes species under \"cf.\"", sep = "\n")
    # Finding species labels to rename
    spp_to_rename <- vector()
    for(i in seq_along(datset)){
       if(length(names(datset[[i]])[adjust_cf[[i]]]) > 0){
         spp_to_rename[[i]] <- names(datset[[i]])[adjust_cf[[i]]]
         spp_to_rename <- spp_to_rename[!is.na(spp_to_rename)]
      }
    }
    cat(spp_to_rename, "WERE RENAMED TO", gsub("_cf_", "_cf", spp_to_rename), "", sep = "\n")
    spp_labels <- lapply(datset, function(x) gsub("_cf_", "_cf", names(x)))
    for(i in seq_along(datset)){
      names(datset[[i]]) <- spp_labels[[i]]
    }
  }

  if(any(unlist(adjust_aff))){
    cat("Any gene dataset includes species under \"aff.\"", sep = "\n")
    # Finding species labels to rename
    spp_to_rename <- vector()
    for(i in seq_along(datset)){
      if(length(names(datset[[i]])[adjust_aff[[i]]]) > 0){
        spp_to_rename[[i]] <- names(datset[[i]])[adjust_aff[[i]]]
        spp_to_rename <- spp_to_rename[!is.na(spp_to_rename)]
      }
    }
    cat(spp_to_rename, "WERE RENAMED TO", gsub("_aff_", "_aff", spp_to_rename), "", sep = "\n")
    spp_labels <- lapply(datset, function(x) gsub("_aff_", "_aff", names(x)))
    for(i in seq_along(datset)){
      names(datset[[i]]) <- spp_labels[[i]]
    }
  }

  # Now running genecomp function in a for loop

  cat(cat("Matching first the gene", names(datset[1]), "with:",
          paste0(names(datset[-1]), "...", collapse = " ")), "", sep = "\n")

  # Creating an empty list to fill in during the loop iteration
  datsetcomp <- list()
  for (i in 1:(numberdatset-1)) {

    if (numberdatset == 2) {
      cat(cat("Gene comparison will exclude sequence set from", names(datset[1]),
              "that is not in", paste0(names(datset[i+1]), "...", collapse = " ")), "",
          sep = "\n")
    } else {
      cat(cat("Gene comparison will exclude sequence set from", names(datset[1]),
              "that is not in",
              paste0(names(datset[i+1]), "...", collapse = " ")), "", sep = "\n")
    }
    # Looping over genes
    datsetcomp[[1]] <- .genecomp(datset[[1]], datset[[i+1]],
                                 data = datset,
                                 loop = i,
                                 shortaxlabel = shortaxlabel,
                                 missdata = missdata,
                                 outgroup = outgroup)
    datset[[1]] <- datsetcomp[[1]]
  }

  cat(cat("Matched result of gene", names(datset[1]),
          "is again matched with",
          paste0(names(datset[-1]), "...", collapse = " ")), "", sep = "\n")

  for (i in 2:numberdatset) {

    if (numberdatset == 2) {
      cat(cat("Gene comparison will exclude sequence set from", names(datset[i]),
              "that is not in", paste0(names(datset[1]), "...", collapse=" ")), "", sep = "\n")
    } else {
      cat(cat("Gene comparison will exclude sequence set from", names(datset[i]),
              "that is not in the matched result of",
              paste0(names(datset[1]), "...", collapse=" ")), "", sep = "\n")
    }
    # Looping over genes
    datset[[i]] <- .genecomp(datset[[i]], datset[[1]],
                             data = datset,
                             loop = i-1,
                             shortaxlabel = shortaxlabel,
                             missdata = missdata,
                             outgroup = outgroup)
  }
  datsetcomp <- datset

  cat("Full gene match is finished!", "",
      sep="\n")

  return(datsetcomp)
}
