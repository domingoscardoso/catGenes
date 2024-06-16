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
#'              outgroup = NULL,
#'              verbose = TRUE)
#'
#' @param ... a list of NEXUS-formatted gene datasets as read by ape's \code{\link{read.nexus.data}}
#' or at least two individually ape-read objects of NEXUS-formatted gene datasets.
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
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' gene matching search will not be printed in the console in full.
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
#'                       outgroup = "Abutilon_costicalyx",
#'                       verbose = TRUE)
#'
#' outgrouptaxa <- c("Abutilon_costicalyx", "Abutilon_itatiaie")
#' catdf <- catfullGenes(Gaya,
#'                       shortaxlabel = TRUE,
#'                       missdata = FALSE,
#'                       outgroup = outgrouptaxa,
#'                       verbose = TRUE)
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom stats setNames na.omit
#'
#' @export

catfullGenes <- function(...,
                         shortaxlabel = TRUE,
                         missdata = TRUE,
                         outgroup = NULL,
                         verbose = TRUE) {

  # Loading all individual genes into a single named list ####
  datset <- .namedlist(...)

  if (length(datset) == 1) {
    datset <- datset[[1]]
  }

  numberdatset <- length(datset)

  if (numberdatset == 1 | numberdatset == 0) {
    stop(paste0("You must provide at least TWO gene datasets in the following format:\n",
                "datset = gene1, gene2, gene3...\n",
                "or a list of genes in a single vector.\n\n"),
         "Find help also at DBOSLab-UFBA\n",
         "(Domingos Cardoso; cardosobot@gmail.com)")
  }


  cf <- lapply(datset, function(x) grepl("_cf_", names(x)))
  aff <- lapply(datset, function(x) grepl("_aff_", names(x)))
  spp_temp <- lapply(datset, function(x) gsub("_aff_|_cf_", " ", names(x)))
  infraspp <- lapply(spp_temp, function(x) grepl("[[:upper:]][[:lower:]]+_[[:lower:]]+_[[:lower:]]+",
                                                 x))

  if (any(unlist(cf))|any(unlist(aff))|any(unlist(infraspp))) {

    nr <- .namesTorename(datset,
                         cf = cf,
                         aff = aff,
                         infraspp = infraspp)

    # Adjusting species labels when they have cf or aff ####
    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    datset <- .adjustnames(datset,
                           cf = cf,
                           aff = aff,
                           infraspp = infraspp)
  }

  if (missdata == FALSE) {
    datset_temp <- datset
  }

  # Now running genecomp function in a for loop ####
  if (verbose) {
    message(paste0("Matching first the gene ", names(datset[1]), " with:\n",
                   paste0(names(datset[-1]), "...", collapse = " ")), "\n")
  }
  # Creating an empty list to fill in during the loop iteration ####
  datsetcomp <- list()
  for (i in 1:(numberdatset-1)) {
    if (verbose) {
      if (numberdatset == 2) {
        message(paste0("Gene comparison will exclude sequence set from ",
                       names(datset[1]), " that is not in ",
                       paste0(names(datset[i+1]), "...", collapse = " ")), "\n")
      } else {
        message(paste0("Gene comparison will exclude sequence set from ",
                       names(datset[1]), " that is not in ",
                       paste0(names(datset[i+1]), "...", collapse = " ")), "\n")
      }
    }
    # Looping over genes  ####
    datsetcomp[[1]] <- .genecomp(datset[[1]], datset[[i+1]],
                                 data = datset,
                                 loop = i,
                                 shortaxlabel = shortaxlabel,
                                 missdata = missdata,
                                 outgroup = outgroup,
                                 verbose = verbose)
    datset[[1]] <- datsetcomp[[1]]
  }
  if (verbose) {
    message(paste0("Matched result of gene ", names(datset[1]),
                   " is again matched with:\n",
                   paste0(names(datset[-1]), "...", collapse = " ")), "\n")
  }
  for (i in 2:numberdatset) {
    if (numberdatset == 2) {
      if (verbose) {
        message(paste0("Gene comparison will exclude sequence set from ",
                       names(datset[i]), " that is not in ",
                       paste0(names(datset[1]), "...", collapse=" ")), "\n")
      } else {
        message(paste0("Gene comparison will exclude sequence set from ",
                       names(datset[i]), " that is not in the matched result of ",
                       paste0(names(datset[1]), "...", collapse=" ")), "\n")
      }
    }
    # Looping over genes
    datset[[i]] <- .genecomp(datset[[i]], datset[[1]],
                             data = datset,
                             loop = i-1,
                             shortaxlabel = shortaxlabel,
                             missdata = missdata,
                             outgroup = outgroup,
                             verbose = verbose)
  }

  if (missdata == FALSE) {
    for (j in seq_along(datset)) {
      for (i in seq_along(datset[[j]][["species"]])) {
        g <- grepl(datset[[j]][["species"]][i], names(datset_temp[[j]]))
        if (any(g)) {
          n <- names(datset_temp[[j]])[g]
          datset[[j]][["species"]][i] <- n
        }
      }
    }
  }

  if (any(unlist(cf))|any(unlist(aff))|any(unlist(infraspp))) {
    # Putting back the names under cf. and aff. ####
    # Adjusting names with infraspecific taxa
    datset <- .namesback(datset,
                         cf = cf,
                         aff = aff,
                         infraspp = infraspp,
                         rename_cf = nr[["rename_cf"]],
                         rename_aff = nr[["rename_aff"]],
                         rename_infraspp = nr[["rename_infraspp"]],
                         shortaxlabel = shortaxlabel,
                         multispp = FALSE)
  }


  # Removing empty, gap-only columns ####
  if (missdata == FALSE) {
    datset <- .delGaps(datset)
  }

  if (verbose) {
    message("Full gene match is finished!", "",
            sep="\n")
  }

  return(datset)
}
