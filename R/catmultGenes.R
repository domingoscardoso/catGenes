#' Compares gene datasets for combined phylogenetic analysis when species are
#' duplicated or represented by multiple accessions in one DNA alignment
#'
#' @author Domingos Cardoso and Quezia Cavalcante
#'
#' @description Compares a list of "n" gene datasets (individual DNA alignments)
#' and makes them with the same number of taxa, ready for combined, multigene
#' phylogenetic analysis. This function is best designed for concatenating DNA
#' alignments where species have duplicated sequences (multiple accessions) from
#' different collections. Then, make sure the species are labeled with both the
#' scientific name and the same identifying number throughout each DNA alignment.
#' Identifying numbers could be collector surname and associated collection number,
#' or an accession number for the isolated DNA from which each gene was sequenced.
#' During the comparison across DNA alignments for concatenation, the function will
#' consider that any species is represented by multiple sequences and so in each
#' individual gene dataset species will fully matched if they have exact scientific
#' name and associated identifying number.
#'
#' @usage
#' catmultGenes(\dots,
#'              maxspp = TRUE,
#'              shortaxlabel = TRUE,
#'              missdata = TRUE,
#'              outgroup = NULL)
#'
#' @param ... a list of NEXUS-formatted gene datasets as read by ape's \code{\link{read.nexus.data}}
#' or at least two individually ape-read objects of NEXUS-formatted gene datasets.
#'
#' @param maxspp Logical, if \code{FALSE} any species never duplicated with
#' multiple accessions in any indvidual DNA alignment might end either duplicated
#' or deleted, depending on the chosen missdata argument. We recomend to maxspp =
#' \code{TRUE} so as to maximize the taxon coverage. This means that if the species
#' is not duplicated in any individual dataset, it will always be kept in the final
#' concatenated dataset no matter each sequence for that species across the individual
#' dataset were generated from distinct collections or accessions.
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
#' @seealso \code{\link{dropSeq}}
#' @seealso \code{\link{nexusdframe}}
#' @seealso \code{\link{phylipdframe}}
#' @seealso \code{\link{fastadframe}}
#'
#' @examples \dontrun{
#' data(Luetzelburgia)
#' catdf <- catmultGenes(Luetzelburgia,
#'                       maxspp = TRUE,
#'                       shortaxlabel = TRUE,
#'                       missdata = TRUE)
#'
#' outgrouptaxa <- c("Vataireopsis_araroba", "Vataireopsis_speciosa")
#' catdf <- catmultGenes(Luetzelburgia,
#'                       maxspp = FALSE,
#'                       shortaxlabel = TRUE,
#'                       missdata = FALSE,
#'                       outgroup = outgrouptaxa)
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom stats setNames na.omit
#'
#' @export

catmultGenes <- function(...,
                         maxspp = TRUE,
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

  cf <- lapply(datset, function(x) grepl("_cf_", names(x)))
  aff <- lapply(datset, function(x) grepl("_aff_", names(x)))
  spp_temp <- lapply(datset, function(x) gsub("_aff_|_cf_", " ", names(x)))
  infraspp <- lapply(spp_temp, function(x) grepl("[[:upper:]][[:lower:]]+_[[:lower:]]+_[[:lower:]]+",
                                                x))

  if(any(unlist(cf))|any(unlist(aff))|any(unlist(infraspp))){

    nr <- .namesTorename(datset,
                         cf = cf,
                         aff = aff,
                         infraspp = infraspp)

    # Adjusting species labels when they have cf. or aff.
    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    datset <- .adjustnames(datset,
                           cf = cf,
                           aff = aff,
                           infraspp = infraspp)
  }

  # Shortening the taxon labels (keeping just the scientific names) in species
  # not duplicated with multiple accessions so as to maximize the taxon coverage in
  # the final concatenatenated dataset.
  if(maxspp){
    datset_temp <- datset
    spp_labels <- lapply(datset_temp, function(x) gsub("(_[^_]+)_.*", "\\1", names(x)))
    for (i in seq_along(datset_temp)) {
      names(datset_temp[[i]]) <- spp_labels[[i]]
    }
    dup_temp <- list()
    dup_spp <- list()
    nondup_spp_temp <- list()
    for (i in seq_along(datset_temp)) {

      dup_temp[[i]] <- c(duplicated(names(datset_temp[[i]]), fromLast = TRUE) |
                          duplicated(names(datset_temp[[i]])))
      dup_spp[[i]] <- unique(names(datset_temp[[i]])[dup_temp[[i]]])
      nondup_spp_temp[[i]] <- names(datset_temp[[i]])[!dup_temp[[i]]]
    }

    dup_spp <- unique(unlist(dup_spp))
    nondup_spp <- list()
    spp_labels <- list()
    for (i in seq_along(datset_temp)) {

      nondup_spp[[i]] <- nondup_spp_temp[[i]][!nondup_spp_temp[[i]] %in% dup_spp]

      spp_labels[[i]]  <- names(datset[[i]])

      spp_labels[[i]][names(datset_temp[[i]]) %in% nondup_spp[[i]]] <- nondup_spp[[i]]

      names(datset[[i]]) <- spp_labels[[i]]
    }
  }


  # Now running genecompmult function in a for loop

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
    datsetcomp[[1]] <- .genecompmult(datset[[1]], datset[[i+1]],
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
    datset[[i]] <- .genecompmult(datset[[i]], datset[[1]],
                                 data = datset,
                                 loop = i-1,
                                 shortaxlabel = shortaxlabel,
                                 missdata = missdata,
                                 outgroup = outgroup)

  }

  if(any(unlist(cf))|any(unlist(aff))|any(unlist(infraspp))){
    # Putting back the names under cf. and aff.
    # Adjusting names with infraspecific taxa
    datset <- .namesback(datset,
                         cf = cf,
                         aff = aff,
                         infraspp = infraspp,
                         rename_cf = nr[["rename_cf"]],
                         rename_aff = nr[["rename_aff"]],
                         rename_infraspp = nr[["rename_infraspp"]],
                         shortaxlabel = shortaxlabel,
                         multispp = TRUE)
  }

  cat("Full gene match is finished!", "",
      sep="\n")

  return(datset)
}
