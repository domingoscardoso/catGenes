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
#'              outgroup = NULL,
#'              verbose = TRUE)
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
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' gene matching search will not be printed in the console in full.
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
#'                       missdata = TRUE,
#'                       verbose = TRUE)
#'
#' outgrouptaxa <- c("Vataireopsis_araroba", "Vataireopsis_speciosa")
#' catdf <- catmultGenes(Luetzelburgia,
#'                       maxspp = FALSE,
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

catmultGenes <- function(...,
                         maxspp = TRUE,
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
         "Find help also with:\n",
         "Domingos Cardoso (JBRJ; cardosobot@gmail.com")
  }

  spp_labels_original <- lapply(datset, function(x) names(x))

  cf <- lapply(datset, function(x) grepl("_cf_", names(x)))
  aff <- lapply(datset, function(x) grepl("_aff_", names(x)))
  spp_temp <- lapply(datset, function(x) gsub("_aff_|_cf_", " ", names(x)))
  infraspp <- lapply(spp_temp, function(x) grepl("[[:upper:]][[:lower:]]+_[[:lower:]]+_[[:lower:]]+", x))


  if (any(unlist(infraspp))) {
    infranames <- unique(as.vector(unlist(spp_labels_original))[as.vector(unlist(infraspp))])
    n <- as.vector(unlist(spp_labels_original))[!as.vector(unlist(infraspp))]
    nu <- unique(gsub("(_[^_]+).*", "\\1", n))
    g <- grepl(paste(nu, collapse = "|"), infranames)

    ni <- infranames[g]
    nn <- unique(n[grepl(paste(gsub("(_[^_]+).*", "\\1", infranames[g]), collapse = "|"), n)])

    if (any(g)) {
      stop(paste0("The following accessions are identified at infraspecific level:\n",
                  ni,
                  "\n\nBUT there are accessions of the same species that are NOT as well fully identified with infraspecific taxa...\n",
                  nn,
                  "\n\nYou should do so!\n\n"),
           "Find help also at DBOSLab-UFBA\n",
           "(Domingos Cardoso; cardosobot@gmail.com)")
    }
  }


  if (any(unlist(cf))|any(unlist(aff))|any(unlist(infraspp))) {

    nr <- .namesTorename(datset,
                         cf = cf,
                         aff = aff,
                         infraspp = infraspp)

    # Adjusting species labels when they have cf. or aff. ####
    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    datset <- .adjustnames(datset,
                           cf = cf,
                           aff = aff,
                           infraspp = infraspp)
  }


  # Stoping if the dataset do not include multiple accession ####
  spp_temp <- lapply(datset, function(x) gsub("(_[^_]+)_.*", "\\1", names(x)))
  dup <- lapply(spp_temp, function(x) duplicated(x))

  if (!any(unlist(dup))) {
    stop(paste0("The loaded alignments do not include species duplicated, with multiple accessions.\n",
                "Please use the function catfullGenes.\n\n"),
         "Find help also at DBOSLab-UFBA\n",
         "(Domingos Cardoso; cardosobot@gmail.com)")
  }


  # Shortening the taxon labels (keeping just the scientific names) in species ####
  # not duplicated with multiple accessions so as to maximize the taxon coverage in
  # the final concatenatenated dataset.
  if (maxspp) {
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


  # Now running genecompmult function in a for loop ####
  if (verbose) {
    message(paste0("Matching first the gene ", names(datset[1]), " with:\n",
                   paste0(names(datset[-1]), "...", collapse = " ")), "\n")
  }
  # Creating an empty list to fill in during the loop iteration
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
    # Looping over genes ####
    datsetcomp[[1]] <- .genecompmult(datset[[1]], datset[[i+1]],
                                     data = datset,
                                     loop = i,
                                     shortaxlabel = shortaxlabel,
                                     missdata = missdata,
                                     outgroup = outgroup,
                                     verbose = verbose)
    datset[[1]] <- datsetcomp[[1]]
  }

  if (verbose) {
    message(paste0("Matched result of gene",
                   names(datset[1]), " is again matched with:\n",
                   paste0(names(datset[-1]), "...", collapse = " ")), "\n")
  }

  for (i in 2:numberdatset) {
    if (verbose) {
      if (numberdatset == 2) {
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
    datset[[i]] <- .genecompmult(datset[[i]], datset[[1]],
                                 data = datset,
                                 loop = i-1,
                                 shortaxlabel = shortaxlabel,
                                 missdata = missdata,
                                 outgroup = outgroup,
                                 verbose = verbose)
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
                         multispp = TRUE)
  }

  # This is to insert original names back when using the arguments ####
  # maxspp = T & shortaxlabel = F
  if (maxspp == T & shortaxlabel == F) {
    for (i in seq_along(datset)) {
      grepl_temp <- !datset[[i]][["species"]] %in% spp_labels_original[[i]]
      ntemp <- datset[[i]][["species"]][grepl_temp]
      suppressWarnings({
        for (j in seq_along(ntemp)) {
          n <- spp_labels_original[[i]][grepl(ntemp[j], spp_labels_original[[i]])]
          if (length(n) > 1) {
            n <- n[gsub("(_[^_]+)_.*", "\\1", n) %in% ntemp[j]]
          }
          if (length(n) == 0) {
            ntemp[j] <- ntemp[j]
          } else {
            ntemp[j] <- n
          }
        }
      })
      datset[[i]][["species"]][grepl_temp] <- ntemp
    }
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
