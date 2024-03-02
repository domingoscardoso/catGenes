# Auxiliary function to get original names with cf., aff., or infraspecific taxa

# Used inside catmultGenes, catfullGenes and dropSeq functions

# Author: Domingos Cardoso

.namesTorename <- function(datset,
                           cf = NULL,
                           aff = NULL,
                           infraspp = NULL) {

  if (!is.data.frame(datset[[1]])) {
    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))) {
      # Finding species labels to rename
      rename_cf <- list()
      for (i in seq_along(datset)) {
        if (length(names(datset[[i]])[cf[[i]]]) > 0) {
          rename_cf[[i]] <- names(datset[[i]])[cf[[i]]]
          rename_cf <- rename_cf[!is.na(rename_cf)]
        }
      }
    }

    if (any(unlist(aff))) {
      # Finding species labels to rename
      rename_aff <- list()
      for (i in seq_along(datset)) {
        if (length(names(datset[[i]])[aff[[i]]]) > 0) {
          rename_aff[[i]] <- names(datset[[i]])[aff[[i]]]
          rename_aff <- rename_aff[!is.na(rename_aff)]
        }
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      # Finding species labels to rename
      rename_infraspp <- list()
      for (i in seq_along(datset)) {
        if (length(names(datset[[i]])[infraspp[[i]]]) > 0) {
          rename_infraspp[[i]] <- names(datset[[i]])[infraspp[[i]]]
          rename_infraspp <- rename_infraspp[!is.na(rename_infraspp)]
        }
      }
    }

  }

  if (is.data.frame(datset[[1]])) {
    # Adjusting names in dataset already passed through the catmultGenes

    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))){
      # Finding species labels to rename
      rename_cf <- list()
      for (i in seq_along(datset)) {
        if (length(datset[[i]][[1]][cf[[i]]]) > 0) {
          rename_cf[[i]] <- datset[[i]][[1]][cf[[i]]]
          rename_cf <- rename_cf[!is.na(rename_cf)]
        }
      }
    }

    if (any(unlist(aff))) {
      # Finding species labels to rename
      rename_aff <- list()
      for (i in seq_along(datset)) {
        if (length(datset[[i]][[1]][aff[[i]]]) > 0) {
          rename_aff[[i]] <- datset[[i]][[1]][aff[[i]]]
          rename_aff <- rename_aff[!is.na(rename_aff)]
        }
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      # Finding species labels to rename
      rename_infraspp <- list()
      for (i in seq_along(datset)) {
        if (length(datset[[i]][[1]][infraspp[[i]]]) > 0) {
          rename_infraspp[[i]] <- datset[[i]][[1]][infraspp[[i]]]
          rename_infraspp <- rename_infraspp[!is.na(rename_infraspp)]
        }
      }
    }

  }

  orignames <- list(
    if (any(unlist(cf))) {
      rename_cf
    },
    if (any(unlist(aff))) {
      rename_aff
    },
    if (any(unlist(infraspp))) {
      rename_infraspp
    })

  names(orignames) <- c("rename_cf", "rename_aff", "rename_infraspp")

  return(orignames)
}
