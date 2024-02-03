# Auxiliary function to adjust names with cf., aff, or infraspecific taxa

# Used inside catmultGenes, catfullGenes and dropSeq functions

# Author: Domingos Cardoso

.adjustnames <- function(datset,
                         cf = NULL,
                         aff = NULL,
                         infraspp = NULL) {

  if (!is.data.frame(datset[[1]])) {
    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))) {
      spp_labels <- lapply(datset, function(x) gsub("_cf_", "_cf", names(x)))
      for (i in seq_along(datset)) {
        names(datset[[i]]) <- spp_labels[[i]]
      }
    }

    if (any(unlist(aff))) {
      spp_labels <- lapply(datset, function(x) gsub("_aff_", "_aff", names(x)))
      for (i in seq_along(datset)) {
        names(datset[[i]]) <- spp_labels[[i]]
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      spp_labels <- list()
      for (i in seq_along(datset)) {
        spp_labels[[i]]  <- names(datset[[i]])

        spp_labels[[i]][infraspp[[i]]] <- gsub("(_[^_]+)_", "\\1",
                                               spp_labels[[i]][infraspp[[i]]])
        names(datset[[i]]) <- spp_labels[[i]]
      }
    }

  }

  if (is.data.frame(datset[[1]])) {
    # Adjusting names in dataset already passed through the catmultGenes
    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))){
      spp_labels <- lapply(datset, function(x) gsub("_cf_", "_cf", x[[1]]))
      for (i in seq_along(datset)) {
        datset[[i]][[1]] <- spp_labels[[i]]
      }
    }

    if (any(unlist(aff))) {
      spp_labels <- lapply(datset, function(x) gsub("_aff_", "_aff", x[[1]]))
      for (i in seq_along(datset)) {
        datset[[i]][[1]] <- spp_labels[[i]]
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      spp_labels <- list()
      for (i in seq_along(datset)) {
        spp_labels[[i]]  <- datset[[i]][[1]]

        spp_labels[[i]][infraspp[[i]]] <- sub("(_[^_]+)_", "\\1",
                                              spp_labels[[i]][infraspp[[i]]])
        datset[[i]][[1]] <- spp_labels[[i]]
      }
    }

  }

  return(datset)
}
