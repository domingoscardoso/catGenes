# Function to adjust names with cf., aff, or infraspecific taxa
# Used inside catmultGenes, catfullGenes and dropSeq functions

# Author: Domingos Cardoso

.adjustnames <- function(datset,
                         adjust_cf = NULL,
                         adjust_aff = NULL,
                         infra_spp = NULL) {

  # Adjusting species labels when they have cf or aff

  if(any(unlist(adjust_cf))){
    #cat("Any gene dataset includes species under \"cf.\"", sep = "\n")
    # Finding species labels to rename
    spp_to_rename_cf <- list()
    for(i in seq_along(datset)){
      if(length(names(datset[[i]])[adjust_cf[[i]]]) > 0){
        spp_to_rename_cf[[i]] <- names(datset[[i]])[adjust_cf[[i]]]
        spp_to_rename_cf <- spp_to_rename_cf[!is.na(spp_to_rename_cf)]
      }
    }
    #cat(unique(unlist(spp_to_rename_cf)), "WERE RENAMED TO", gsub("_cf_", "_cf", unique(unlist(spp_to_rename_cf))), "", sep = "\n")
    spp_labels <- lapply(datset, function(x) gsub("_cf_", "_cf", names(x)))
    for(i in seq_along(datset)){
      names(datset[[i]]) <- spp_labels[[i]]
    }
  }

  if(any(unlist(adjust_aff))){
    #cat("Any gene dataset includes species under \"aff.\"", sep = "\n")
    # Finding species labels to rename
    spp_to_rename_aff <- list()
    for(i in seq_along(datset)){
      if(length(names(datset[[i]])[adjust_aff[[i]]]) > 0){
        spp_to_rename_aff[[i]] <- names(datset[[i]])[adjust_aff[[i]]]
        spp_to_rename_aff <- spp_to_rename_aff[!is.na(spp_to_rename_aff)]
      }
    }
    #cat(unique(unlist(spp_to_rename_aff)), "WERE RENAMED TO", gsub("_aff_", "_aff", unique(unlist(spp_to_rename_aff))), "", sep = "\n")
    spp_labels <- lapply(datset, function(x) gsub("_aff_", "_aff", names(x)))
    for(i in seq_along(datset)){
      names(datset[[i]]) <- spp_labels[[i]]
    }
  }

  # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
  if(any(unlist(infra_spp))){
    #cat("Any gene dataset includes infraspecific taxa", sep = "\n")
    # Finding species labels to rename
    infraspp_to_rename <- list()
    for(i in seq_along(datset)){
      if(length(names(datset[[i]])[infra_spp[[i]]]) > 0){
        infraspp_to_rename[[i]] <- names(datset[[i]])[infra_spp[[i]]]
        infraspp_to_rename <- infraspp_to_rename[!is.na(infraspp_to_rename)]
      }
    }
    spp_labels <- list()
    for (i in seq_along(datset)){
      spp_labels[[i]]  <- names(datset[[i]])

      spp_labels[[i]][infra_spp[[i]]] <- gsub("(_[^_]+)_", "\\1",
                                              spp_labels[[i]][infra_spp[[i]]])
      names(datset[[i]]) <- spp_labels[[i]]
    }
  }

  return(datset)
}
