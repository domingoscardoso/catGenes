# Function to adjust names with cf., aff, or infraspecific taxa
# Used inside catmultGenes and dropSeq functions

# Author: Domingos Cardoso

.namesback <- function(datset,
                       shortaxlabel = TRUE) {

  # Putting back the names under cf. and aff.
  if(any(unlist(adjust_cf))){
    names_temp_orig <- unique(unlist(spp_to_rename_cf))
    if(shortaxlabel){
      names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
    }
    names_temp <- gsub("_cf_", "_cf", names_temp_orig)
    n = 0
    for (i in names_temp){
      n = n + 1
      spp_labels <- lapply(datset, function(x) gsub(i, names_temp_orig[n], x[[1]]))

      for(j in seq_along(datset)){
        datset[[j]][[1]] <- spp_labels[[j]]
      }
    }
  }
  if(any(unlist(adjust_aff))){
    names_temp_orig <- unique(unlist(spp_to_rename_aff))
    if(shortaxlabel){
      names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
    }
    names_temp <- gsub("_aff_", "_aff", names_temp_orig)
    n = 0
    for (i in names_temp){
      n = n + 1
      spp_labels <- lapply(datset, function(x) gsub(i, names_temp_orig[n], x[[1]]))

      for(j in seq_along(datset)){
        datset[[j]][[1]] <- spp_labels[[j]]
      }
    }
  }

  # Adjusting names with infraspecific taxa
  if(any(unlist(infra_spp))){
    names_temp_orig <- unique(unlist(infraspp_to_rename))
    if(shortaxlabel){
      names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
    }
    names_temp <- gsub("(_[^_]+)_", "\\1", names_temp_orig)
    n = 0
    for (i in names_temp){
      n = n + 1
      spp_labels <- lapply(datset, function(x) gsub(i, names_temp_orig[n], x[[1]]))

      for(j in seq_along(datset)){
        datset[[j]][[1]] <- spp_labels[[j]]
      }
    }
  }

  return(datset)
}
