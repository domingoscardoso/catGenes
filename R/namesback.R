# Function to adjust names with cf., aff, or infraspecific taxa
# Used inside catmultGenes and dropSeq functions

# Author: Domingos Cardoso

.namesback <- function(datset,
                       cf = NULL,
                       aff = NULL,
                       infraspp = NULL,
                       rename_cf = NULL,
                       rename_aff = NULL,
                       rename_infraspp = NULL,
                       shortaxlabel = TRUE,
                       multispp = TRUE) {

  if(multispp){
    # Putting back the names under cf. and aff.
    if(any(unlist(cf))){
      names_temp_orig <- unique(unlist(rename_cf))
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
    if(any(unlist(aff))){
      names_temp_orig <- unique(unlist(rename_aff))
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
    if(any(unlist(infraspp))){
      names_temp_orig <- unique(unlist(rename_infraspp))
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
  }

  if(!multispp){
    # Putting back the names under cf. and aff.
    if(any(unlist(cf))){
      names_temp_orig <- unique(unlist(rename_cf))
      if(shortaxlabel){
        names_temp_orig <- gsub("(_[^_]+)_.*", "\\1", names_temp_orig)
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
    if(any(unlist(aff))){
      names_temp_orig <- unique(unlist(rename_aff))
      if(shortaxlabel){
        names_temp_orig <- gsub("(_[^_]+)_.*", "\\1", names_temp_orig)
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
    if(any(unlist(infraspp))){
      names_temp_orig <- unique(unlist(rename_infraspp))
      if(shortaxlabel){
        names_temp_orig <- gsub("(_[^_]+)_.*", "\\1", names_temp_orig)
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

  }

  return(datset)
}
