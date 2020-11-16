# Function to adjust names with cf., aff, or infraspecific taxa
# Used inside catmultGenes and dropSeq functions

# Author: Domingos Cardoso

.namesback <- function(datset,
<<<<<<< HEAD
                       cf = NULL,
                       aff = NULL,
                       infraspp = NULL,
                       rename_cf = NULL,
                       rename_aff = NULL,
                       rename_infraspp = NULL,
=======
                       adjust_cf = NULL,
                       adjust_aff = NULL,
                       infra_spp = NULL,
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
                       shortaxlabel = TRUE,
                       multispp = TRUE) {

  if(multispp){
    # Putting back the names under cf. and aff.
<<<<<<< HEAD
    if(any(unlist(cf))){
      names_temp_orig <- unique(unlist(rename_cf))
=======
    if(any(unlist(adjust_cf))){
      names_temp_orig <- unique(unlist(spp_to_rename_cf))
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
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
<<<<<<< HEAD
    if(any(unlist(aff))){
      names_temp_orig <- unique(unlist(rename_aff))
=======
    if(any(unlist(adjust_aff))){
      names_temp_orig <- unique(unlist(spp_to_rename_aff))
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
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
<<<<<<< HEAD
    if(any(unlist(infraspp))){
      names_temp_orig <- unique(unlist(rename_infraspp))
=======
    if(any(unlist(infra_spp))){
      names_temp_orig <- unique(unlist(infraspp_to_rename))
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
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
<<<<<<< HEAD
    if(any(unlist(cf))){
      names_temp_orig <- unique(unlist(rename_cf))
=======
    if(any(unlist(adjust_cf))){
      names_temp_orig <- unique(unlist(spp_to_rename_cf))
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
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
<<<<<<< HEAD
    if(any(unlist(aff))){
      names_temp_orig <- unique(unlist(rename_aff))
=======
    if(any(unlist(adjust_aff))){
      names_temp_orig <- unique(unlist(spp_to_rename_aff))
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
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
<<<<<<< HEAD
    if(any(unlist(infraspp))){
      names_temp_orig <- unique(unlist(rename_infraspp))
=======
    if(any(unlist(infra_spp))){
      names_temp_orig <- unique(unlist(infraspp_to_rename))
>>>>>>> 1b9658920d6f23287175a314e9fb37660c5a3601
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
