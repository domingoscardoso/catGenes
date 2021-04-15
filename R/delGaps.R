# Auxiliar function to delete gap-only columns

# Used inside the functions catfullGenes and catmultGenes

# Author: Domingos Cardoso

.delGaps <- function(x) {

  dgaps <- list()
  dgaps_temp <- list()

  for (j in seq_along(x)) {
    temp <- list()
    for (i in seq_along(x[[j]][["sequence"]])) {
      temp[[i]] <-  which(grepl("-", strsplit(x[[j]][["sequence"]], "")[[i]]) == T)
    }
    dgaps_temp[[j]] <- temp

    # List of all columns to be removed
    dgaps[[j]] <- .intersectAll(dgaps_temp[[j]])
  }

  # Removing all gap-only columns
  for (j in seq_along(x)) {
    if (length(dgaps[[j]]) != 0) {

      for (i in seq_along(x[[j]][["sequence"]])) {

        s <- strsplit(x[[j]][["sequence"]], "")[[i]][-dgaps[[j]]]

        x[[j]][["sequence"]][i] <- paste(s, collapse = "")
      }
    }
  }

return(x)
}
