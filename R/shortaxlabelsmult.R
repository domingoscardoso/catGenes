# Auxiliary function to delete all extra identifiers of the DNA sequence (e.g. GenBank numbers)
# so as to maintain just the scientific names and associated collector number

# Used inside the function catmultGenes

# Author: Domingos Cardoso

.shortaxlabelsmult <- function(x) {

  colnames(x) <- c("species", "sequence")
  x$species <- gsub("(_[^_]+_[^_]+)_.*", "\\1", x$species)

  return(x)
}
