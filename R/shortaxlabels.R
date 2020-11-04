# Function to delete all identifiers of the DNA sequence (e.g. collector and GenBank numbers)
# so as to maintain just the scientific names

# Used inside genecomp function

# Author: Domingos Cardoso

.shortaxlabels <- function(x) {

  colnames(x) <- c("species", "sequence")
  x$species <- gsub("(_[^_]+)_.*", "\\1", x$species)

  return(x)
}
