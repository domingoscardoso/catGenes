# Function to create a named list of the input individual genes.

# Used inside the functions catfullGenes, catmultGenes, writeNexus, and writePhylip.

# Author: Domingos Cardoso

.namedlist <- function(...) {
  nms <- sapply(as.list(substitute(list(...))), deparse)[-1]
  setNames(list(...), nms)
}
