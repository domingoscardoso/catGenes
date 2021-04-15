# Auxiliary function for returning TRUE if all values are equal and FALSE if
# it contains different values

# Used inside the function dropSeq

# Author: Domingos Cardoso

equalnumb <- function(x) {
  res <- FALSE
  x <- na.omit(as.vector(x))
  if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
  return(res)

}
