# Auxiliary function to obtain the intersection of a list of vectors

# Used inside the function delGaps

# Author: Domingos Cardoso

.intersectAll <- function(...) {
  args <- list(...)
  nargs <- length(args[[1]])
  if (nargs <= 1) {
    if (nargs == 1 && is.list(unlist(args[[1]][[1]]))) {
      do.call("intersectAll", unlist(args[[1]][[1]]))
    } else {
      stop("Cannot evaluate intersection fewer than 2 arguments")
    }
  } else if (nargs == 2) {
    intersect(unlist(args[[1]][[1]]), unlist(args[[1]][[2]]))
  } else {
    intersect(unlist(args[[1]][[1]]), .intersectAll(args[[1]][-1]))
  }
}
