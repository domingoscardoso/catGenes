# Small functions to evaluate user input for the main functions and
# return meaningful errors.
# Author: Domingos Cardoso

#_______________________________________________________________________________
# Check if the dir input is "character" type and if it has a "/" in the end
.arg_check_dir <- function(x) {
  # Check classes
  class_x <- class(x)
  if (!"character" %in% class_x) {
    stop(paste0("The argument dir should be a character, not '", class_x, "'."),
         call. = FALSE)
  }
  if (grepl("[/]$", x)) {
    x <- gsub("[/]$", "", x)
  }
  return(x)
}


#_______________________________________________________________________________
# Check path
.arg_check_inputdf <- function(inputdf, gb.colnames) {

  temp <- inputdf[gb.colnames]
  tf <- lapply(temp, function(x) duplicated(na.omit(x)))
  temp <- lapply(seq_along(temp), function(i) na.omit(temp[[i]])[tf[[i]]])
  n <- unlist(lapply(seq_along(temp), function(i) length(temp[[i]])))
  gb.colnames_dup <- gb.colnames[n > 0]
  list_dup <- vector()
  for (i in seq_along(gb.colnames_dup)) {
    list_dup[i] <- paste0(gb.colnames_dup[i], ": ",
                          paste0(temp[n > 0][[i]], collapse = ", "))
  }

  if (any(n > 0)) {
    stop("The following accessions numbers are duplicated in your dataset:\n\n",
         paste0(list_dup, collapse = "\n\n"),
         call. = FALSE)
  }

}

