#' Abbreviate genus name
#'
#' @author Domingos Cardoso
#'
#' @description This function abbreviates genus names. In cases where multiple
#' genera share similarities in their first, second, and third letters, the
#' function will generate an abbreviation, extending up to the fourth letter for
#' disambiguation.
#'
#' @usage
#' abbrevGen(tiplabels = NULL)
#'
#' @param tiplabels A vector of genus names or tip labels to be abbreviated.
#'
#' @param abbrevfull A logical value indicating whether to fully abbreviate genus
#' names, if set to \code{TRUE}. If all genera have different initial letters,
#' abbreviations will be until the first letter. If at least two different genera
#' are similar in the first letter, then all genera will be abbreviated until the
#' second letter.  If at least two different genera  are similar in the first two
#' letters, then all genera will be abbreviated until the third letter.
#'
#' @param abbrevmult A logical value indicating whether to fully abbreviate genus
#' names in different ways depending on how they differ to each other based on
#' the initial letters. If set to \code{TRUE}, in cases where multiple genera
#' share similarities in their first, second, and third letters, the function
#' will generate an abbreviation, extending up to the fourth letter for
#' disambiguation.
#'
#' @return A dataframe.
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(Harpalyce_bayes_tree)
#'
#' df <- abbrevGen(tiplabels = Harpalyce_bayes_tree@phylo$tip.label)
#'}
#'
#' @export
#'

abbrevGen <- function(tiplabels = NULL,
                      abbrevfull = TRUE,
                      abbrevmult = FALSE) {
  tiplabels <- gsub("_", " ", tiplabels)
  genus_vector <- gsub("\\s.*|_.*", "", tiplabels)

  abbrev_df <- data.frame(original_tiplabels = tiplabels,
                          original_genus = genus_vector,
                          abbreviation = NA,
                          stringsAsFactors = FALSE)

  genera <- unique(genus_vector)

  if (abbrevfull) {
    # Check if all genera have different first letters
    all_diff_first_letters <- length(unique(substr(genera, 1, 1))) == length(genera)

    # Check if there are at least two genera with the same first two letters
    same_two_letters <- any(duplicated(substr(genera, 1, 2)))

    # Determine the abbreviation length based on the conditions
    abbreviation_length <- ifelse(all_diff_first_letters, 1,
                                  ifelse(same_two_letters, 3, 2))

    # Loop through each genus in the vector
    for (i in seq_along(genera)) {
      # Abbreviate the genus name based on the determined length
      abbreviation <- paste0(substr(genera[i], 1, abbreviation_length), ".")
      abbrev_df$abbreviation[genus_vector %in% genera[i]] <- abbreviation
    }

  } else if (abbrevmult) {

    for (i in seq_along(genera)) {

      # Find other genera with the same starting letter
      same_letter_genera <- grep(paste0("^", substr(genera[i], 1, 1)),
                                 genera, value = TRUE)

      # Determine the abbreviation based on the number of genera with the same
      # starting letter.
      if (length(same_letter_genera) == 1) {
        abbreviation <- paste0(substr(genera[i], 1, 1), ".")
      } else if (length(same_letter_genera) > 1) {
        abbreviation <- paste0(substr(genera[i], 1, 2), ".")
      }
      abbrev_df$abbreviation[genus_vector %in% genera[i]] <- abbreviation
    }

    unique_rows <- abbrev_df[!duplicated(abbrev_df$original_genus), ]

    temp <- .abbrev.further(unique_rows, nletters = 3)
    temp <- .abbrev.further(temp, nletters = 4)

    for (i in seq_along(temp$original_genus)) {
      tf <- abbrev_df$original_genus %in% temp$original_genus[i]
      tftf <- !abbrev_df$abbreviation[tf] %in% temp$abbreviation[i]
      if (any(tftf)) {
        abbrev_df$abbreviation[tf] <- temp$abbreviation[i]
      }
    }

  }

  abbrev_df$abbrev_tiplabels <- paste(abbrev_df$abbreviation,
                                      sub(".*? ", "", abbrev_df$original_tiplabels))
  return(abbrev_df)
}


.abbrev.further <- function(genera, nletters) {

  tf <- duplicated(genera$abbreviation)
  if (any(tf)) {
    temp <- genera$abbreviation[tf]
    tf <- genera$abbreviation %in% temp
    genera$abbreviation[tf] <- paste0(substr(genera$original_genus[tf], 1,
                                             nletters), ".")
  }
  return(genera)
}

