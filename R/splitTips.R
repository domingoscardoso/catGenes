#' Parse tip labels into a dataframe with multiple columns for each content
#'
#' @author Domingos Cardoso
#'
#' @description This function is designed to split tip labels into a dataframe
#' object. Each column in the dataframe corresponds to the genus, species,
#' infraspecific name, and associated voucher and GenBank information, whenever
#' available. If doubtful particles such as 'aff.' or 'cf.' are provided, they
#' will also be separated into a distinct column.
#'
#' @usage
#' splitTips(tiplabels = NULL)
#'
#' @param tiplabels A vector of tip labels from a phylogeny.
#'
#' @return A dataframe.
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(Harpalyce_bayes_tree)
#'
#' df <- splitTips(tiplabels = Harpalyce_bayes_tree@phylo$tip.label)
#'}
#'
#' @importFrom tibble add_column
#' @importFrom magrittr %>%
#'
#' @export
#'

splitTips <- function(tiplabels = NULL) {

  if (any(grepl("_", tiplabels))) {
    tiplabels <- gsub("_", " ", tiplabels)
  }

  # Split the vector using strsplit
  split_vector <- strsplit(tiplabels, " ")
  # Find the maximum length of the split vectors
  max_length <- max(sapply(split_vector, length))
  # Pad the shorter vectors with NA
  padded_vectors <- lapply(split_vector,
                           function(x) c(x, rep(NA, max_length - length(x))))
  # Convert the list to a data frame
  df <- as.data.frame(do.call(rbind, padded_vectors))

  names(df)[1] <- "genus"

  tf <- grepl("^aff$|^aff[.]$|^cf$|^cf[.]$", df[[2]])
  if (any(tf)) {
    # Add a new column after the first column
    df <- cbind(df[, 1, drop = FALSE], doubtID = tf, df[, -1])
    # Find rows where the condition is TRUE
    rows_to_shift <- which(df$doubtID)
    # Loop through each specified row and shift values to the left
    for (i in rows_to_shift) {
      true_index <- which(df[i, ] == TRUE)[1]
      df[i, (true_index):(ncol(df) - 1)] <- df[i, (true_index + 1):ncol(df)]
      df[i, ncol(df)] <- NA
    }
    df$doubtID[which(df$doubtID == FALSE)] <- NA
    names(df)[3] <- "species"
  } else {
    names(df)[2] <- "species"
  }

  # Define the character to search for GenBank accession number
  search_character <- "^([[:upper:]]){2}([[:digit:]]){6}$"

  # Check for the specific character in each row and move it to a newly created
  # last column.
  df$genbank <- NA  # Create a new column to store the matched values

  for (i in 1:nrow(df)) {
    if (any(grepl(search_character, df[i, ]))) {
      match_index <- which(grepl(search_character, df[i, ]))[1]
      df[i, "genbank"] <- df[i, match_index]
      df[i, match_index] <- NA
    }
  }

  # Remove all empty columns
  rm <- colSums(is.na(df))<nrow(df)
  n <- names(rm)[which(rm == TRUE)][-1]
  df <- df[, rm]

  tf_genbank <- "genbank" %in% names(df)
  if (tf_genbank) {
    genbank <- df$genbank
    df <- df[, !names(df) %in% "genbank"]
  }

  ll <- which(names(df) %in% "species") + 1

  tf <- grepl("^[[:lower:]]+$", df[[ll]])
  if (any(tf)) {
    # Add a new column after the first column
    df <- df %>% tibble::add_column(infrasp = tf, .before = ll)
    # Find rows where the condition is TRUE
    rows_to_shift <- which(df$infrasp)
    # Loop through each specified row and shift values to the left
    for (i in rows_to_shift) {
      true_index <- which(df[i, ] == TRUE)[1]
      df[i, (true_index):(ncol(df) - 1)] <- df[i, (true_index + 1):ncol(df)]
      df[i, ncol(df)] <- NA
    }
    df$infrasp[which(df$infrasp == FALSE)] <- NA
  }

  tf <- grepl("^V[[:digit:]]", names(df))

  names(df)[which(tf)][1] <- "voucher"

  # Concatenate the selected columns and create a new column 'concatenated_result'
  df$voucher <- apply(df[which(tf)], 1, function(x) paste0(x, collapse = " "))

  df$voucher <- gsub("(\\sNA){1,}$.*|(^NA\\s.*)", "", df$voucher)
  df$voucher <- gsub("\\s", "\\\u00ad", df$voucher) # add a proper hyphen
  df$voucher <- gsub("^$", NA, df$voucher)

  tf <- grepl("^V[[:digit:]]", names(df))
  df <- df[, which(tf == FALSE)]

  tf <- !"genbank" %in% names(df)
  if (tf & tf_genbank) {
    df <- cbind(df, genbank)
  }

  df <- df %>% tibble::add_column(tiplabels = tiplabels, .before = "genus")

  return(df)
}


# Regular expressions
# https://evoldyn.gitlab.io/evomics-2018/ref-sheets/R_strings.pdf

# # remove all after second space
# d$species[tf] <- sub("^(\\S*\\s+\\S+).*", "\\1", d$tiplabels[tf])
# # remove all before second space
# d$voucher[tf] <- sub("\\w+\\s\\w+\\s", "", d$tiplabels[tf])
#
# d$genus <- gsub("\\s.*", "", d$species) # remove all after first space
# d$species <- sub(".*? ", "", d$species) # remove all before first space
# d$voucher <- sub("\\s+[^ ]+$", "", d$voucher) # remove all after last space
# d$voucher <- gsub("\\s", "\\\u00ad", d$voucher) # add a proper hyphen
