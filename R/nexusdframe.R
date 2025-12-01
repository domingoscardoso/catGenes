#' Writes a NEXUS-formatted DNA alignment from a dataframe-formatted DNA alignment
#'
#' @author Domingos Cardoso
#'
#' @description Writes \code{data.frame} formatted DNA alignment or \code{list}
#' formatted NEXUS file as originally imported with \code{\link{ape}}'s function
#' \code{\link{read.nexus.data}} into a NEXUS-formatted file. It is useful for
#' writing each gene dataset from within the resulting list of compared gene datasets,
#' after running the concatenating functions \code{\link{catfullGenes}} and
#' \code{\link{catmultGenes}}. The function is also useful for saving the original
#' list-formatted NEXUS object as read by \code{\link{read.nexus.data}}, after
#' making specific changes in such original individual alignment (e.g. corrections
#' of species names).
#'
#' @usage
#' nexusdframe(x, file,
#'             dropmisseq = TRUE,
#'             endgaps.to.miss = TRUE)
#'
#' @param x The object to be written, any two-column-sized \code{data.frame} where
#' the first column contains the taxon names and the second column the DNA sequence.
#' Otherwise, the object may be a list-formatted NEXUS file as originally
#' imported with \code{\link{ape}}'s function \code{\link{read.nexus.data}}.
#'
#' @param file Either a character string naming a file or a \code{\link{connection}}
#'  open for writing.
#'
#' @param dropmisseq Logical, if \code{FALSE} the function will not drop species
#' with empty DNA sequence. After running the concatenating function
#' \code{\link{catmultGenes}} using missdata = \code{TRUE}, and then using
#' \code{\link{dropSeq}} to remove duplicated accessions of the same species,
#' you might find useful to keep dropmisseq = \code{TRUE} so as to save each
#' individual DNA alignment by also removing species that fully miss the sequence
#' data.
#'
#' @param endgaps.to.miss Logical, if \code{FALSE} the function will not replace
#' terminal GAPs into missing character (?).
#'
#' @seealso \code{\link{catfullGenes}}
#' @seealso \code{\link{catmultGenes}}
#' @seealso \code{\link{dropSeq}}
#'
#' @examples \dontrun{
#' data(Gaya)
#' catdf <- catfullGenes(Gaya,
#'                       multiaccessions = FALSE,
#'                       shortaxlabel = TRUE,
#'                       missdata = FALSE,
#'                       outgroup = "Abutilon_costicalyx")
#'
#' ITS <- catdf[[1]]
#' petLpsbE <- catdf[[2]]
#' rpl16 <- catdf[[3]]
#'
#' nexusdframe(ITS, file = "filename.nex")
#' nexusdframe(petLpsbE, file = "filename.nex")
#' nexusdframe(rpl16, file = "filename.nex")
#' }
#'
#' @importFrom stringr str_trunc str_pad str_extract_all
#' @importFrom tidyr unite
#' @importFrom stats setNames
#' @importFrom magrittr "%>%"
#'
#' @export

nexusdframe <- function(x, file,
                        dropmisseq = TRUE,
                        endgaps.to.miss = TRUE) {

  if (inherits(x, "list")) {
    for (i in 1:length(x)) {
      x[[i]] <- paste(x[[i]], collapse = "")
      x[[i]] <- toupper(x[[i]])
    }
    spp <- stats::setNames(data.frame(names(x)), "species")
    seqs <- stats::setNames(data.frame(unlist(x, use.names = FALSE)), "sequence")
    x <- cbind(spp, seqs)
  }

  if (!is.data.frame(x)) {
    x <- data.frame(x)
  }

  ncolumns <- length(x)
  if (ncolumns == 1) {
    stop("You must provide a two-column-sized data.frame containing the taxon names and DNA sequences.\n\n",
         "Find help also with:\n",
         "Domingos Cardoso (JBRJ; cardosobot@gmail.com")
  }

  if (names(x)[1] != "species" | names(x)[2] != "sequence") {
    names(x) <- c("species", "sequence")
  }

  if (dropmisseq) {
    # Get number of missing data "N" and "?" in each sequence
    missdata <- vector()
    missdataN <- vector()
    misstotal_temp <- list()
    numbchar <- nchar(x[1, 2])
    for (i in seq_along(x$sequence)) {
      missdata <- length(stringr::str_extract_all(x$sequence[i], "[?]", simplify = FALSE)[[1]])
      missdataN <- length(stringr::str_extract_all(x$sequence[i], "N", simplify = FALSE)[[1]])
      misstotal_temp[i] <- missdata + missdataN
    }
    # Dropping species with empty seqs
    x <- x[!unlist(misstotal_temp) %in% numbchar, ]
  }

  # Replace terminal GAPs into missing character (?)
  if (endgaps.to.miss) {
    x <- .replace_terminal_gaps(x)
  }

  writtenby <- paste("[catGenes (DBOSLab-UFBA), ", date(), "]\n\n", sep = "")
  nexus <- paste("#NEXUS", "", sep="\n")
  begindata <- paste("BEGIN DATA;")
  dimname <- paste("DIMENSIONS")
  ntax <- paste0("NTAX=", length(rownames(x)))
  numbchar <- paste0("NCHAR=", nchar(as.character(x[1, 2])), ";")
  dimensions <- paste(dimname, ntax, numbchar, sep = " ")
  format <- paste("FORMAT DATATYPE=DNA GAP=- MISSING=?;")
  matrix <- paste("MATRIX", "", sep = "\n")
  n <- nchar(x[1, 2])

  # Calculating the space between the taxon labels and corresponding DNA sequence
  vector.list1 <- vector("list")
  for (i in x$species) {
    vector.list1[[i]] <- nchar(i)
  }
  numtaxlab <- as.data.frame(unlist(vector.list1, use.names = FALSE))
  colnames(numtaxlab) <- "numtaxlab"
  row.names(numtaxlab) <- x$species

  x$species <- x$species %>%
    stringr::str_trunc(max(numtaxlab$numtaxlab, na.rm = FALSE) + 5) %>% # Just increase this last number if we want to add more space
    stringr::str_pad(max(numtaxlab$numtaxlab, na.rm = FALSE) + 5, "right")
  x <- tidyr::unite(x, "unamed", colnames(x), sep = "")
  end <- paste(";","END;", sep = "\n")
  x[nrow(x) + 1, ] <- end

  # Creating the scale of char length
  x <- rbind(x[rep(1, 1), ], x)
  npad <- nchar(sub("\\s.*", "", x[1, ]))
  x[1, ] <- sub(".*?\\s", "", x[1, ])
  x[1, ] <- paste0(paste(rep(" ", npad), collapse = ""), x[1, ])
  x <- rbind(x[rep(1, 1), ], x)
  x[1, ] <- gsub("^","[", x[1, ])
  x[1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]","", x[1, ])

  x[2, ] <- gsub("^","[", x[2, ])
  x[2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]","", x[2, ])

  nbr <- paste(seq(0, n+20, by = 10))
  nbrseq <- vector()
  dotseq <- vector()
  for (i in seq_along(nbr)) {
    if (nbr[i] == 0) {
      nbrseq[i] <- paste(1, paste(rep(" ", 8-nchar(nbr[i])), collapse = ""), collapse = "")
      dotseq[i] <- paste(".", paste(rep(" ", 8-nchar(nbr[1])), collapse = ""), collapse = "")
    } else {
      nbrseq[i] <- paste(nbr[i], paste(rep(" ", 9-nchar(nbr[i])), collapse = ""), collapse = "")
      dotseq[i] <- paste(".", paste(rep(" ", 9-nchar(nbr[1])), collapse = ""), collapse = "")
    }
  }

  x[1, ] <- paste0(x[1, ], paste0(paste(nbrseq, collapse = ""), "]"))
  x[2, ] <- paste0(x[2, ], paste0(paste(dotseq, collapse = ""), "]"))

  colnames(x) <- paste(nexus,
                       writtenby,
                       begindata,
                       dimensions,
                       format,
                       matrix,
                       sep = "\n")

  zz <- file(file, "w")
  write.table(x, zz,
              append = FALSE, quote = FALSE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE)
  close(zz)
}
