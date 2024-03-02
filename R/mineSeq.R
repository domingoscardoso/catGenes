#' Read and download DNA sequences from GenBank
#'
#' @author Domingos Cardoso
#'
#' @description An \code{\link{ape}}-based function to connect with the
#' [GenBank](https://www.ncbi.nlm.nih.gov/genbank) database, read nucleotide
#' sequences using accession numbers given as arguments, and write them in a
#' fasta format file.
#'
#' @usage
#' mineSeq(inputdf = NULL,
#'         gb.colnames = NULL,
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         dir = "RESULTS_mineSeq",
#'         filename = "GenBanK_seqs")
#'
#' @param inputdf A dataframe object containing the taxon names in a 'Species'
#' column, the voucher information in 'Voucher' column, and the GenBank accessions
#' for each genes in separate columns named by the corresponding gene.
#'
#' @param gb.colnames A vector with column names within the \code{inputdf}
#' dataframe corresponding to each gene, where the GenBank accession numbers are
#' listed.
#'
#' A vector with column names in the "inputdf" dataframe corresponding to each gene, identified by GenBank accession numbers.
#'
#' @param as.character a logical controlling whether to return the sequences as
#' an object of class "DNAbin" (the default).

#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the GenBank search will not be printed in the console in full.
#'
#' @param save Logical, if \code{TRUE}, the edited tree will be saved on disk.
#'
#' @param dir Pathway to the computer's directory, where the mined DNA sequences
#' in a fasta format file will be saved provided that the argument \code{save}
#' is set up in \code{TRUE}. The default is to create a directory named
#' **RESULTS_mineSeq** and the sequences will be saved within a subfolder named
#' after the current date.
#'
#' @param filename Name of the output file to be saved. The default is to
#' create a file entitled **GenBanK_seqs**.
#'
#' @return A list of DNA sequences made of vectors of class "DNAbin", or of
#' single characters (if as.character = TRUE) with two attributes (species and
#' description).
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(Harpalyce_bayes_tree)
#'
#' mineSeq(inputdf = NULL,
#'         gb.colnames = NULL,
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         dir = "RESULTS_mineSeq",
#'         filename = "GenBanK_seqs")
#'}
#'
#' @importFrom ape read.GenBank write.dna
#'
#' @export
#'

mineSeq <- function(inputdf = NULL,
                    gb.colnames = NULL,
                    as.character = FALSE,
                    verbose = TRUE,
                    save = TRUE,
                    dir = "RESULTS_mineSeq",
                    filename = "GenBanK_seqs") {

  # Create folder to save mined GenBank seqs
  filenames <- vector()
  for (i in seq_along(gb.colnames)) {
    if (!grepl(gb.colnames[i], filename)) {
      filenames[i] <- paste0(filename, "_", gb.colnames[i])
    }
  }
  foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  if (!dir.exists(foldername)) {
    dir.create(foldername)
  }

  # Adjusting vouchers and species columns
  inputdf$Species <- gsub("[.]|(^\\s){1,}|(\\s$){1,}", "", inputdf$Species)
  inputdf$Species <- gsub("(\\s){1,}", "_", inputdf$Species)

  tf <- is.na(inputdf$Voucher)
  if (any(tf)) {
    inputdf$Voucher[tf] <- "Unvouchered"
  }
  inputdf$Voucher <- gsub("\\s[(].*|[:].*", "", inputdf$Voucher)
  inputdf$Voucher <- gsub("s[.]n[.]", "SN", inputdf$Voucher)
  inputdf$Voucher <- gsub("[/]|\\s|[-]", "", inputdf$Voucher)

  # Downloading seqs
  seqs <- list()
  for (i in seq_along(gb.colnames)) {
    if (verbose) {
      message(paste0("Mining GenBanK sequences for '", gb.colnames[i], "'..."))
    }
    seqs[[i]] <- ape::read.GenBank(access.nb = na.omit(inputdf[[gb.colnames[i]]]),
                                   species.names = T,
                                   as.character = as.character)

    # attr(seqs[[i]] , "species")
    # attr(seqs[[i]], "description")

    # Renaming accessions from GenBank
    for (l in seq_along(names(seqs[[i]]))){
      acc <- inputdf[[gb.colnames[i]]] %in% names(seqs[[i]])[l]
      names(seqs[[i]])[l] <- paste0(inputdf$Species[acc], "_",
                                    inputdf$Voucher[acc], "_",
                                    names(seqs[[i]])[l])
    }

    if (verbose) {
      message(paste0("GenBanK sequences for '", gb.colnames[i], "' fully downloaded!"))
    }

    # Exporting DNA matrix in fasta format
    if (save) {
      if (verbose) {
        message(paste0("GenBanK sequences for '", gb.colnames[i], "' saved on disk."))
      }
      ape::write.dna(seqs[[i]], paste0(foldername, "/", filenames[i], ".fasta"),
                     format = "fasta")
    }

  }

  names(seqs) <- gb.colnames

  return(seqs)
}

