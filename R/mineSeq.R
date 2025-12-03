#' Read and download DNA sequences from GenBank
#'
#' @author Domingos Cardoso
#'
#' @description An [ape](https://cran.r-project.org/web/packages/ape/index.html)-based
#' function to connect with the [GenBank](https://www.ncbi.nlm.nih.gov/genbank)
#' database, read nucleotide sequences using accession numbers, and write them
#' in a fasta format file.
#'
#' @usage
#' mineSeq(inputdf = NULL,
#'         gb.colnames = NULL,
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         filename = "GenBanK_seqs",
#'         dir = "RESULTS_mineSeq")
#'
#' @param inputdf A dataframe object containing the taxon names in a 'Species'
#' column, the voucher information in 'Voucher' column, and the GenBank accessions
#' for each genes in separate columns named by the corresponding gene. If the
#' columns 'Species' and 'Voucher' are not provided in the dataframe, then the
#' function will consider the taxonomy of the retrieved sequences as originally
#' available in GenBank.
#'
#' @param gb.colnames A vector with column names within the \code{inputdf}
#' dataframe corresponding to each gene, where the GenBank accession numbers are
#' listed.
#'
#' @param as.character A logical controlling whether to return the sequences as
#' an object of class "DNAbin" (the default).
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the GenBank search will not be printed in the console in full.
#'
#' @param save Logical, if \code{TRUE}, the mined sequences will be saved on disk.
#'
#' @param dir The path to the directory where the mined DNA sequences in a
#' fasta format file will be saved provided that the argument \code{save}
#' is set up in \code{TRUE}. The default is to create a directory named
#' **RESULTS_mineSeq** and the sequences will be saved within a subfolder named
#' after the current date.
#'
#' @param filename Name of the output file to be saved. The default is to create
#' a file entitled **GenBanK_seqs**.
#'
#' @return A list of DNA sequences made of vectors of class 'DNAbin', or of
#' single characters (if as.character = TRUE) with two attributes (species and
#' description).
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(GenBank_accessions)
#'
#' mineSeq(inputdf = GenBank_accessions,
#'         gb.colnames = c("ETS", "ITS", "matK", "petBpetD", "trnTF", "Xdh"),
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         filename = "GenBanK_seqs",
#'         dir = "RESULTS_mineSeq")
#'}
#'
#' @importFrom ape read.GenBank write.dna
#' @importFrom flora remove.authors
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#'
#' @export
#'

mineSeq <- function(inputdf = NULL,
                    gb.colnames = NULL,
                    as.character = FALSE,
                    verbose = TRUE,
                    save = TRUE,
                    filename = "GenBanK_seqs",
                    dir = "RESULTS_mineSeq") {

  # dir check
  dir <- .arg_check_dir(dir)

  # inputdf check
  .arg_check_inputdf(inputdf, gb.colnames)

  # Create folder to save mined GenBank seqs
  if (save) {
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
  }
  # Adjusting vouchers and species columns
  inputdf <- .tax_voucher_adjust(inputdf = inputdf,
                                 gb.colnames = gb.colnames)

  # Downloading seqs
  seqs <- list()
  for (i in seq_along(gb.colnames)) {
    if (verbose) {
      message(paste0("Mining GenBanK sequences for '", gb.colnames[i], "'..."))
    }
    seqs[[i]] <- ape::read.GenBank(access.nb = na.omit(inputdf[[gb.colnames[i]]]),
                                   species.names = T,
                                   as.character = as.character)
    # attr(seqs[[i]], "species")
    # attr(seqs[[i]], "description")

    # Renaming accessions from GenBank
    if ("Species" %in% names(inputdf) &
        "Voucher" %in% names(inputdf)) {
      for (l in seq_along(names(seqs[[i]]))){
        acc <- inputdf[[gb.colnames[i]]] %in% names(seqs[[i]])[l]
        if (!is.na(inputdf$Species[acc])) {
          names(seqs[[i]])[l] <- paste0(inputdf$Species[acc], "_",
                                        inputdf$Voucher[acc], "_",
                                        names(seqs[[i]])[l])
        } else {
          names(seqs[[i]])[l] <- paste0(attr(seqs[[i]], "species")[l], "_",
                                        inputdf$Voucher[acc], "_",
                                        names(seqs[[i]])[l])
        }
      }
    }
    if ("Species" %in% names(inputdf) &
        !"Voucher" %in% names(inputdf)) {
      for (l in seq_along(names(seqs[[i]]))){
        acc <- inputdf[[gb.colnames[i]]] %in% names(seqs[[i]])[l]
        if (!is.na(inputdf$Species[acc])) {
          names(seqs[[i]])[l] <- paste0(inputdf$Species[acc], "_",
                                        names(seqs[[i]])[l])
        } else {
          names(seqs[[i]])[l] <- paste0(attr(seqs[[i]], "species")[l], "_",
                                        names(seqs[[i]])[l])
        }
      }
    }
    if (!"Species" %in% names(inputdf) &
        !"Voucher" %in% names(inputdf)) {
      for (l in seq_along(names(seqs[[i]]))){
        acc <- inputdf[[gb.colnames[i]]] %in% names(seqs[[i]])[l]
        names(seqs[[i]])[l] <- paste0(attr(seqs[[i]], "species")[l], "_",
                                      names(seqs[[i]])[l])
      }
      names(seqs[[i]]) <- gsub("[.]", "", names(seqs[[i]]))
    }
    if (!"Species" %in% names(inputdf) &
        "Voucher" %in% names(inputdf)) {
      for (l in seq_along(names(seqs[[i]]))){
        acc <- inputdf[[gb.colnames[i]]] %in% names(seqs[[i]])[l]
        names(seqs[[i]])[l] <- paste0(attr(seqs[[i]], "species")[l], "_",
                                      inputdf$Voucher[acc], "_",
                                      names(seqs[[i]])[l])
      }
      names(seqs[[i]]) <- gsub("[.]", "", names(seqs[[i]]))
    }

    if (verbose) {
      message(paste0("GenBanK sequences for '", gb.colnames[i],
                     "' fully downloaded!"))
    }

    # Exporting DNA matrix in fasta format
    if (save) {
      if (verbose) {
        message(paste0("GenBanK sequences for '", gb.colnames[i],
                       "' saved on disk."))
      }
      ape::write.dna(seqs[[i]], paste0(foldername, "/", filenames[i], ".fasta"),
                     format = "fasta")
    }

  }

  names(seqs) <- gb.colnames

  return(seqs)
}
