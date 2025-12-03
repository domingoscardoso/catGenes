#' Convert any format of DNA alignment into another format
#'
#' @author Domingos Cardoso
#'
#' @description Convert one or multiple DNA alignments in FASTA, NEXUS or PHYLIP
#' format into each other.
#'
#' @usage
#' convertAlign(filepath = NULL,
#'              format = NULL,
#'              rmfiles = FALSE,
#'              dir = NULL)
#'
#' @param filepath Path to the directory where the DNA alignments are stored.
#'
#' @param format Define either "NEXUS", "FASTA" or "PHYLIP" format for writing
#' the resulting newly converted files.
#'
#' @param rmfiles Logical, if \code{TRUE}, the original input file(s) are removed
#' from the directory and only the newly converted files are kept.
#'
#' @param dir The path to the directory where the newly converted files should
#' be saved. If no directory is given, then the files will saved in a subfolder
#' named after the current date under the folder named **RESULTS_convertAlign**.
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(GenBank_accessions)
#'
#' todaydate <- format(Sys.time(), "%d%b%Y")
#' folder_mined_seqs <- paste0("RESULTS_mineSeq/", todaydate)
#' folder_aligned_seqs <- paste0("RESULTS_alignSeqs/", todaydate)
#'
#' mineSeq(inputdf = GenBank_accessions,
#'         gb.colnames = c("ITS", "matK"),
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         dir = "RESULTS_mineSeq",
#'         filename = "GenBanK_seqs")
#'
#' alignSeqs(filepath = folder_mined_seqs,
#'           method = "ClustalW",
#'           gapOpening = "default",
#'           format = "NEXUS",
#'           verbose = TRUE,
#'           dir = "RESULTS_alignSeqs")
#'
#' convertAlign(filepath = folder_aligned_seqs,
#'           format = "PHYLIP",
#'           rmfiles = FALSE,
#'           verbose = TRUE,
#'           dir = "RESULTS_convertAlign")
#'}
#'
#' @importFrom ape read.nexus.data
#' @importFrom stats setNames
#'
#' @export
#'
convertAlign <- function(filepath = NULL,
                         format = NULL,
                         rmfiles = FALSE,
                         verbose = TRUE,
                         dir = "RESULTS_convertAlign") {

  inputfiles <- list.files(filepath)
  inputfiles_names <- gsub("[.].*", "", inputfiles)

  # Create folder to save the converted DNA sequences
  foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  if (!dir.exists(foldername)) {
    dir.create(foldername)
  }

  for (i in seq_along(inputfiles)) {

    temp <- readLines(paste0(filepath, "/", inputfiles[i]))

    if (strsplit(temp[1], '')[[1]][1] == "#") {
      if (format == "NEXUS") {
        message(paste0("The input file ", inputfiles[i],
                       "\nis already in NEXUS format!\n"),
                "This file was not converted.\n\n")
        next
      }
      # Convert a NEXUS file ####
      filetoconvert <- ape::read.nexus.data(paste0(filepath, "/", inputfiles[i]))
      for (j in seq_along(filetoconvert)) {
        filetoconvert[[j]] <- paste(filetoconvert[[j]], collapse = "")
        filetoconvert[[j]] <- toupper(filetoconvert[[j]])
      }
      spp <- stats::setNames(data.frame(names(filetoconvert)), "species")
      seqs <- stats::setNames(data.frame(unlist(filetoconvert, use.names = FALSE)), "sequence")
      convertedfile <- cbind(spp, seqs)

    } else if (strsplit(temp[1], '')[[1]][1] == ">") {
      if (format == "FASTA") {
        message(paste0("The input file ", inputfiles[i],
                       "\nis already in FASTA format!\n"),
                "This file was not converted.\n\n")
        next
      }
      # Convert a FASTA file ####
      filetoconvert <- readLines(paste0(filepath, "/", inputfiles[i]))
      tf <- grepl(">", filetoconvert)
      convertedfile <- data.frame(species = gsub(">", "", filetoconvert[tf]),
                                  sequence = filetoconvert[!tf])

    } else if (strsplit(temp[1], '')[[1]][1] != ">" &
               strsplit(temp[1], '')[[1]][1] != "#") {
      if (format == "PHYLIP") {
        message(paste0("The input file ", inputfiles[i],
                       "\nis already in PHYLIP format!\n"),
                "This file was not converted.\n\n")
        next
      }
      # Convert a PHYLIP file ####
      filetoconvert <- readLines(paste0(filepath, "/", inputfiles[i]))
      filetoconvert <- filetoconvert[-1]
      convertedfile <- data.frame(species = gsub("\\s.*", "", filetoconvert),
                                  sequence = gsub(".*\\s", "", filetoconvert))
    }

    if (format == "NEXUS") {
      nexusdframe(convertedfile, paste0(foldername, "/",
                                        inputfiles_names[i],".nex"))
      # Remove original file from the directory
      if (rmfiles) {
        unlink(paste0(filepath, "/", inputfiles[i]))
      }

    } else if (format == "PHYLIP") {
      phylipdframe(convertedfile,
                   paste0(foldername, "/", inputfiles_names[i],".phy"))
      # Remove original file from the directory
      if (rmfiles) {
        unlink(paste0(filepath, "/", inputfiles[i]))
      }

    } else if (format == "FASTA") {
      fastadframe(convertedfile, paste0(foldername, "/",
                                        inputfiles_names[i],".fasta"))
      # Remove original file from the directory
      if (rmfiles) {
        unlink(paste0(filepath, "/", inputfiles[i]))
      }
    }

  }

}

