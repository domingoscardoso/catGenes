#' Automated multiple sequence alignment
#'
#' @author Domingos Cardoso
#'
#' @description Perform automated multiple sequence alignment with
#' [msa](https://bioconductor.org/packages/release/bioc/html/msa.html) package
#' based either on [ClustalW](https://doi.org/10.1093/bioinformatics/btm404)
#' or [Muscle](https://doi.org/10.1186/1471-2105-5-113) algorithms. The function
#' uses one or multiple FASTA-formatted files to perform alignments and may save
#' the aligned sequences in FASTA, NEXUS or PHYLIP format.
#'
#' @usage
#' alignSeqs(filepath = GenBank_accessions,
#'           method = NULL,
#'           gapOpening = "default",
#'           format = "NEXUS",
#'           verbose = TRUE,
#'           dir = "RESULTS_alignSeqs",
#'           filename = "aligned_seqs")
#'
#' @param filepath Path to the directory where the FASTA-formatted DNA alignments
#' are stored.
#'
#' @param method Specifies the multiple sequence alignment to be used. Currently,
#' "ClustalW" and "Muscle" are supported.
#'
#' @param gapOpening Gap opening penalty; the defaults are specific to the
#' algorithm (see \code{\link{msaClustalW}} and \code{\link{msaMuscle)}}. Note
#' that the sign of this parameter is ignored. The sign is automatically adjusted
#' such that the called algorithm penalizes gaps instead of rewarding them.
#'
#' @param format Define either "NEXUS", "FASTA" or "PHYLIP" for writing the
#' resulting aligned DNA sequences in such formats. The default is to save the
#' aligned sequences in a NEXUS-formatted file.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the multiple sequence alignment will not be printed in the console in full.
#'
#' @param dir The path to the directory where the mined DNA sequences
#' in a fasta format file will be saved provided that the argument \code{save}
#' is set up in \code{TRUE}. The default is to create a directory named
#' **RESULTS_alignSeqs** and the sequences will be saved within a subfolder named
#' after the current date.
#'
#' @param filename A name or a vector of names of the output file(s) to be saved.
#' The default is to create output file(s) named based on the original name of the
#' input file(s) but also including an identifier suffix "aligned".
#'
#' @seealso \code{\link{mineSeq}}
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' data(GenBank_accessions)
#'
# 'dir <- "RESULTS_mineSeq"
# 'todaydate <- format(Sys.time(), "%d%b%Y")
#' folder_name_mined_seqs <- paste0(dir, "/", todaydate)
#'
#' mineSeq(inputdf = GenBank_accessions,
#'         gb.colnames = c("ETS", "ITS", "matK", "petBpetD", "trnTF", "Xdh"),
#'         as.character = FALSE,
#'         verbose = TRUE,
#'         save = TRUE,
#'         dir = dir,
#'         filename = "GenBanK_seqs")
#'
#' alignSeqs(filepath = folder_name_mined_seqs,
#'           method = "ClustalW",
#'           gapOpening = "default",
#'           format = "NEXUS",
#'           verbose = TRUE,
#'           dir = "RESULTS_alignSeqs",
#'           filename = "aligned_seqs")
#'}
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom msa msa msaConvert
#'
#' @export
#'

alignSeqs <- function(filepath = NULL,
                      method = NULL,
                      gapOpening = "default",
                      format = "NEXUS",
                      verbose = TRUE,
                      dir = "RESULTS_alignSeqs",
                      filename = NULL) {

  fasta_files <- list.files(filepath)

  if (!is.null(filename)) {
    name_files <- filename
  } else {
    name_files <- gsub(".*_|[.]fasta", "", fasta_files)
  }

  # Create folder to save the aligned DNA sequences
  foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  if (!dir.exists(foldername)) {
    dir.create(foldername)
  }

  # Align all DNA matrices automatically and save them in the directory
  for (i in seq_along(fasta_files)) {

    mySequences <- suppressWarnings(Biostrings::readDNAStringSet(paste0(filepath,
                                                       "/", fasta_files[i])))

    myAlignment <- msa::msa(mySequences,
                            method = method,
                            gapOpening = ifelse(gapOpening == "default",
                                                "default", gapOpening))
    myAlignment <- msa::msaConvert(myAlignment, type = "seqinr::alignment")

    myAlignment <- data.frame(species = myAlignment[["nam"]],
                              sequence = myAlignment[["seq"]])

    if (format == "NEXUS") {
      nexusdframe(myAlignment,
                  file = ifelse(!is.null(filename),
                                paste0(foldername, "/", name_files[i], ".nex"),
                                paste0(foldername, "/", name_files[i], "_aligned.nex")))
    } else if (format == "FASTA") {
      fastadframe(myAlignment,
                  file = ifelse(!is.null(filename),
                                paste0(foldername, "/", name_files[i], ".fasta"),
                                paste0(foldername, "/", name_files[i], "_aligned.fasta")))
    } else if (format == "PHYLIP") {
      phylipdframe(myAlignment,
                   file = ifelse(!is.null(filename),
                                 paste0(foldername, "/", name_files[i], ".phy"),
                                 paste0(foldername, "/", name_files[i], "_aligned.phy")))
    }

    message(paste0("Aligned ", name_files[i], " sequences already saved on disk!"))
  }

}
