#' Read and download targeted loci from plastome sequences in GenBank
#'
#' @author Domingos Cardoso
#'
#' @description A function built on the [genbankr](https://bioconductor.org/packages/release/bioc/html/genbankr.html)
#' package, designed to establish a connection with the GenBank database. This
#' function reads plastome sequences using provided accession numbers, extracting
#' and formatting any specified targeted loci, and finally writing them in a
#' fasta file format.
#'
#' @usage
#' minePlastome(genbank = NULL,
#'              taxon = NULL,
#'              voucher = NULL,
#'              CDS = TRUE,
#'              genes = NULL,
#'              verbose = TRUE,
#'              dir = "RESULTS_minePlastome")
#'
#' @param genbank A vector comprising the GenBank accession numbers specifically
#' corresponding to the plastome sequence targeted for locus mining.
#'
#' @param taxon A vector containing the taxon name linked to the plastome sequence.
#' In the absence of this information, the function will default to the existing
#' nomenclature linked to the plastome, as originally provided in GenBank.
#'
#' @param voucher A vector containing relevant voucher information linked to the
#' plastome sequence. If this information is supplied, the function will promptly
#' append it immediately following the taxon name of the downloaded targeted
#' sequence.
#'
#' @param CDS a logical controlling whether the targeted loci are protein coding
#' genes, otherwise the function understands that entered gene names are e.g.
#' intron or intergenic spacer regions.
#'
#' @param genes A vector of one or more gene names as annotated in GenBank.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the GenBank search will not be printed in the console in full.
#'
#' @param dir Pathway to the computer's directory, where the mined DNA sequences
#' in a fasta format file will be saved. The default is to create a directory
#' named **RESULTS_minePlastome** and the sequences will be saved within a
#' subfolder named after the current date.
#'
#' @return A fasta format file of DNA sequences saved on disk.
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#' library(dplyr)
#'
#' data(GenBank_accessions)
#'
#' GenBank_plastomes <- GenBank_accessions %>%
#'   filter(!is.na(Plastome)) %>%
#'   select(c("Species", "Voucher", "Plastome"))
#'
#' minePlastome(genbank = GenBank_plastomes$Plastome,
#'              taxon = GenBank_plastomes$Species,
#'              voucher = GenBank_plastomes$Voucher,
#'              CDS = TRUE,
#'              genes = c("matK", "rbcL"),
#'              verbose = TRUE,
#'              dir = "RESULTS_minePlastome")
#'}
#'
#' @importFrom genbankr GBAccession readGenBank getSeq
#' @importFrom flora remove.authors
#'
#' @export
#'

minePlastome <- function(genbank = NULL,
                         taxon = NULL,
                         voucher = NULL,
                         CDS = TRUE,
                         genes = NULL,
                         verbose = TRUE,
                         dir = "RESULTS_minePlastome") {

  # Create folder to save mined GenBank seqs from plastomes
  foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  if (!dir.exists(foldername)) {
    dir.create(foldername)
  }

  # Adjusting vouchers and species columns
  temp <- .tax_voucher_adjust(inputdf = NULL,
                              taxon = taxon,
                              voucher = voucher,
                              genbank = genbank)
  taxon = temp[[1]]
  voucher = temp[[2]]
  genbank = temp[[3]]

  list() -> plastomes_info -> plastomes_seq
  for (i in seq_along(genbank)) {

    # Skipping error in for-loop
    #https://stackoverflow.com/questions/14748557/skipping-error-in-for-loop
    tryCatch({

      gba <- genbankr::GBAccession(genbank[i])

      plastomes_info[[i]] <- genbankr::readGenBank(gba, partial = TRUE,
                                                   verbose = verbose)

      if (!is.null(plastomes_info[[i]])) {

        plastomes_seq[[i]] <- genbankr::getSeq(plastomes_info[[i]])

        for (j in seq_along(genes)) {

          if (CDS) {
            all_genes <- plastomes_info[[i]]@cds@elementMetadata@listData[["gene"]]
            tf <- all_genes %in% genes[j]
            if (any(tf)) {
              start <- plastomes_info[[i]]@cds@ranges@start[tf]
              end <- plastomes_info[[i]]@cds@ranges@width[tf]-1
            } else {
              stop(paste0("You must provide a correct gene name,\n",
                          "which may be any of the following:\n\n",
                          paste0(all_genes, collapse = ", ")),"\n\n",
                   "Find help also at DBOSLab-UFBA\n",
                   "(Domingos Cardoso; cardosobot@gmail.com)")
            }
          } else {
            all_genes <- plastomes_info[[i]]@other_features@elementMetadata@listData[["gene"]]
            tf <- all_genes %in% genes[j]
            if (any(tf)) {
              start <-  plastomes_info[[i]]@other_features@ranges@start[tf]
              end <- plastomes_info[[i]]@other_features@ranges@width[tf]-1
            } else {
              stop(paste0("You must provide a correct gene name,\n",
                          "which may be any of the following:\n\n",
                          paste0(all_genes, collapse = ", ")),"\n\n",
                   "Find help also at DBOSLab-UFBA\n",
                   "(Domingos Cardoso; cardosobot@gmail.com)")
            }
          }

          # get specific range of seqs
          targeted_seq <- substr(plastomes_seq[[i]], start, start+end)

          taxon_temp <- gsub(" ", "_", names(plastomes_seq[[i]]))
          genbank_temp <- gsub("_", "", genbank[i])
          write(
            paste0(">",
                   ifelse(is.null(taxon), taxon_temp, taxon[i]),
                   "_",
                   ifelse(is.null(voucher), "", paste0(voucher[i], "_")),
                   genbank_temp, "\n",
                   seq_revcompl(targeted_seq)),
            paste0(foldername, "/", genes[j], "_from_plastome.fasta"),
            append = TRUE
          )
          if (verbose) {
            message("'", paste0(genes[j], "' of ", names(plastomes_seq[[i]]),
                                " ", genbank_temp, " plastome retrieved!"))
          }
        }
      }
    }, error = function(e){cat(paste0("ERROR retrieving ", genbank[i], ":"),
                               conditionMessage(e), "\n")})
  }
}
