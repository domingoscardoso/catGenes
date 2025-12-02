#' Read and download targeted loci from plastome sequences in GenBank
#'
#' @author Domingos Cardoso
#'
#' @description A function built on the [rentrez](https://docs.ropensci.org/rentrez/)
#' and [geneviewer](https://nvelden.github.io/geneviewer/)
#' packages, designed to establish a connection with the GenBank database, donwload,
#' and parse plastomes. This function downloads plastome sequences using provided
#' accession numbers, extracting and formatting any specified targeted loci, and
#' finally writing them in a fasta file format.
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
#' @param rm_gb_files Logical, if \code{TRUE}, the downloaded .gb files from
#' GenBank will be removed from the directory after extracting the targeted loci.
#' The default is \code{FALSE}, keeping the original .gb files.
#'
#' @param verbose Logical, if \code{FALSE}, a message showing each step during
#' the GenBank search will not be printed in the console in full.
#'
#' @param dir The path to the directory where the mined DNA sequences
#' in a fasta format file will be saved. The default is to create a directory
#' named **RESULTS_minePlastome** and the sequences will be saved within a
#' subfolder named after the current date.
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
#'              rm_gb_files = FALSE,
#'              verbose = TRUE,
#'              dir = "RESULTS_minePlastome")
#'}
#'
#' @importFrom rentrez entrez_fetch
#' @importFrom geneviewer read_gbk gbk_features_to_df
#' @importFrom flora remove.authors
#'
#' @export
#'

minePlastome <- function(genbank = NULL,
                         taxon = NULL,
                         voucher = NULL,
                         CDS = TRUE,
                         genes = NULL,
                         rm_gb_files = FALSE,
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

      filename <- paste0(genbank[i], "_plastome.gb")

      record <- rentrez::entrez_fetch(db = "nucleotide",
                                      id = genbank[i],
                                      rettype = "gb",
                                      retmode = "text")
      write(record, file = paste0(foldername, "/", filename))

      # change the path to where you have saved the file
      file_path <- paste0(foldername, "/", filename)
      gbk <- geneviewer::read_gbk(file_path)

      cds_df <- geneviewer::gbk_features_to_df(
        gbk,
        feature = "CDS",
        keys = c("region", "gene", "protein_id", "gene_kind", "product"),
        process_region = TRUE
      )

      gene_df <- geneviewer::gbk_features_to_df(
        gbk,
        feature = "gene",
        keys = c("region", "gene", "gene_kind"),
        process_region = TRUE
      )

        for (j in seq_along(genes)) {
          if (CDS) {
            all_genes <- cds_df$gene
            tf <- all_genes %in% genes[j]
            if (any(tf)) {
              start <- cds_df$start[tf]
              end <- cds_df$end[tf]
              strand <- cds_df$strand[tf]
            } else {
              stop(paste0("You must provide a correct gene name,\n",
                          "which may be any of the following:\n\n",
                          paste0(all_genes, collapse = ", ")),"\n\n",
                   "Find help also with:\n",
                   "(Domingos Cardoso (JBRJ; cardosobot@gmail.com)")
            }
          } else {
            all_genes <- gene_df$gene
            tf <- all_genes %in% genes[j]
            if (any(tf)) {
              start <- gene_df$start[tf]
              end <- gene_df$end[tf]
              strand <- gene_df$strand[tf]
            } else {
              stop(paste0("You must provide a correct gene name,\n",
                          "which may be any of the following:\n\n",
                          paste0(all_genes, collapse = ", ")),"\n\n",
                   "Find help also with:\n",
                   "(Domingos Cardoso (JBRJ; cardosobot@gmail.com)")
            }
          }

          # get specific range of seqs
          targeted_seq <- substr(gbk[[1]][["ORIGIN"]], start, end)

          taxon_temp <- gsub(" ", "_", gbk[[1]][["FEATURES"]][["source"]][[1]][["organism"]])
          genbank_temp <- gsub("_", "", genbank[i])
          write(
            paste0(">",
                   ifelse(is.null(taxon), taxon_temp, taxon[i]),
                   "_",
                   ifelse(is.null(voucher), "", paste0(voucher[i], "_")),
                   genbank_temp, "\n",
                   .seq_revcompl(targeted_seq, strand)),
            paste0(foldername, "/", genes[j], "_from_plastome.fasta"),
            append = TRUE
          )
          if (verbose) {
            message("'", paste0(genes[j], "' of ", taxon_temp,
                                " ", genbank_temp, " plastome retrieved!"))
          }
        }

      # Remove .gb files if rm_gb_files is TRUE
      if (rm_gb_files) {
        file.remove(file_path)
        if (verbose) {
          message("Removed .gb file: ", filename)
        }
      }

    }, error = function(e){cat(paste0("ERROR retrieving ", genbank[i], ":"),
                               conditionMessage(e), "\n")})
  }

  # Final message about file removal
  if (rm_gb_files && verbose) {
    message("\nAll .gb files have been removed from directory.")
    message("Only extracted loci in FASTA format are kept.")
  } else if (!rm_gb_files && verbose) {
    message("\nAll .gb files are kept in directory for reference.")
  }

}
