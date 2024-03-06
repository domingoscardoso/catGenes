#' Read and download targeted loci from mitochondrial genomes in GenBank
#'
#' @author Domingos Cardoso
#'
#' @description A function built on the [genbankr](https://bioconductor.org/packages/release/bioc/html/genbankr.html)
#' package, designed to establish a connection with the GenBank database. This
#' function reads mitochondrial genomes using provided accession numbers, extracting
#' and formatting any specified targeted loci, and finally writing them in a
#' fasta file format.
#'
#' @usage
#' mineMitochondrion(genbank = NULL,
#'                   taxon = NULL,
#'                   voucher = NULL,
#'                   CDS = TRUE,
#'                   genes = NULL,
#'                   verbose = TRUE,
#'                   dir = "RESULTS_mineMitochondrion")
#'
#' @param genbank A vector comprising the GenBank accession numbers specifically
#' corresponding to the mitochondrial genome targeted for locus mining.
#'
#' @param taxon A vector containing the taxon name linked to the mitochondrial
#' genome. In the absence of this information, the function will default to the
#' existing nomenclature linked to the mitochondrial genome, as originally
#' provided in GenBank.
#'
#' @param voucher A vector containing relevant voucher information linked to the
#' mitochondrial genome. If this information is supplied, the function will promptly
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
#' named **RESULTS_mineMitochondrion** and the sequences will be saved within a
#' subfolder named after the current date.
#'
#' @return A fasta format file of DNA sequences saved on disk.
#'
#' @examples
#' \dontrun{
#' library(catGenes)
#'
#' mineMitochondrion(genbank = c("MN356196", "NC_008549"),
#'                   CDS = TRUE,
#'                   genes = c("COX1", "COX2", "ND4L"),
#'                   verbose = TRUE,
#'                   dir = "RESULTS_mineMitochondrion")
#'}
#'
#' @importFrom genbankr GBAccession readGenBank getSeq
#'
#' @export
#'

mineMitochondrion <- function(genbank = NULL,
                              taxon = NULL,
                              voucher = NULL,
                              CDS = TRUE,
                              genes = NULL,
                              verbose = TRUE,
                              dir = "RESULTS_mineMitochondrion") {

  # Create folder to save mined GenBank seqs from mitochondrial genomes
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

  list() -> mitochondrion_info -> mitochondrion_seq
  for (i in seq_along(genbank)) {

    # Skipping error in for-loop
    #https://stackoverflow.com/questions/14748557/skipping-error-in-for-loop
    tryCatch({

      gba <- genbankr::GBAccession(genbank[i])

      mitochondrion_info[[i]] <- genbankr::readGenBank(gba, partial = TRUE,
                                                   verbose = verbose)

      if (!is.null(mitochondrion_info[[i]])) {

        mitochondrion_seq[[i]] <- genbankr::getSeq(mitochondrion_info[[i]])

        for (j in seq_along(genes)) {

          if (CDS) {
            all_genes <- mitochondrion_info[[i]]@cds@elementMetadata@listData[["gene"]]
            tf <- all_genes %in% genes[j]
            if (any(tf)) {
              start <- mitochondrion_info[[i]]@cds@ranges@start[tf]
              end <- mitochondrion_info[[i]]@cds@ranges@width[tf]-1
            } else {
              stop(paste0("You must provide a correct gene name,\n",
                          "which may be any of the following:\n\n",
                          paste0(all_genes, collapse = ", ")),"\n\n",
                   "Find help also at DBOSLab-UFBA\n",
                   "(Domingos Cardoso; cardosobot@gmail.com)")
            }
          } else {
            all_genes <- mitochondrion_info[[i]]@other_features@elementMetadata@listData[["gene"]]
            tf <- all_genes %in% genes[j]
            if (any(tf)) {
              start <-  mitochondrion_info[[i]]@other_features@ranges@start[tf]
              end <- mitochondrion_info[[i]]@other_features@ranges@width[tf]-1
            } else {
              stop(paste0("You must provide a correct gene name,\n",
                          "which may be any of the following:\n\n",
                          paste0(all_genes, collapse = ", ")),"\n\n",
                   "Find help also at DBOSLab-UFBA\n",
                   "(Domingos Cardoso; cardosobot@gmail.com)")
            }
          }

          # get specific range of seqs
          targeted_seq <- substr(mitochondrion_seq[[i]], start, start+end)

          taxon_temp <- gsub(" ", "_", names(mitochondrion_seq[[i]]))
          genbank_temp <- gsub("_", "", genbank[i])
          write(
            paste0(">",
                   ifelse(is.null(taxon), taxon_temp, taxon[i]),
                   "_",
                   ifelse(is.null(voucher), "", paste0(voucher[i], "_")),
                   genbank_temp, "\n",
                   seq_revcompl(targeted_seq)),
            paste0(foldername, "/", genes[j], "_from_mitochondrion.fasta"),
            append = TRUE
          )
          if (verbose) {
            message("'", paste0(genes[j], "' of ", names(mitochondrion_seq[[i]]),
                                " ", genbank_temp, " mitochondrion retrieved!"))
          }
        }
      }
    }, error = function(e){cat(paste0("ERROR retrieving ", genbank[i], ":"),
                               conditionMessage(e), "\n")})
  }
}


