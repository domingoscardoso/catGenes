#' Mine DNA sequences from GenBank using taxonomic queries
#'
#' @author Domingos Cardoso
#'
#' @description A comprehensive function built on the [rentrez](https://docs.ropensci.org/rentrez/)
#' package to search, download, and process DNA sequences from GenBank using
#' taxonomic query terms. This function performs multiple operations including
#' taxon name cleaning, voucher information extraction, duplicate removal,
#' and plastome separation. It iss particularly useful for phylogenetic studies
#' requiring large-scale sequence retrieval with consistent naming conventions.
#'
#' @usage
#' mineTaxa(term = NULL,
#'          db = "nucleotide",
#'          filename = "mined_seqs_by_taxon.fasta",
#'          clean_taxa = TRUE,
#'          add_voucher = TRUE,
#'          original_query = FALSE,
#'          plastome_apart = TRUE,
#'          rm_duplicated = TRUE,
#'          retmax = 4000,
#'          save = TRUE,
#'          verbose = TRUE,
#'          dir = "RESULTS_mineTaxa")
#'
#' @param term Character string specifying the search query for GenBank.
#' Should follow NCBI Entrez search syntax (e.g., "Solanaceae\\[Organism\\] AND matk",
#' "Arabidopsis\\[Organism\\] AND rbcL\\[Gene\\]"). See NCBI's
#' [Entrez Help](https://www.ncbi.nlm.nih.gov/books/NBK3837/) for detailed syntax.
#'
#' @param db Character string specifying the NCBI database to search.
#' Default is "nucleotide". Other options include "protein", "popset", etc.
#' See \code{\link[rentrez]{entrez_search}} for details.
#'
#' @param filename Character string for the output FASTA file name.
#' Default is "mined_seqs_by_taxon.fasta".
#'
#' @param clean_taxa Logical. If \code{TRUE} (default), taxon names are cleaned
#' and standardized (removes authors, subsp./var. indicators, special characters).
#' If \code{FALSE}, original GenBank names are preserved.
#'
#' @param add_voucher Logical. If \code{TRUE} (default), voucher specimen information
#' (when available) is extracted and appended to sequence names in the format:
#' Genus_species_voucher_accession. If \code{FALSE}, voucher information is omitted.
#'
#' @param original_query Logical. If \code{TRUE}, a copy of the original unprocessed
#' FASTA file is saved with "_ORIGINAL_QUERY" suffix. If \code{FALSE} (default),
#' only the cleaned file is kept. Ignored if \code{clean_taxa = FALSE}.
#'
#' @param plastome_apart Logical. If \code{TRUE} (default), complete plastome sequences
#' are separated into a separate file with "PLASTOMES_" prefix. Useful for
#' distinguishing whole plastomes from individual gene sequences.
#'
#' @param rm_duplicated Logical. If \code{TRUE} (default), removes duplicate sequences
#' from the same species, keeping only the longest sequence. Applied separately
#' to plastomes and other sequences.
#'
#' @param retmax Numeric. Maximum number of sequences to retrieve from GenBank.
#' Default is 4000. Increase for large queries, but note NCBI rate limits.
#'
#' @param save Logical. If \code{TRUE} (default), results are saved to disk in
#' the specified directory. If \code{FALSE}, results are returned as R objects
#' without saving.
#'
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed
#' to the console. If \code{FALSE}, function runs silently.
#'
#' @param dir Character string specifying the directory path for saving results.
#' Default is "RESULTS_mineTaxa". A subdirectory with current date will be created
#' within this directory.
#'
#' @return If \code{save = TRUE}, returns an invisible list of DNA sequences in
#' \code{DNAbin} format (if \code{clean_taxa = TRUE}) or raw FASTA text (if
#' \code{clean_taxa = FALSE}). If \code{save = FALSE}, returns the processed
#' sequences as an R object. Files are written to disk in either case when
#' \code{save = TRUE}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Connects to GenBank using the specified search term
#'   \item Downloads sequences (up to \code{retmax} records)
#'   \item If \code{clean_taxa = TRUE}:
#'   \itemize{
#'     \item Extracts and cleans organism names from each accession
#'     \item Optionally adds voucher information
#'     \item Standardizes naming format: Genus_species_voucher_accession
#'     \item Separates plastomes if \code{plastome_apart = TRUE}
#'     \item Removes duplicates if \code{rm_duplicated = TRUE}
#'   }
#'   \item Saves results to specified directory with date-stamped subfolder
#' }
#'
#' Taxon name cleaning includes:
#' \itemize{
#'   \item Removing periods and hyphens
#'   \item Converting "var.", "subsp.", "f." to simple spaces
#'   \item Converting spaces to underscores
#'   \item Trimming "sp." designations
#'   \item Removing trailing underscores
#'   \item Handling hybrid cultivars
#' }
#'
#' @note
#' \itemize{
#'   \item Large queries may take considerable time due to NCBI rate limiting
#'   \item Voucher extraction may not work for all GenBank entries (depends on annotation)
#'   \item For very large datasets (>10,000 sequences), consider splitting queries
#'   \item Always verify taxon name cleaning results for your specific taxonomic group
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic search for matK sequences in Solanaceae (excluding Solanum)
#' result1 <- mineTaxa(
#'   term = "(Solanaceae[Organism] AND matk NOT Solanum[Organism])",
#'   filename = "Solanaceae_outSolanum.fasta"
#' )
#'
#' # Example 2: Search for rbcL sequences in Fabaceae with custom settings
#' result2 <- mineTaxa(
#'   term = "Fabaceae[Organism] AND rbcL[Gene]",
#'   filename = "Fabaceae_rbcL.fasta",
#'   clean_taxa = TRUE,
#'   add_voucher = FALSE,          # Don't add voucher info
#'   original_query = TRUE,        # Keep original unprocessed file
#'   plastome_apart = FALSE,       # Don't separate plastomes
#'   rm_duplicated = FALSE,        # Keep all sequences
#'   retmax = 1000,
#'   verbose = TRUE
#' )
#'
#' # Example 3: Protein database search
#' result3 <- mineTaxa(
#'   term = "Rubisco[Protein] AND plants[Organism]",
#'   db = "protein",
#'   filename = "plant_Rubisco.fasta",
#'   clean_taxa = FALSE           # Don't clean names for proteins
#' )
#'
#' # Example 4: Return results without saving to disk
#' result4 <- mineTaxa(
#'   term = "Arabidopsis[Organism] AND chloroplast[Title]",
#'   save = FALSE,
#'   verbose = FALSE
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[rentrez]{entrez_search}} for NCBI search syntax
#'   \item \code{\link[ape]{read.FASTA}} for FASTA file reading
#'   \item \code{\link{minePlastome}} for targeted plastome mining
#' }
#'
#' @importFrom rentrez entrez_search entrez_fetch entrez_summary
#' @importFrom ape read.FASTA write.FASTA
#' @importFrom utils head tail
#'
#' @export
#'
mineTaxa <- function(term = NULL,
                     db = "nucleotide",
                     filename = "mined_seqs_by_taxon.fasta",
                     clean_taxa = TRUE,
                     add_voucher = TRUE,
                     original_query = FALSE,
                     plastome_apart = TRUE,
                     rm_duplicated = TRUE,
                     retmax = 4000,
                     save = TRUE,
                     verbose = TRUE,
                     dir = "RESULTS_mineTaxa") {

  # Validate inputs
  if (is.null(term) || term == "") {
    stop("Search term cannot be empty. Please provide a valid NCBI query.")
  }

  if (!is.character(term)) {
    stop("Search term must be a character string.")
  }

  if (!is.numeric(retmax) || retmax <= 0) {
    stop("retmax must be a positive number.")
  }

  # Create folder to save mined taxa in GenBank
  if (save) {
    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername)
    }
  }

  if (verbose) {
    message(paste0("Searching GenBank for: '", term, "'"))
  }

  # Perform GenBank search
  search <- rentrez::entrez_search(db = db,
                                   term = term,
                                   use_history = TRUE)

  if (search$count == 0) {
    stop("No sequences found for the specified search term.")
  }

  if (verbose) {
    message(paste0("Found ", search$count, " records. Downloading sequences..."))
  }

  # Fetch sequences
  base_query <- rentrez::entrez_fetch(db = db,
                                      web_history = search$web_history,
                                      retmax = retmax,
                                      rettype = "fasta")

  # just to see the first 1000 characters (for debugging)
  # cat(strwrap(substr(base_query, 1, 10000)), sep="\n")

  # Save original query if requested
  if (save) {
    if (verbose) {
      message(paste0("DNA successfully retrieved!\n"),
              paste0("FASTA file '", filename, "' saved within ", foldername))
    }
    write(base_query, file = paste0(foldername, "/", filename))
  }

  # Process sequences if cleaning is requested
  if (clean_taxa) {

    if (verbose) {
      message("Processing and cleaning taxon names...")
    }

    # Read original downloaded fasta back into R
    fasta_query <- ape::read.FASTA(paste0(foldername, "/", filename))
    full_labels <- names(fasta_query)

    # Get accession numbers and original taxon names
    accessions <- gsub("[.].*", "", names(fasta_query))
    metadata <- list()

    # Process each accession
    for (i in seq_along(accessions)) {

      metadata[[i]] <- .get_metadata(acc = accessions[i],
                                     db = db,
                                     verbose = verbose)

      if (verbose) {
        message(paste0(i,"/", length(accessions),
                       " '", metadata[[i]]$organism,"'",
                       " extracted from accession ",
                       accessions[i]))
      }

      # Replace original taxon names by the cleaned names + accessions
      if (add_voucher) {
        if (metadata[[i]]$voucher != "Unvouchered" && verbose) {
          message(paste0("Voucher '", metadata[[i]]$voucher, "' added!"))
        }
        names(fasta_query)[i] <- paste0(metadata[[i]]$organism, "_",
                                        metadata[[i]]$voucher, "_",
                                        accessions[i])
      } else {
        names(fasta_query)[i] <- paste0(metadata[[i]]$organism, "_",
                                        accessions[i])
      }
    }

    # Handle original query file
    if (original_query && save) {
      file.copy(from = paste0(foldername, "/", filename),
                to = gsub("[.]fasta", "_ORIGINAL_QUERY.fasta",
                          paste0(foldername, "/", filename)))
      if (verbose) {
        message("Original query saved as: ",
                gsub("[.]fasta", "_ORIGINAL_QUERY.fasta", filename))
      }
    } else if (save) {
      unlink(paste0(foldername, "/", filename))
    }

    # Save processed sequences
    if (save) {
      fasta_query_temp <- fasta_query

      # Separate plastomes if requested
      if (plastome_apart) {
        tf <- grepl("complete genome", full_labels, ignore.case = TRUE)
        if (any(tf)) {
          temp <- names(fasta_query_temp)[tf]
          plastomes <- fasta_query_temp[names(fasta_query_temp) %in% temp]
          fasta_query_temp <- fasta_query_temp[!names(fasta_query_temp) %in% temp]

          # Remove duplicate plastomes if requested
          if (rm_duplicated && length(plastomes) > 0) {
            plastomes <- .rm_duplicated_seqs(plastomes)
            if (verbose) {
              message(paste0(length(plastomes),
                             " unique plastome sequences kept after removing duplicates."))
            }
          }

          # Save plastomes separately in GenBank format
          if (length(plastomes) > 0) {
            # Create subdirectory for plastome GB files
            plastome_dir <- paste0(foldername, "/PLASTOME_GB_FILES")
            if (!dir.exists(plastome_dir)) {
              dir.create(plastome_dir)
            }

            if (verbose) {
              message("Downloading plastomes in GenBank format...")
            }

            for (i in seq_along(plastomes)) {
              # Extract accession from the cleaned name
              # Format is: Genus_species_voucher_accession or Genus_species_accession
              full_name <- names(plastomes)[i]

              # Extract accession - it's the last component after the last underscore
              accession <- sub(".*_([A-Z]+[0-9]+\\.[0-9]+)$", "\\1", full_name)

              # If pattern doesn't match, try alternative pattern (without version number)
              if (accession == full_name) {
                accession <- sub(".*_([A-Z]+[0-9]+)$", "\\1", full_name)
              }

              # If still no match, try to extract from original labels
              if (accession == full_name) {
                # Get the original label before cleaning
                original_idx <- which(grepl(accession, full_labels, fixed = TRUE))[1]
                if (!is.na(original_idx)) {
                  # Extract accession from original GenBank header
                  accession <- gsub("^.*gb\\|([A-Z]+[0-9]+\\.[0-9]+)\\|.*", "\\1",
                                    full_labels[original_idx])
                  if (accession == full_labels[original_idx]) {
                    accession <- gsub("^.*gb\\|([A-Z]+[0-9]+)\\|.*", "\\1",
                                      full_labels[original_idx])
                  }
                }
              }

              if (!is.na(accession) && accession != full_name) {
                tryCatch({
                  # Download GenBank format file
                  filename_gb <- paste0(accession, "_plastome.gb")
                  record <- rentrez::entrez_fetch(db = "nucleotide",
                                                  id = accession,
                                                  rettype = "gb",
                                                  retmode = "text")
                  write(record, file = paste0(plastome_dir, "/", filename_gb))

                  if (verbose && i %% 10 == 0) {
                    message(paste0("Downloaded ", i, "/", length(plastomes),
                                   " plastomes in GenBank format"))
                  }

                  # Be polite to NCBI servers
                  Sys.sleep(0.2)

                }, error = function(e) {
                  if (verbose) {
                    message(paste0("Failed to download GenBank file for ", accession,
                                   ": ", conditionMessage(e)))
                  }
                })
              } else {
                if (verbose) {
                  message(paste0("Could not extract valid accession from: ", full_name))
                }
              }
            }
            if (verbose) {
              message(paste0("Plastome sequences saved:"))
              message(paste0("  - GenBank (.gb) files in: ", plastome_dir))
              message(paste0("Total plastomes downloaded: ", length(plastomes)))
            }
          }
        }
      }

      # Remove duplicate sequences if requested
      if (rm_duplicated && length(fasta_query_temp) > 0) {
        original_count <- length(fasta_query_temp)
        fasta_query_temp <- .rm_duplicated_seqs(fasta_query_temp)
        removed_count <- original_count - length(fasta_query_temp)
        if (verbose && removed_count > 0) {
          message(paste0("Removed ", removed_count,
                         " duplicate sequences, keeping longest for each species."))
        }
      }

      # Save main sequence file
      ape::write.FASTA(fasta_query_temp,
                       paste0(foldername, "/", filename))

      if (verbose) {
        message("\nProcessing complete!")
        message(paste0("Total sequences processed: ", length(fasta_query)))
        if (exists("plastomes")) {
          message(paste0("Plastomes separated: ", length(plastomes)))
        }
        message(paste0("Final sequences saved: ", length(fasta_query_temp)))
        message(paste0("Results saved in: ", foldername))
      }

      return(invisible(fasta_query))
    } else {
      return(fasta_query)
    }

  } else {
    # Return unprocessed results
    if (save && verbose) {
      message("Sequences saved without cleaning.")
    }
    if (save) {
      return(invisible(base_query))
    } else {
      return(base_query)
    }
  }
}


#-------------------------------------------------------------------------------
# Internal function to get metadata from GenBank accessions
.get_metadata <- function(acc,
                          db = db,
                          verbose = verbose) {

  # Search for the accession
  search_result <- rentrez::entrez_search(db = db,
                                          term = acc)

  # Fetch full record (summary) using the retrieved ID
  summary_record <- rentrez::entrez_summary(db = db,
                                            id = search_result$ids[[1]])

  # Fetch full GenBank record for voucher extraction
  gb_record <- rentrez::entrez_fetch(db = db,
                                     id = acc,
                                     rettype = "gb",
                                     retmode = "text")

  # Extract voucher information using regex
  voucher_info <- sub(".*\\s+/specimen_voucher=\"([^\"]+)\".*", "\\1", gb_record)

  # Check if voucher extraction failed
  if (grepl("LOCUS ", voucher_info) || voucher_info == gb_record) {
    voucher_info <- "Unvouchered"
  } else {
    # Clean voucher string
    voucher_info <- gsub(".*[:]", "", voucher_info)
    voucher_info <- gsub("\\set al[.]|\\s[&]\\sal[.]", "", voucher_info)
    voucher_info <- sub(".*?\\b(\\w+)\\b.*?(\\d+).*", "\\1\\2", voucher_info)
    voucher_info <- gsub("\\s[(].*|[/]|[.]|[<]|[>]|[:]|[-]|[,]", "", voucher_info)
    voucher_info <- gsub("\\s", "", voucher_info)
    voucher_info <- gsub("personalcollection", "", voucher_info)

    # Final check
    if (voucher_info == "" || nchar(voucher_info) < 2) {
      voucher_info <- "Unvouchered"
    }
  }

  # Create metadata dataframe
  metadata <- data.frame(organism = summary_record$organism,
                         voucher = voucher_info)

  # Clean the taxon names
  if (verbose) {
    message(paste0("Cleaning taxon name for ", acc, "... "))
  }

  # Apply cleaning steps
  metadata$organism <- gsub("[.]|[-]", "", metadata$organism)
  metadata$organism <- gsub("\\svar\\s|\\ssubsp\\s|\\sf\\s", " ", metadata$organism)
  metadata$organism <- gsub("\\s", "_", metadata$organism)
  metadata$organism <- sub("(_sp_).*", "\\1", metadata$organism)
  metadata$organism <- gsub("_$", "", metadata$organism)
  metadata$organism <- gsub("_x_", "_", metadata$organism)
  metadata$organism <- gsub("hybrid_cultivar", "sp", metadata$organism)

  return(metadata)
}


#-------------------------------------------------------------------------------
# Internal function to remove duplicated sequences by keeping just the largest one
.rm_duplicated_seqs <- function(seqs) {

  if (length(seqs) == 0) {
    return(seqs)
  }

  # Extract species names (Genus_species or Genus_species_subspecies)
  temp <- sub("(^[A-Za-z]+_[a-z]+(?:_[a-z]+)?)_.*$", "\\1", names(seqs))

  # Check for duplicates
  tf <- duplicated(temp) | duplicated(temp, fromLast = TRUE)

  if (any(tf)) {
    # Calculate sequence lengths
    sequence_lengths <- sapply(seqs, length)
    names(sequence_lengths) <- temp

    # Get species names
    species_names <- names(sequence_lengths)

    # Compute the max value per species
    max_values <- tapply(sequence_lengths, species_names, max)

    # Identify which entries are the maximum per species
    is_max <- sequence_lengths == max_values[species_names]

    # Ensure only the first occurrence of max values is kept
    tf <- !duplicated(paste(species_names, sequence_lengths)) & is_max

    # Filter sequences
    keep_names <- names(seqs)[tf]
    seqs <- seqs[names(seqs) %in% keep_names]
  }

  return(seqs)
}
