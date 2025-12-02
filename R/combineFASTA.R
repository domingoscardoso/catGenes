#' Combine multiple FASTA files into a single file
#'
#' @author Domingos Cardoso
#'
#' @description This function reads multiple FASTA files and combines them into a
#' single FASTA file. It saves the combined sequences to a specified output directory.
#'
#' @usage
#' combineFASTA(input_files = NULL,
#'              output_file = "combined_sequences.fasta",
#'              save = TRUE,
#'              verbose = TRUE,
#'              dir = "RESULTS_combineFASTA")
#'
#' @param input_files Character vector specifying the paths to FASTA files to be
#' combined. This parameter is required - users must specify exactly which files
#' they want to combine.
#'
#' @param output_file Character string specifying the name of the output combined
#' FASTA file. Default is "combined_sequences.fasta".
#'
#' @param save Logical. If \code{TRUE} (default), the combined sequences are saved
#' to disk in the specified directory. If \code{FALSE}, the combined sequences are
#' only returned as an R object without saving.
#'
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed
#' to the console.
#'
#' @param dir Character string specifying the directory path for saving results.
#' Default is "RESULTS_combineFASTA". A subdirectory with current date will be
#' created within this directory when \code{save = TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sequences}: DNAbin object with all combined sequences
#'   \item \code{summary}: Data frame with summary statistics
#'   \item \code{output_path}: Path to the saved combined FASTA file (if saved)
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that all specified input files exist
#'   \item Creates an output directory with date stamp if \code{save = TRUE}
#'   \item Reads each specified FASTA file using \code{ape::read.FASTA()}
#'   \item Combines all sequences into a single DNAbin object
#'   \item Saves the combined sequences to the output file if \code{save = TRUE}
#'   \item Returns summary statistics
#' }
#' Note: This function does not remove duplicate sequences. All sequences from
#' all input files are included in the output.
#'
#' @examples
#' \dontrun{
#' # Basic usage: combine specific FASTA files
#' result <- combineFASTA(input_files = c("file1.fasta", "file2.fasta", "file3.fasta"))
#'
#' # Combine files with custom output name
#' result <- combineFASTA(
#'   input_files = c("data/gene1.fasta", "data/gene2.fasta"),
#'   output_file = "all_genes.fasta"
#' )
#'
#' # Return results without saving to disk
#' result <- combineFASTA(
#'   input_files = c("temp1.fasta", "temp2.fasta"),
#'   save = FALSE,
#'   verbose = TRUE
#' )
#'
#' # Access results
#' sequences <- result$sequences
#' summary_stats <- result$summary
#' print(summary_stats)
#' }
#'
#' @importFrom ape read.FASTA write.FASTA
#' @importFrom utils head tail
#'
#' @export
#'
combineFASTA <- function(input_files = NULL,
                         output_file = "combined_sequences.fasta",
                         save = TRUE,
                         verbose = TRUE,
                         dir = "RESULTS_combineFASTA") {

  # Validate inputs
  if (is.null(input_files) || length(input_files) == 0) {
    stop("input_files parameter is required. Please specify which FASTA files to combine.")
  }

  if (!is.character(input_files)) {
    stop("input_files must be a character vector of file paths.")
  }

  # Check that all files exist
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0) {
    stop("The following files do not exist:\n  ",
         paste(missing_files, collapse = "\n  "))
  }

  # Create output directory if saving
  if (save) {
    foldername <- paste0(dir, "/", format(Sys.time(), "%d%b%Y"))
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    if (!dir.exists(foldername)) {
      dir.create(foldername, recursive = TRUE)
    }
    output_path <- file.path(foldername, output_file)
  } else {
    output_path <- NULL
    if (verbose) {
      message("Running in memory-only mode (save = FALSE). Results will not be saved to disk.")
    }
  }

  if (verbose) {
    message("Processing ", length(input_files), " FASTA file(s):")
    for (file in input_files) {
      message("  - ", basename(file))
    }
    if (save) {
      message("Output will be saved to: ", output_path)
    }
  }

  # Initialize storage - create empty DNAbin object
  all_sequences <- list()
  seq_names <- character()
  file_info <- data.frame(
    filename = character(),
    sequences = numeric(),
    stringsAsFactors = FALSE
  )

  # Read and combine sequences
  total_sequences <- 0

  for (i in seq_along(input_files)) {
    file <- input_files[i]

    if (verbose) {
      message("\nProcessing: ", basename(file), " (", i, "/", length(input_files), ")")
    }

    tryCatch({
      # Read sequences
      sequences <- ape::read.FASTA(file)

      if (is.null(sequences) || length(sequences) == 0) {
        if (verbose) {
          message("  No sequences found in this file.")
        }
        next
      }

      # Store sequences - ensure we're working with DNAbin objects

        # If it's a single DNAbin sequence
        all_sequences <- c(all_sequences, list(sequences))
        seq_names <- c(seq_names, names(sequences))

      # Calculate sequence lengths for this file
      seq_lengths <- sapply(sequences, length)

      # Record file info
      file_info <- rbind(file_info, data.frame(
        filename = basename(file),
        sequences = length(sequences),
        stringsAsFactors = FALSE
      ))

      total_sequences <- total_sequences + length(sequences)

      if (verbose) {
        message("  Added ", length(sequences), " sequence(s) from this file.")
      }

    }, error = function(e) {
      warning("Error reading file '", basename(file), "': ", e$message)
    })
  }

  if (length(all_sequences) == 0) {
    stop("No sequences could be read from the specified files.")
  }

  # Convert to proper DNAbin object
  # First, ensure we have a proper DNAbin structure
  if (verbose) {
    message("\nSuccessfully read ", total_sequences, " sequences from ",
            nrow(file_info), " file(s).")
    message("Note: Duplicate sequences are NOT removed. All sequences are included.")
  }

  # Name the sequences properly
  if (length(seq_names) == length(all_sequences)) {
    names(all_sequences) <- seq_names
  }

  # Save to file if requested
  if (save) {
    if (verbose) {
      message("\nWriting combined sequences to: ", output_path)
    }

    # If all_sequences is a list of DNAbin objects
    merged_sequences <- do.call(c, all_sequences)

    # Try different approaches to write FASTA
      # Use write.FASTA if it's a proper DNAbin
        ape::write.FASTA(merged_sequences, output_path)

  }

  # Create summary statistics
  final_lengths <- sapply(all_sequences, length)
  summary_stats <- data.frame(
    total_files = nrow(file_info),
    total_sequences = length(all_sequences),
    output_file = ifelse(save, output_file, "Not saved to disk"),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("\n", .strrep("-", 50))
    message("COMBINATION COMPLETE")
    message(.strrep("-", 50))
    message("Summary:")
    message("  Input files: ", summary_stats$total_files)
    message("  Total sequences: ", summary_stats$total_sequences)
    if (save) {
      message("  Output file: ", output_path)
    } else {
      message("  Output: Returned as R object (not saved to disk)")
    }
    message(.strrep("-", 50))
  }

  # Ensure proper DNAbin class
  if (inherits(all_sequences, "list") && length(all_sequences) > 0) {
    if (inherits(all_sequences[[1]], "DNAbin")) {
      class(all_sequences) <- c("DNAbin", "list")
    }
  }

  # Return results
  result <- list(
    sequences = all_sequences,
    summary = summary_stats,
    file_info = file_info
  )

  if (save) {
    result$output_path <- output_path
  }

  class(result) <- "combinedFASTA"

  return(invisible(result))
}


# Helper function for repeated strings
.strrep <- function(x, times) {
  paste(rep(x, times), collapse = "")
}

