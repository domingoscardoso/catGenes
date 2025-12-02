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

  # Initialize storage
  all_sequences <- list()
  file_info <- data.frame(
    filename = character(),
    sequences = numeric(),
    min_length = numeric(),
    max_length = numeric(),
    mean_length = numeric(),
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

      if (length(sequences) == 0) {
        if (verbose) {
          message("  No sequences found in this file.")
        }
        next
      }

      # Store sequences
      all_sequences <- c(all_sequences, sequences)

      # Calculate sequence lengths for this file
      seq_lengths <- sapply(sequences, length)

      # Record file info
      file_info <- rbind(file_info, data.frame(
        filename = basename(file),
        sequences = length(sequences),
        min_length = min(seq_lengths),
        max_length = max(seq_lengths),
        mean_length = mean(seq_lengths),
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

  if (verbose) {
    message("\nSuccessfully read ", total_sequences, " sequences from ",
            nrow(file_info), " file(s).")
    message("Note: Duplicate sequences are NOT removed. All sequences are included.")
  }

  # Save to file if requested
  if (save) {
    if (verbose) {
      message("\nWriting combined sequences to: ", output_path)
    }

    ape::write.FASTA(all_sequences, output_path)

    if (verbose) {
      message("File successfully saved.")
    }
  }

  # Create summary statistics
  final_lengths <- sapply(all_sequences, length)
  summary_stats <- data.frame(
    total_files = nrow(file_info),
    total_sequences = length(all_sequences),
    min_sequence_length = min(final_lengths),
    max_sequence_length = max(final_lengths),
    mean_sequence_length = mean(final_lengths),
    median_sequence_length = median(final_lengths),
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
    message("  Sequence length range: ",
            summary_stats$min_sequence_length, " - ",
            summary_stats$max_sequence_length)
    message("  Mean sequence length: ",
            round(summary_stats$mean_sequence_length, 1))
    if (save) {
      message("  Output file: ", output_path)
    } else {
      message("  Output: Returned as R object (not saved to disk)")
    }
    message(.strrep("-", 50))
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


# Print method for combineFASTA results
print.combinedFASTA <- function(x, ...) {
  cat("FASTA Combination Results\n")
  cat(.strrep("=", 40), "\n")
  if ("output_path" %in% names(x)) {
    cat("Output file:", x$output_path, "\n")
  } else {
    cat("Output: Not saved to disk\n")
  }
  cat("Total sequences:", x$summary$total_sequences, "\n")
  cat("Sequence length range:",
      x$summary$min_sequence_length, "-",
      x$summary$max_sequence_length, "\n")
  cat("Mean length:", round(x$summary$mean_sequence_length, 1), "\n")

  if (!is.null(x$file_info) && nrow(x$file_info) > 0) {
    cat("\nInput files:\n")
    print(x$file_info[, c("filename", "sequences", "mean_length")])
  }

  invisible(x)
}


# Helper function for repeated strings (compatible with older R versions)
.strrep <- function(x, times) {
  paste(rep(x, times), collapse = "")
}
