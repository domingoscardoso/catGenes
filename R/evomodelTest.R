#' Automated Evolutionary Model Selection and MrBayes Command Generation
#'
#' @author Domingos Cardoso
#'
#' @description This function performs DNA substitution model selection using
#' the modelTest function from [phangorn](https://klausvigo.github.io/phangorn/index.html)
#' package, identifies the best model according to AIC and BIC criteria, and generates
#' corresponding MrBayes commands for phylogenetic analysis.
#'
#' @usage
#' evomodelTest(nexus_file_path = NULL,
#'              model_criteria = "AIC",
#'              models_to_test = "all",
#'              include_G = TRUE,
#'              include_I = TRUE,
#'              mc.cores = 4,
#'              verbose = TRUE,
#'              dir = "RESULTS_evomodelTest")
#'
#' @param nexus_file_path Path to the NEXUS alignment file(s).
#' @param model_criteria Selection criterion: "AIC" or "BIC" (default: "BIC").
#' @param models_to_test Which models to test: "all" for all models, "standard"
#' for JC/F81/K80/HKY/SYM/GTR only, or a character vector of specific models.
#' @param include_G Should Gamma rate heterogeneity be tested? (default: TRUE).
#' @param include_I Should invariant sites be tested? (default: TRUE).
#' @param mc.cores Number of cores for parallel computation (default: 4).
#' @param verbose Print progress messages? (default: TRUE).
#' @param dir The path to the directory where the model selection results will
#' be saved. The default is to create a directory named **RESULTS_evomodelTest**
#' and the resulting analysis will be saved within a subfolder named after the
#' current date (format: "ddMonYYYY", e.g., "02Dec2025").
#'
#' @return Invisibly returns a list containing all results. Saves a detailed
#'         report to a text file in the results directory.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' evomodelTest("alignment.nex")
#'
#' # Specify AIC criterion and partition 1
#' evomodelTest("alignment.nex", model_criteria = "AIC", partition = "1")
#'
#' # Test only standard models
#' evomodelTest("alignment.nex", models_to_test = "standard")
#' }
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom phangorn read.phyDat modelTest pml.control
#' @importFrom ape base.freq as.DNAbin
#'
#' @export
#'
evomodelTest <- function(nexus_file_path = NULL,
                         model_criteria = "AIC",
                         models_to_test = "all",
                         include_G = TRUE,
                         include_I = TRUE,
                         mc.cores = 4,
                         verbose = TRUE,
                         dir = "RESULTS_evomodelTest") {

  nexus_files <- sort(list.files(nexus_file_path, full.names = TRUE))

  # Validate inputs
  if (is.null(nexus_file_path)) {
    stop("Please provide the path to your NEXUS alignment file(s).")
  }

  if (length(nexus_files) == 0) {
    stop(paste("No files found within the folder:", nexus_file_path))
  }

  if (!model_criteria %in% c("AIC", "BIC", "AICc")) {
    warning("model_criteria should be 'AIC', 'BIC', or 'AICc'. Using 'AIC'.")
    model_criteria <- "AIC"
  }

  # Create results directory structure
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  foldername <- file.path(dir, format(Sys.time(), "%d%b%Y"))
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  if (!dir.exists(foldername)) {
    dir.create(foldername, recursive = TRUE)
  }

  # Initialize storage
  selected_models <- vector("character", length(nexus_files))
  all_results <- list()
  file_results <- list()

  # ============================================================
  # PART 1: PROCESS EACH FILE INDIVIDUALLY
  # ============================================================
  for (l in seq_along(nexus_files)) {

    if (verbose) {
      message("[", format(Sys.time(), "%H:%M:%S"), "] Processing file ", l,
              " of ", length(nexus_files), ": ", basename(nexus_files[l]), sep = "")
    }

    # Generate output filename
    base_name <- tools::file_path_sans_ext(basename(nexus_files[l]))
    output_file <- file.path(foldername, paste0(base_name, "_model_selection.txt"))

    # Start writing to file
    sink(output_file)

    # Write comprehensive header with attribution
    cat("=", rep("=", 50), "\n", sep = "")
    cat("EVOLUTIONARY MODEL SELECTION REPORT\n")
    cat("=", rep("=", 50), "\n\n", sep = "")

    cat("ANALYSIS INFORMATION:\n")
    cat(rep("-", 50), "\n", sep = "")
    cat("Analysis performed on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("Input file: ", basename(nexus_files[l]), "\n", sep = "")
    cat("File size: ", file.size(nexus_files[l]), " bytes\n", sep = "")
    cat("Results directory: ", foldername, "\n\n", sep = "")

    cat("SOFTWARE CREDITS:\n")
    cat(rep("-", 50), "\n", sep = "")

    cat("Model selection performed with phangorn's function modelTest\n", sep = "")
    cat("R package: phangorn (Maximum Likelihood phylogenetic analysis in R)\n")
    cat("Authors: Klaus Schliep, Emmanuel Paradis, et al.\n\n")

    cat("This analysis was conducted using the evomodelTest R function,\n")
    cat("developed for integration with the catGenes pipeline by\n")
    cat("Domingos Cardoso (Jardim Botanico do Rio de Janeiro - JBRJ)\n")
    cat("GitHub repository: https://github.com/domingoscardoso/catGenes\n\n")

    cat(rep("-", 50), "\n\n", sep = "")

    # Read and analyze the alignment
    if (verbose) {
      message("[", format(Sys.time(), "%H:%M:%S"), "] Reading NEXUS alignment...", sep = "")
    }

    phy <- phangorn::read.phyDat(nexus_files[l], format = "nexus")

    cat("\nALIGNMENT SUMMARY:\n")
    cat(rep("-", 50), "\n", sep = "")
    cat("Number of taxa:", length(phy), "\n")
    cat("Number of sites:", sum(attr(phy, "weight")), "\n")
    cat("Data type:", attr(phy, "type"), "\n")

    # Calculate some basic statistics
    if (verbose) {
      message("[", format(Sys.time(), "%H:%M:%S"), "] Calculating alignment statistics...", sep = "")
    }

    # Convert to DNAbin for basic stats
    dna <- ape::as.DNAbin(phy)
    cat("Proportion of gaps/missing data:",
        round(sum(ape::base.freq(dna, all = TRUE)["n"] +
                    ape::base.freq(dna, all = TRUE)["-"]) /
                length(dna) * 100, 2), "%\n")

    base_freq <- ape::base.freq(dna)
    cat("Base frequencies:\n")
    cat("  A:", round(base_freq["a"], 4),
        "  C:", round(base_freq["c"], 4),
        "  G:", round(base_freq["g"], 4),
        "  T:", round(base_freq["t"], 4), "\n\n")

    # Run model selection
    if (verbose) {
      message("[", format(Sys.time(), "%H:%M:%S"), "] Running model selection (this may take a while)...", sep = "")
    }

    model_test <- phangorn::modelTest(
      phy,
      model = models_to_test,
      G = include_G,
      I = include_I,
      control = phangorn::pml.control(trace = 0),
      multicore = (mc.cores > 1),
      mc.cores = mc.cores
    )

    # Save raw results
    raw_results_file <- file.path(foldername, paste0(base_name, "_model_test_raw.csv"))
    write.csv(model_test, raw_results_file, row.names = FALSE)

    # Identify best models
    best_aic <- model_test[which.min(model_test$AIC), ]
    best_bic <- model_test[which.min(model_test$BIC), ]
    best_aicc <- model_test[which.min(model_test$AICc), ]

    # Calculate AIC and BIC weights
    aic_weights <- exp(-0.5 * (model_test$AIC - min(model_test$AIC)))
    aic_weights <- aic_weights / sum(aic_weights)
    model_test$AICw <- aic_weights

    bic_weights <- exp(-0.5 * (model_test$BIC - min(model_test$BIC)))
    bic_weights <- bic_weights / sum(bic_weights)
    model_test$BICw <- bic_weights

    # Display results
    cat("\nMODEL SELECTION RESULTS:\n")
    cat(rep("-", 50), "\n", sep = "")

    # Top 5 models by each criterion
    cat("\nTOP 5 MODELS BY AIC:\n")
    cat(rep("-", 70), "\n", sep = "")
    aic_sorted <- model_test[order(model_test$AIC), ]
    top5_aic <- aic_sorted[1:5, ]

    cat(sprintf("%-20s %10s %12s %10s %10s\n",
                "Model", "AIC", "AICw", "BIC", "logLik"))
    cat(rep("-", 70), "\n", sep = "")
    for(i in 1:nrow(top5_aic)) {
      cat(sprintf("%-20s %10.2f %12.6f %10.2f %10.2f\n",
                  top5_aic$Model[i],
                  top5_aic$AIC[i],
                  top5_aic$AICw[i],
                  top5_aic$BIC[i],
                  top5_aic$logLik[i]))
    }
    cat("\n")

    cat("TOP 5 MODELS BY BIC:\n")
    cat(rep("-", 70), "\n", sep = "")
    bic_sorted <- model_test[order(model_test$BIC), ]
    top5_bic <- bic_sorted[1:5, ]

    cat(sprintf("%-20s %10s %12s %10s %10s\n",
                "Model", "BIC", "BICw", "AIC", "logLik"))
    cat(rep("-", 70), "\n", sep = "")
    for(i in 1:nrow(top5_bic)) {
      cat(sprintf("%-20s %10.2f %12.6f %10.2f %10.2f\n",
                  top5_bic$Model[i],
                  top5_bic$BIC[i],
                  top5_bic$BICw[i],
                  top5_bic$AIC[i],
                  top5_bic$logLik[i]))
    }
    cat("\n")

    cat("BEST MODEL BY EACH CRITERION:\n")
    cat(rep("-", 50), "\n", sep = "")
    cat("AIC  :", best_aic$Model, "(ΔAIC = 0.000)\n")
    cat("AICc :", best_aicc$Model, "(ΔAICc = 0.000)\n")
    cat("BIC  :", best_bic$Model, "(ΔBIC = 0.000)\n\n")

    cat("AIC WEIGHTS (relative support):\n")
    cat(rep("-", 50), "\n", sep = "")
    cat("1.", best_aic$Model, ":", round(best_aic$AICw, 4), "\n")
    second_aic <- model_test[order(model_test$AIC), ][2, ]
    if (!is.null(second_aic) && nrow(second_aic) > 0) {
      cat("2.", second_aic$Model, ":", round(second_aic$AICw, 4),
          "(evidence ratio:", round(best_aic$AICw / second_aic$AICw, 2), ")\n")
    }
    cat("\n\n")

    # ============================================================
    # Generate MrBayes commands for INDIVIDUAL file
    # ============================================================
    if (verbose) {
      cat("[", format(Sys.time(), "%H:%M:%S"), "] Generating MrBayes commands...\n\n", sep = "")
    }

    cat("MRBAYES COMMANDS FOR SELECTED MODEL:\n")
    cat(rep("-", 50), "\n", sep = "")

    # Generate commands for best model by selected criterion
    if (model_criteria == "AIC") {
      selected_model <- best_aic$Model
      cat("Using best AIC model:", selected_model, "\n\n")
    } else if (model_criteria == "BIC") {
      selected_model <- best_bic$Model
      cat("Using best BIC model:", selected_model, "\n\n")
    } else {
      selected_model <- best_aicc$Model
      cat("Using best AICc model:", selected_model, "\n\n")
    }
    selected_models[l] <- selected_model

    # Create a vector with this single model for evomodelCmds
    # Position in vector = partition number (will be 1 for single file)
    model_vector_single <- selected_model

    # Generate commands using evomodelCmds
    mrbayes_commands <- evomodelCmds(model_vector = model_vector_single)

    # ============================================================
    # Write MrBayes blocks based on evomodelCmds output
    # ============================================================
    cat("RECOMMENDED MRBAYES BLOCKS:\n")
    cat(rep("-", 50), "\n\n", sep = "")

    cat(rep("=", 50), "\n", sep = "")
    cat("MRBAYES BLOCK FOR THIS PARTITION\n")
    cat(rep("=", 50), "\n\n", sep = "")

    cat("begin mrbayes;\n")
    cat("  # ================================================\n")
    cat("  # MODEL SETTINGS (from evomodelTest)\n")
    cat("  # ================================================\n")
    cat("  # Selected model: ", selected_model, "\n", sep = "")
    cat("  # Selection criterion: ", model_criteria, "\n", sep = "")

    # Use the commands from evomodelCmds - they're already optimized
    # For single partition, evomodelCmds creates commands like "applyto=(1)"
    for (cmd in mrbayes_commands$all_commands) {
      cat("  ", cmd, "\n", sep = "")
    }

    cat("  \n")
    cat("  # ================================================\n")
    cat("  # ADDITIONAL MRBAYES SETTINGS\n")
    cat("  # ================================================\n")
    cat("  [ Base frequencies (empirical) ]\n")
    cat("  prset applyto=(all) statefreqpr=dirichlet(1,1,1,1);\n")
    cat("  \n")
    cat("  [ Allow different rates across partitions if multiple exist ]\n")
    cat("  prset ratepr=variable;\n")
    cat("  \n")
    cat("  [ Unlink parameters between partitions ]\n")
    cat("  unlink statefreqpr=(all) revmatpr=(all) shapepr=(all) pinvarpr=(all);\n")
    cat("  \n")
    cat("  [ MCMC settings (adjust as needed) ]\n")
    cat("  mcmc ngen=2000000 samplefreq=1000 printfreq=1000\n")
    cat("       nchains=4 nruns=2 savebrlens=yes;\n")
    cat("  \n")
    cat("  [ Convergence diagnostics ]\n")
    cat("  sumt burnin=500;\n")
    cat("  sump burnin=500;\n")
    cat("end;\n\n")

    # Store alternative models for comparison
    other_models <- list(
      "AIC" = best_aic$Model,
      "BIC" = best_bic$Model,
      "AICc" = best_aicc$Model
    )

    # Show alternative models if they differ
    if (any(sapply(other_models, function(x) x != selected_model))) {
      cat("\nALTERNATIVE MODELS (other best models):\n")
      cat(rep("-", 40), "\n", sep = "")
      for (criterion in names(other_models)) {
        if (other_models[[criterion]] != selected_model) {
          alt_commands <- evomodelCmds(model_vector = other_models[[criterion]])
          cat("# ", criterion, ": ", other_models[[criterion]], "\n", sep = "")
          cat("# lset:", alt_commands$lset_cmd[1], "\n")
          cat("# prset:", alt_commands$prset_cmd[1], "\n\n")
        }
      }
    }

    # Close the output file
    sink()

    # Store results for this file
    file_results[[l]] <- list(
      alignment = phy,
      model_test = model_test,
      best_models = list(
        AIC = best_aic,
        BIC = best_bic,
        AICc = best_aicc
      ),
      selected_model = selected_model,
      mrbayes_commands = mrbayes_commands,
      output_files = list(
        report = output_file,
        raw_data = raw_results_file
      ),
      settings = list(
        nexus_file = basename(nexus_files[l]),
        criteria = model_criteria,
        timestamp = timestamp
      )
    )
    names(file_results)[l] <- basename(nexus_files[l])

    # Print summary to console
    if (verbose) {
      message("[", format(Sys.time(), "%H:%M:%S"), "] Analysis complete for: ",
              basename(nexus_files[l]), sep = "")
      message("  Best model (by ", model_criteria, "): ", selected_model, sep = "")
      message("  Results saved to: ", output_file, "\n", sep = "")
    }
  } # End of for loop

  # ============================================================
  # PART 2: CREATE COMBINED MRBAYES COMMANDS FOR MULTIPLE FILES
  # ============================================================
  if (length(selected_models) > 1) {
    # Generate optimized commands for all partitions together
    combined_mrbayes_commands <- evomodelCmds(model_vector = selected_models)

    # ============================================================
    # CREATE SEPARATE TXT FILE WITH COMBINED MRBAYES BLOCKS
    # ============================================================
    compact_output_file <- file.path(foldername, "combined_mrbayes_blocks.txt")

    sink(compact_output_file)

    cat(rep("=", 50), "\n", sep = "")
    cat("COMBINED MRBAYES BLOCKS FOR PARTITIONED ANALYSIS\n")
    cat(rep("=", 50), "\n\n", sep = "")

    cat("ANALYSIS SUMMARY\n")
    cat(rep("-", 50), "\n", sep = "")
    cat("Total partitions: ", length(nexus_files), "\n", sep = "")
    cat("Model selection criterion: ", model_criteria, "\n", sep = "")
    cat("Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("Directory: ", foldername, "\n\n", sep = "")

    cat("PARTITION - MODEL MAPPING\n")
    cat(rep("-", 50), "\n", sep = "")
    for (i in seq_along(nexus_files)) {
      cat("Partition ", i, " (", basename(nexus_files[i]), "): ",
          selected_models[i], "\n", sep = "")
    }
    cat("\n")

    cat("OPTIMIZED MRBAYES COMMANDS\n")
    cat(rep("-", 50), "\n", sep = "")
    cat("# Commands are grouped by identical models to avoid duplication\n")
    cat("# Partitions with the same model share a single command\n\n")

    # Display the grouped commands
    cat("Model Configuration Groups:\n")
    for (i in seq_len(nrow(combined_mrbayes_commands$grouped_commands))) {
      group <- combined_mrbayes_commands$grouped_commands[i, ]
      cat("# Group ", i, ": Partitions ", group$partitions,
          " -> ", group$model, "\n", sep = "")
    }
    cat("\n")

    # The actual MrBayes commands
    for (cmd in combined_mrbayes_commands$all_commands) {
      cat(cmd, "\n")
    }

    cat("\n")
    cat("COMPLETE MRBAYES BLOCK FOR ALL PARTITIONS\n")
    cat(rep("-", 50), "\n\n", sep = "")

    cat("# Copy and paste this entire block into your NEXUS file\n\n")

    cat("begin mrbayes;\n")
    cat("  \n")
    cat("  # ================================================\n")
    cat("  # PARTITION DEFINITIONS (ADJUST TO YOUR DATA)\n")
    cat("  # ================================================\n")
    for (i in seq_along(nexus_files)) {
      cat("  # charset part", i, " = [YOUR_RANGE_HERE];  # ",
          basename(nexus_files[i]), "\n", sep = "")
    }
    cat("  # partition analysis = ", length(nexus_files), ": ",
        paste0("part", seq_along(nexus_files), collapse = ", "), ";\n", sep = "")
    cat("  # set partition = analysis;\n")
    cat("  \n")
    cat("  # ================================================\n")
    cat("  # MODEL SETTINGS (optimized by evomodelTest)\n")
    cat("  # ================================================\n")

    # Add the actual commands
    for (cmd in combined_mrbayes_commands$all_commands) {
      if (!grepl("^#", cmd)) {  # Skip comment lines from evomodelCmds
        cat("  ", cmd, "\n", sep = "")
      }
    }

    cat("  \n")
    cat("  # ================================================\n")
    cat("  # ADDITIONAL MRBAYES SETTINGS\n")
    cat("  # ================================================\n")
    cat("  [ Base frequencies (empirical) ]\n")
    cat("  prset applyto=(all) statefreqpr=dirichlet(1,1,1,1);\n")
    cat("  \n")
    cat("  [ Allow different rates across partitions ]\n")
    cat("  prset ratepr=variable;\n")
    cat("  \n")
    cat("  [ Unlink parameters between partitions ]\n")
    cat("  unlink statefreqpr=(all) revmatpr=(all) shapepr=(all) pinvarpr=(all);\n")
    cat("  \n")
    cat("  [ MCMC settings (adjust as needed) ]\n")
    cat("  mcmc ngen=2000000 samplefreq=1000 printfreq=1000\n")
    cat("       nchains=4 nruns=2 savebrlens=yes;\n")
    cat("  \n")
    cat("  [ Convergence diagnostics ]\n")
    cat("  sumt burnin=500;\n")
    cat("  sump burnin=500;\n")
    cat("end;\n")

    sink()

    if (verbose) {
      message("\n[", format(Sys.time(), "%H:%M:%S"),
              "] Created combined MrBayes blocks file!", sep = "")
      message("  Saved to: ", compact_output_file, sep = "")
      message("\n  OPTIMIZED COMMAND GROUPS:")
      for (i in seq_len(nrow(combined_mrbayes_commands$grouped_commands))) {
        group <- combined_mrbayes_commands$grouped_commands[i, ]
        message("    Group ", i, ": lset applyto=(", group$partitions,
                ") for ", group$model, sep = "")
      }
    }

    # Add combined results to the return object
    all_results <- list(
      individual_files = file_results,
      combined_analysis = list(
        selected_models = selected_models,
        mrbayes_commands = combined_mrbayes_commands,
        compact_output_file = compact_output_file,
        partition_mapping = setNames(selected_models,
                                     paste0("Partition ", seq_along(selected_models)))
      )
    )
  } else {
    # For single file, return simpler structure
    all_results <- file_results[[1]]
  }

  # Return results invisibly
  invisible(all_results)
}


#' Convert model names into efficient MrBayes commands
#'
#' @param model_vector A character vector of model names where the position/index
#'        corresponds to the MrBayes partition number (e.g., `c("GTR+G+I", "HKY+G", "GTR+G+I")`).
#'        Partition 1 gets the first model, Partition 2 the second, etc.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{grouped_commands}: A data frame showing partition grouping.
#'     \item \code{all_commands}: A character vector of unique `lset` and `prset` commands.
#'     \item \code{partition_map}: A named vector showing final partition assignments.
#'   }
#'
#' @examples
#' # Partitions 1 and 3 share the same model, partition 2 is different
#' models <- c("GTR+G(4)+I", "HKY+G(4)", "GTR+G(4)+I")
#' cmds <- evomodelCmds(model_vector = models)
#' cat(cmds$all_commands, sep = "\n")
#'
#' @export
#'
evomodelCmds <- function(model_vector) {

  # Input validation
  if (!is.character(model_vector) || length(model_vector) == 0) {
    stop("'model_vector' must be a non-empty character vector.")
  }

  # Helper function: Parse a single model (same as before)
  .parse_single_model <- function(model_name) {
    model <- as.character(model_name)

    # Extract components
    has_gamma <- grepl("\\+G", model)
    has_invariant <- grepl("\\+I", model)
    ncat <- ifelse(grepl("G\\(([0-9]+)\\)", model),
                   gsub(".*G\\(([0-9]+)\\).*", "\\1", model), "4")

    # Determine rates parameter
    if (has_gamma & has_invariant) {
      rates <- paste0("invgamma ngammacat=", ncat)
    } else if (has_gamma) {
      rates <- paste0("gamma ngammacat=", ncat)
    } else if (has_invariant) {
      rates <- "propinv"
    } else {
      rates <- "equal"
    }

    # Determine nst and revmatpr based on model
    model_base <- gsub("\\+.*", "", model)

    revmatpr_list <- list(
      JC = list(nst = 1, revmatpr = "fixed(1,1,1,1,1,1)"),
      F81 = list(nst = 1, revmatpr = "fixed(1,1,1,1,1,1)"),
      K80 = list(nst = 2, revmatpr = "fixed(1,2,1,1,2,1)"),
      HKY = list(nst = 2, revmatpr = "fixed(1,2,1,1,2,1)"),
      TrN = list(nst = 6, revmatpr = "fixed(1,2,1,1,3,1)"),
      TrNe = list(nst = 6, revmatpr = "fixed(1,2,1,1,3,1)"),
      TIM1 = list(nst = 6, revmatpr = "fixed(1,2,3,1,2,4)"),
      TIM1e = list(nst = 6, revmatpr = "fixed(1,2,3,1,2,4)"),
      TIM2 = list(nst = 6, revmatpr = "fixed(1,2,1,1,3,1)"),
      TIM2e = list(nst = 6, revmatpr = "fixed(1,2,1,1,3,1)"),
      TIM3 = list(nst = 6, revmatpr = "fixed(1,2,1,3,2,4)"),
      TIM3e = list(nst = 6, revmatpr = "fixed(1,2,1,3,2,4)"),
      TPM1 = list(nst = 6, revmatpr = "fixed(1,2,3,2,1,3)"),
      TPM1u = list(nst = 6, revmatpr = "fixed(1,2,3,2,1,3)"),
      TPM2 = list(nst = 6, revmatpr = "fixed(1,2,3,2,1,3)"),
      TPM2u = list(nst = 6, revmatpr = "fixed(1,2,3,2,1,3)"),
      TPM3 = list(nst = 6, revmatpr = "fixed(1,2,1,3,2,3)"),
      TPM3u = list(nst = 6, revmatpr = "fixed(1,2,1,3,2,3)"),
      K81 = list(nst = 6, revmatpr = "fixed(1,2,1,1,2,1)"),
      TVM = list(nst = 6, revmatpr = "fixed(1,2,3,4,5,2)"),
      TVMe = list(nst = 6, revmatpr = "fixed(1,2,3,4,5,2)"),
      SYM = list(nst = 6, revmatpr = "fixed(1,2,3,4,5,6)"),
      GTR = list(nst = 6, revmatpr = "fixed(1,2,3,4,5,6)")
    )

    matched_model <- NULL
    for (pattern in names(revmatpr_list)) {
      if (grepl(paste0("^", pattern), model_base)) {
        matched_model <- revmatpr_list[[pattern]]
        break
      }
    }

    if (is.null(matched_model)) {
      warning(paste("Model", model_name, "not recognized. Defaulting to GTR."))
      matched_model <- revmatpr_list[["GTR"]]
    }

    # Return a unique model signature for grouping
    return(list(
      nst = matched_model$nst,
      rates = rates,
      revmatpr = matched_model$revmatpr,
      model_name = model,
      signature = paste(matched_model$nst, rates, matched_model$revmatpr, sep = "|")
    ))
  }

  # Step 1: Parse all models
  parsed_list <- lapply(model_vector, .parse_single_model)

  # Create data frame with partition numbers as positions
  df <- data.frame(
    partition = seq_along(model_vector),  # Position = partition number
    model = model_vector,
    nst = sapply(parsed_list, function(x) x$nst),
    rates = sapply(parsed_list, function(x) x$rates),
    revmatpr = sapply(parsed_list, function(x) x$revmatpr),
    signature = sapply(parsed_list, function(x) x$signature),
    stringsAsFactors = FALSE
  )

  # Step 2: Group partitions by identical model signature
  # Create grouping by identical model configurations
  df$group_id <- as.numeric(factor(df$signature, levels = unique(df$signature)))

  # Create partition groupings like "1,3" for partitions 1 and 3 with same model
  partition_groups <- tapply(df$partition, df$group_id, function(x) {
    if (length(x) == 1) {
      return(as.character(x))
    } else {
      return(paste(x, collapse = ","))
    }
  })

  # Get unique model configurations per group
  unique_configs <- df[!duplicated(df$group_id),
                       c("group_id", "model", "nst", "rates", "revmatpr", "signature")]

  # Combine into final grouped data frame
  grouped_df <- data.frame(
    partitions = unname(partition_groups),
    unique_configs[, c("model", "nst", "rates", "revmatpr", "signature")],
    stringsAsFactors = FALSE
  )

  # Step 3: Generate MrBayes commands
  grouped_df$lset_cmd <- paste0("lset applyto=(", grouped_df$partitions, ") nst=",
                                grouped_df$nst, " rates=", grouped_df$rates, ";")
  grouped_df$prset_cmd <- paste0("prset applyto=(", grouped_df$partitions, ") revmatpr=",
                                 grouped_df$revmatpr, ";")

  # Create clean command vector
  all_cmds <- c(
    "# Model settings from partitioned analysis",
    paste0("# Generated for ", length(model_vector), " partition(s)"),
    grouped_df$lset_cmd,
    grouped_df$prset_cmd
  )

  # Create a partition map for reference
  partition_map <- setNames(model_vector, seq_along(model_vector))

  return(list(
    grouped_commands = grouped_df,
    all_commands = all_cmds,
    partition_map = partition_map,
    raw_data = df
  ))
}
