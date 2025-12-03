#' Automated Evolutionary Model Selection and MrBayes Command Generation
#'
#' @author Domingos Cardoso
#'
#' @description This function performs DNA substitution model selection using phangorn::modelTest,
#' identifies the best model according to AIC and BIC criteria, and generates
#' corresponding MrBayes commands for phylogenetic analysis.
#'
#' @param nexus_file_path Path to the NEXUS alignment file (required)
#' @param model_criteria Selection criterion: "AIC" or "BIC" (default: "BIC")
#' @param partition MrBayes partition number to apply model to (default: "all")
#' @param models_to_test Which models to test: "all" for all models,
#'                       "standard" for JC/F81/K80/HKY/SYM/GTR only,
#'                       or a character vector of specific models
#' @param include_G Should Gamma rate heterogeneity be tested? (default: TRUE)
#' @param include_I Should invariant sites be tested? (default: TRUE)
#' @param mc.cores Number of cores for parallel computation (default: 4)
#' @param verbose Print progress messages? (default: TRUE)
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
                         model_criteria = "BIC",
                         partition = "all",
                         models_to_test = "all",
                         include_G = TRUE,
                         include_I = TRUE,
                         mc.cores = 4,
                         verbose = TRUE,
                         dir = "RESULTS_evomodelTest") {

  # Validate inputs
  if (is.null(nexus_file_path)) {
    stop("Please provide the path to your NEXUS alignment file.")
  }

  if (!file.exists(nexus_file_path)) {
    stop(paste("File not found:", nexus_file_path))
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

  # Generate output filename
  base_name <- tools::file_path_sans_ext(basename(nexus_file_path))
  output_file <- file.path(foldername, paste0(base_name, "_model_selection.txt"))

  # Start writing to file
  sink(output_file)

  # Write comprehensive header with attribution
  cat("=", rep("=", 50), "\n", sep = "")
  cat("EVOLUTIONARY MODEL SELECTION REPORT\n")
  cat("=", rep("=", 50), "\n", sep = "")

  cat("ANALYSIS INFORMATION:\n")
  cat(rep("-", 50), "\n", sep = "")
  cat("Analysis performed on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("Input file: ", nexus_file_path, "\n", sep = "")
  cat("File size: ", file.size(nexus_file_path), " bytes\n", sep = "")
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
    message("[", format(Sys.time(), "%H:%M:%S"), "] Reading NEXUS alignment...\n", sep = "")
  }

  phy <- phangorn::read.phyDat(nexus_file_path, format = "nexus")

  cat("\nALIGNMENT SUMMARY:\n")
  cat(rep("-", 50), "\n", sep = "")
  cat("Number of taxa:", length(phy), "\n")
  cat("Number of sites:", sum(attr(phy, "weight")), "\n")
  cat("Data type:", attr(phy, "type"), "\n")

  # Calculate some basic statistics
  if (verbose) {
    message("[", format(Sys.time(), "%H:%M:%S"), "] Calculating alignment statistics...\n", sep = "")
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
    message("[", format(Sys.time(), "%H:%M:%S"), "] Running model selection (this may take a while)...\n", sep = "")
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

  # Calculate AIC and BIC weights FIRST, before using them in tables
  aic_weights <- exp(-0.5 * (model_test$AIC - min(model_test$AIC)))
  aic_weights <- aic_weights / sum(aic_weights)
  model_test$AICw <- aic_weights

  bic_weights <- exp(-0.5 * (model_test$BIC - min(model_test$BIC)))
  bic_weights <- bic_weights / sum(bic_weights)
  model_test$BICw <- bic_weights

  # Display results
  cat("\nMODEL SELECTION RESULTS:\n")
  cat(rep("-", 50), "\n", sep = "")

  # Top 5 models by each criterion - FIXED to keep everything in same row
  cat("\nTOP 5 MODELS BY AIC:\n")
  cat(rep("-", 70), "\n", sep = "")
  aic_sorted <- model_test[order(model_test$AIC), ]
  top5_aic <- aic_sorted[1:5, ]

  # Custom printing to keep everything in one row per model
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
  # Find second best by AIC
  second_aic <- model_test[order(model_test$AIC), ][2, ]
  if (!is.null(second_aic) && nrow(second_aic) > 0) {
    cat("2.", second_aic$Model, ":", round(second_aic$AICw, 4),
        "(evidence ratio:", round(best_aic$AICw / second_aic$AICw, 2), ")\n")
  }
  cat("\n\n")

  # Generate MrBayes commands
  if (verbose) {
    message("[", format(Sys.time(), "%H:%M:%S"), "] Generating MrBayes commands...\n\n", sep = "")
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

  mrbayes_commands <- model_to_mrbayes(
    model_name = selected_model,
    partition = partition
  )

  # Also generate for the other best models for comparison
  other_models <- list(
    "AIC" = best_aic$Model,
    "BIC" = best_bic$Model,
    "AICc" = best_aicc$Model
  )

  # Write MrBayes blocks
  cat("RECOMMENDED MRBAYES BLOCKS:\n")
  cat(rep("-", 50), "\n\n", sep = "")

  # Basic block for selected model
  cat("=", rep("=", 50), "\n", sep = "")
  cat("MRBAYES BLOCK FOR", selected_model, "\n")
  cat("=", rep("=", 50), "\n", sep = "")

  cat("begin mrbayes;\n")
  cat("  [ Model settings from automated model selection ]\n")
  cat("  ", mrbayes_commands$lset, "\n")
  cat("  ", mrbayes_commands$prset, "\n")
  cat("  \n")
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

  # Close the output file
  sink()

  # Also print summary to console
  if (verbose) {
    message("\n[", format(Sys.time(), "%H:%M:%S"), "] Analysis complete!\n", sep = "")
    message("Results saved to: ", output_file, "\n")
    message("Raw results saved to: ", raw_results_file, "\n")
    message("\nSUMMARY:\n")
    message("Best model (by ", model_criteria, "): ", selected_model, "\n", sep = "")
    message("MrBayes commands generated for partition: ", partition, "\n", sep = "")
  }

  # Return invisible list with all results
  invisible(list(
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
      nexus_file = nexus_file_path,
      criteria = model_criteria,
      partition = partition,
      timestamp = timestamp
    )
  ))
}

#' Convert model name to MrBayes commands
#'
#' @param model_name Model name from phangorn::modelTest
#' @param partition MrBayes partition specification
#'
#' @return List with lset and prset commands
#'
#' @keywords internal
#'
model_to_mrbayes <- function(model_name,
                             partition = "all") {

  model <- as.character(model_name)

  # Extract components
  has_gamma <- grepl("\\+G", model)
  has_invariant <- grepl("\\+I", model)
  ncat <- ifelse(grepl("G\\(([0-9]+)\\)", model),
                 gsub(".*G\\(([0-9]+)\\).*", "\\1", model), "4")
  gamma_categories = ncat

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
  model_base <- gsub("\\+.*", "", model)  # Remove +G+I etc.

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

  # Find the model (with or without exact match)
  matched_model <- NULL
  for (pattern in names(revmatpr_list)) {
    if (grepl(paste0("^", pattern), model_base)) {
      matched_model <- revmatpr_list[[pattern]]
      break
    }
  }

  # Default to GTR if no match
  if (is.null(matched_model)) {
    warning(paste("Model", model_name, "not recognized. Defaulting to GTR."))
    matched_model <- revmatpr_list[["GTR"]]
  }

  # Generate commands
  lset_cmd <- paste0("lset applyto=(", partition, ") nst=",
                     matched_model$nst, " rates=", rates, ";")
  prset_cmd <- paste0("prset applyto=(", partition, ") revmatpr=",
                      matched_model$revmatpr, ";")

  return(list(
    lset = lset_cmd,
    prset = prset_cmd,
    model_info = list(
      base_model = model_base,
      has_gamma = has_gamma,
      has_invariant = has_invariant,
      gamma_categories = ncat
    )
  ))
}

