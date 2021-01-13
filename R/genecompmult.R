# Compares two individual gene datasets (DNA alignments) and make them with the same
# number of taxa. This function also needs the shortaxlabelsmult function that makes
# the tip labels with just the scientific name

# Used inside catmultGenes function

# Author: Domingos Cardoso and Quezia Cavalcante

.genecompmult <- function(...,
                         data = NULL,
                         loop = NULL,
                         shortaxlabel = TRUE,
                         missdata = TRUE,
                         outgroup = NULL) {

  twogenes <- list(...)
  names(twogenes)[1] <- names(data[1])
  gene_1 <- twogenes[[1]]
  names(twogenes)[2] <- names(data[loop+1])
  gene_2 <- twogenes[[2]]


  if (class(gene_1) == "list") {
    for (i in 1:length(gene_1)) {
      gene_1[[i]] <- paste(gene_1[[i]], collapse = "")
      gene_1[[i]] <- toupper(gene_1[[i]])
    }
    spp <- stats::setNames(data.frame(names(gene_1)), "species")
    seqs <- stats::setNames(data.frame(unlist(gene_1, use.names = FALSE)), "sequence")
    gene_1 <- cbind(spp, seqs)
  }

  if (class(gene_2) == "list") {

    for (i in 1:length(gene_2)) {
      gene_2[[i]] <- paste(gene_2[[i]], collapse = "")
      gene_2[[i]] <- toupper(gene_2[[i]])
    }
    spp <- stats::setNames(data.frame(names(gene_2)), "species")
    seqs <- stats::setNames(data.frame(unlist(gene_2, use.names = FALSE)), "sequence")
    gene_2 <- cbind(spp, seqs)
  }

  colnames(gene_1) <- c("species", "sequence")
  colnames(gene_2) <- c("species", "sequence")

  if (missdata) {
    cat(cat("Gene comparison will include missing data into",
            names(twogenes[1])), "", sep = "\n")
    cat("PS. The resulting dataset will have sequences sorted alphabetically by taxon.", "", sep = "\n")
    taxlabs <- gsub("(_[^_]+_[^_]+)_.*", "\\1", gene_1$species)
    gene_2 <- .shortaxlabelsmult(gene_2)
    unmatching_seqs <- gene_2[which(!gene_2$species %in% taxlabs), ]
    if (length(unmatching_seqs$species) != 0) {
      gene_1[nrow(gene_1) + length(unmatching_seqs$species), ] <- NA
    }
    gene_1$species <- ifelse(is.na(gene_1$species),
                            as.character(unmatching_seqs$species),
                            as.character(gene_1$species))
    gene_1$sequence <- ifelse(is.na(gene_1$sequence),
                             as.character(paste(rep("?", times = nchar(as.character(gene_1[1,2]))), collapse = "")),
                             as.character(gene_1$sequence))
    gene_1 <- dplyr::arrange(gene_1, species)
  }

  if (missdata == FALSE) {
    # cat(cat("Gene comparison will exclude sequence from", names(twogenes[1]),
    #         "that is not in", names(twogenes[2])), "", sep = "\n")
    if (is.null(outgroup)) {
      cat(cat("You have not provided any outgroup"), "", sep = "\n")
      cat(cat("Fully matching the genes..."), "", sep = "\n")
      taxlabs <- gsub("(_[^_]+_[^_]+)_.*", "\\1", gene_1$species)
      gene_1 <- .shortaxlabelsmult(gene_1)
      gene_2 <- .shortaxlabelsmult(gene_2)
      matching_seqs_2 <- gene_2[which(gene_2$species %in% taxlabs), ]
      gene_2 <- stats::na.omit(matching_seqs_2)
      matching_labs <- gene_2$species
      matching_seqs_1 <- gene_1[which(gene_1$species %in% matching_labs), ]

      gene_1 <- stats::na.omit(matching_seqs_1)

      gene_1 <- dplyr::arrange(gene_1, species)
    }

    if (is.null(outgroup) == FALSE) {
      cat(cat("You have provided outgroup:",
              paste0(outgroup, collapse = ", ")), "", sep = "\n")

      # Finding outgroup in each dataset
      length_1 <- length(which(grepl(paste(outgroup, collapse = "|"),
                                     gene_1$species) == TRUE))
      length_2 <- length(which(grepl(paste(outgroup, collapse = "|"),
                                     gene_2$species) == TRUE))

      if (length_1 == 0 & length_2 == 0) {
        stop("OUTGROUP must match exactly with any taxon name in the DNA alignments.\n",
             "Make sure outgroups are as VECTOR or LIST provided they are more than one.\n",
             "Find help also at DBOSLab-UFBA (Domingos Cardoso; cardosobot@gmail.com)")
      }

      cat(cat("Outgroup is present in", names(twogenes[1]), "or", names(twogenes[2])), "", sep = "\n")

      if (length_1 >= 1) {
        # Saving outgroup from gene_1
        seqs_1_out <- gene_1[which(grepl(paste(outgroup, collapse = "|"),
                                         gene_1$species) == TRUE),]
        seqs_1_out <- dplyr::arrange(seqs_1_out, species)
        # Saving ingroup from gene_1
        seqs_1_in <- gene_1[which(grepl(paste(outgroup, collapse = "|"),
                                        gene_1$species) == FALSE),]
        seqs_1_in <- dplyr::arrange(seqs_1_in, species)
      }

      if (length_1 < 1) {
        # Saving ingroup from gene_1
        seqs_1_in <- gene_1
      }

      if (length_2 >= 1) {
        # Saving outgroup from gene_2
        seqs_2_out <- gene_2[which(grepl(paste(outgroup, collapse = "|"),
                                         gene_2$species) == TRUE), ]
        seqs_2_out <- dplyr::arrange(seqs_2_out, species)
        # Saving ingroup from gene_2
        seqs_2_in <- gene_2[which(grepl(paste(outgroup, collapse = "|"),
                                        gene_2$species) == FALSE), ]
        seqs_2_in <- dplyr::arrange(seqs_2_in, species)
      }

      if (length_2 < 1) {
        # Saving ingroup from gene_2
        seqs_2_in <- gene_2
      }

      # Fully matching ingroup sequences
      taxlabs <- gsub("(_[^_]+_[^_]+)_.*", "\\1", seqs_1_in$species)
      seqs_2_in <- .shortaxlabelsmult(seqs_2_in)
      seqs_1_in <- .shortaxlabelsmult(seqs_1_in)
      matching_seqs_2_in <- seqs_2_in[which(seqs_2_in$species %in% taxlabs), ]
      seqs_2_in <- stats::na.omit(matching_seqs_2_in)
      matching_labs <- seqs_2_in$species
      matching_seqs_1_in <- seqs_1_in[which(seqs_1_in$species %in% matching_labs), ]
      seqs_1_in <- stats::na.omit(matching_seqs_1_in)
      seqs_1_in <- dplyr::arrange(seqs_1_in, species)

      # combining outgroup sequences
      if (length_2 < 1) {
        cat(cat("Outgroup sequences only in", names(twogenes[1])), "", sep = "\n")
        seqsout <- dplyr::arrange(seqs_1_out, species)
      }

      if (length_1 < 1) {
        cat(cat("Outgroup sequences only in", names(twogenes[2])), "", sep = "\n")
        seqs_2_out$sequence <- NA
        seqs_2_out$sequence <- ifelse(is.na(seqs_2_out$sequence),
                                    as.character(paste(rep("?", times = nchar(as.character(seqs_1_in[1,2]))),
                                                       collapse = "")),
                                    as.character(seqs_2_out$sequence))
        seqsout <- dplyr::arrange(seqs_2_out, species)
      }

      if (length_1 >= 1 & length_2 >= 1) {
        cat(cat("Outgroup sequences in both genes", ".", sep = ""), "", sep = "\n")
        taxlabsout <- gsub("(_[^_]+)_.*", "\\1", seqs_1_out$species)
        seqs_2_out <- .shortaxlabels(seqs_2_out)
        unmatching_seqsout <- seqs_2_out[which(!seqs_2_out$species %in% taxlabsout), ]
        if (length(unmatching_seqsout$species) != 0) {
          seqs_1_out[nrow(seqs_1_out) + length(unmatching_seqsout$species),] <- NA
        }
        seqs_1_out$species <- ifelse(is.na(seqs_1_out$species),
                                   as.character(unmatching_seqsout$species),
                                   as.character(seqs_1_out$species))
        seqs_1_out$sequence <- ifelse(is.na(seqs_1_out$sequence),
                                    as.character(paste(rep("?", times = nchar(as.character(seqs_1_out[1,2]))),
                                                       collapse = "")),
                                    as.character(seqs_1_out$sequence))
        seqsout <- dplyr::arrange(seqs_1_out, species)
      }

      cat(cat("Fully matching sequences..."), "", sep = "\n")
      # combining outgroup with ingroup sequences
      gene_1 <- rbind(seqsout, seqs_1_in)

      }
    }

  if (shortaxlabel) {
    cat("Shortening taxon labels...", "", sep = "\n")
    gene_1 <- .shortaxlabelsmult(gene_1)
  }

  cat(cat("Match between", names(twogenes[1]),
          "and", names(twogenes[2]), "is finished!"), "",
      sep = "\n")

  return(gene_1)
}
