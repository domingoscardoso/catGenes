# Auxiliary functions to support main functions
# Author: Domingos Cardoso


#-------------------------------------------------------------------------------
# Auxiliary function to create a named list of the input individual genes
# Used inside the functions catfullGenes, catmultGenes, writeNexus, and writePhylip

.namedlist <- function(...) {
  nms <- sapply(as.list(substitute(list(...))), deparse)[-1]
  setNames(list(...), nms)
}


#-------------------------------------------------------------------------------
# Auxiliary function to delete all identifiers of the DNA sequence (e.g. collector
# and GenBank numbers) so as to maintain just the scientific names
# Used inside the function catfullGenes

.shortaxlabels <- function(x) {

  colnames(x) <- c("species", "sequence")
  x$species <- gsub("(_[^_]+)_.*", "\\1", x$species)

  return(x)
}


#-------------------------------------------------------------------------------
# Auxiliary function to delete all extra identifiers of the DNA sequence (e.g. GenBank numbers)
# so as to maintain just the scientific names and associated collector number
# Used inside the function catmultGenes

.shortaxlabelsmult <- function(x) {

  colnames(x) <- c("species", "sequence")
  x$species <- gsub("(_[^_]+_[^_]+)_.*", "\\1", x$species)

  return(x)
}


#-------------------------------------------------------------------------------
# Auxiliary function to obtain the intersection of a list of vectors
# Used inside the function delGaps

.intersectAll <- function(...) {
  args <- list(...)
  nargs <- length(args[[1]])
  if (nargs <= 1) {
    if (nargs == 1 && is.list(unlist(args[[1]][[1]]))) {
      do.call("intersectAll", unlist(args[[1]][[1]]))
    } else {
      stop("Cannot evaluate intersection fewer than 2 arguments")
    }
  } else if (nargs == 2) {
    intersect(unlist(args[[1]][[1]]), unlist(args[[1]][[2]]))
  } else {
    intersect(unlist(args[[1]][[1]]), .intersectAll(args[[1]][-1]))
  }
}


#-------------------------------------------------------------------------------
# Auxiliary function for returning TRUE if all values are equal and FALSE if it
# contains different values.
# Used inside the function dropSeq

equalnumb <- function(x) {
  res <- FALSE
  x <- na.omit(as.vector(x))
  if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
  return(res)

}


#-------------------------------------------------------------------------------
# Auxiliary function to adjust names with cf., aff, or infraspecific taxa
# Used inside catmultGenes, catfullGenes and dropSeq functions

.adjustnames <- function(datset,
                         cf = NULL,
                         aff = NULL,
                         infraspp = NULL) {

  if (!is.data.frame(datset[[1]])) {
    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))) {
      spp_labels <- lapply(datset, function(x) gsub("_cf_", "_cf", names(x)))
      for (i in seq_along(datset)) {
        names(datset[[i]]) <- spp_labels[[i]]
      }
    }

    if (any(unlist(aff))) {
      spp_labels <- lapply(datset, function(x) gsub("_aff_", "_aff", names(x)))
      for (i in seq_along(datset)) {
        names(datset[[i]]) <- spp_labels[[i]]
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      spp_labels <- list()
      for (i in seq_along(datset)) {
        spp_labels[[i]]  <- names(datset[[i]])

        spp_labels[[i]][infraspp[[i]]] <- gsub("(_[^_]+)_", "\\1",
                                               spp_labels[[i]][infraspp[[i]]])
        names(datset[[i]]) <- spp_labels[[i]]
      }
    }

  }

  if (is.data.frame(datset[[1]])) {
    # Adjusting names in dataset already passed through the catmultGenes
    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))){
      spp_labels <- lapply(datset, function(x) gsub("_cf_", "_cf", x[[1]]))
      for (i in seq_along(datset)) {
        datset[[i]][[1]] <- spp_labels[[i]]
      }
    }

    if (any(unlist(aff))) {
      spp_labels <- lapply(datset, function(x) gsub("_aff_", "_aff", x[[1]]))
      for (i in seq_along(datset)) {
        datset[[i]][[1]] <- spp_labels[[i]]
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      spp_labels <- list()
      for (i in seq_along(datset)) {
        spp_labels[[i]]  <- datset[[i]][[1]]

        spp_labels[[i]][infraspp[[i]]] <- sub("(_[^_]+)_", "\\1",
                                              spp_labels[[i]][infraspp[[i]]])
        datset[[i]][[1]] <- spp_labels[[i]]
      }
    }

  }

  return(datset)
}


#-------------------------------------------------------------------------------
# Secondary function to reverse and complement the DNA sequence originally
# Adjust names of taxon, voucher and genbank accessions
# # Used inside mineSeq, minePlastome and mineMitochondrion functions

.tax_voucher_adjust <- function (inputdf = NULL,
                                 taxon = NULL,
                                 voucher = NULL,
                                 genbank = NULL) {

  # Adjust names within mineSeq
  if (!is.null(inputdf) &
      is.null(taxon) &
      is.null(voucher) &
      is.null(genbank)) {

    if ("Species" %in% names(inputdf)) {
      inputdf$Species <- gsub("[.]|(^\\s){1,}|(\\s$){1,}", "", inputdf$Species)
      inputdf$Species <- gsub("(\\s){1,}", "_", inputdf$Species)
      inputdf$Species <- unlist(lapply(inputdf$Species, flora::remove.authors))
    }
    if ("Voucher" %in% names(inputdf)) {
      tf <- is.na(inputdf$Voucher)
      if (any(tf)) {
        inputdf$Voucher[tf] <- "Unvouchered"
      }
      if (length(which(inputdf$Voucher %in% "Unvouchered")) == length(inputdf$Voucher)) {
        inputdf <- inputdf %>% select(-c("Voucher"))
      } else {
        inputdf$Voucher <- gsub("\\s[(].*|[:].*", "", inputdf$Voucher)
        inputdf$Voucher <- gsub("s[.]n[.]", "SN", inputdf$Voucher)
        inputdf$Voucher <- gsub("[/]|\\s|[-]", "", inputdf$Voucher)
        inputdf$Voucher <- gsub(".*[.]\\s|.*[.]", "", inputdf$Voucher)
      }
    }

    # Clean white space across the gene columns
    inputdf[, gb.colnames] <- apply(inputdf[, gb.colnames],
                                    MARGIN = 2, FUN = function(x) gsub("\\s", "", x))

    return(inputdf)
  }

  # Adjust names within minePlastome and mineMitochondrion
  if (is.null(inputdf)) {
    if (!is.null(taxon)) {
      taxon <- gsub("[.]|(^\\s){1,}|(\\s$){1,}", "", taxon)
      taxon <- gsub("(\\s){1,}", "_", taxon)
      taxon <- unlist(lapply(taxon, flora::remove.authors))
    }
    if (!is.null(voucher)) {
      tf <- is.na(voucher)
      if (any(tf)) {
        voucher[tf] <- "Unvouchered"
      }
      if (length(which(voucher %in% "Unvouchered")) == length(voucher)) {
        voucher <- NA
      } else {
        voucher <- gsub("\\s[(].*|[:].*", "", voucher)
        voucher <- gsub("s[.]n[.]", "SN", voucher)
        voucher <- gsub("[/]|\\s|[-]", "", voucher)
        voucher <- gsub(".*[.]\\s|.*[.]", "", voucher)
      }
    }
    genbank <- gsub("\\s", "", genbank)

    return(list(taxon, voucher, genbank))
  }
}


#-------------------------------------------------------------------------------
# Secondary function to reverse and complement the DNA sequence originally
# retrieved from GenBank.
# Adpated from the original function at:
# https://medium.com/biosyntax/reverse-and-find-complement-sequence-in-r-baf33847aab1

# Used inside the function minePlastome

seq_revcompl <- function(seq) {

  alphabets <- strsplit(seq, split = "")[[1]]
  seq <- rev(alphabets)

  # Check if there's "T" in the sequence
  RNA <- Reduce(`|`, seq == "U")
  complvec <- sapply(seq, function(base) {
    # This makes DNA the default
    # As long as there's no U, the sequence is treated as DNA
    if (RNA) {
      switch(base, "A" = "U", "C" = "G", "G" = "C", "U" = "A")
    } else {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
    }
  })
  return(paste(complvec, collapse = ""))
}


#-------------------------------------------------------------------------------
# Auxiliary function to adjust names with cf., aff, or infraspecific taxa
# Used inside catmultGenes and dropSeq functions

.namesback <- function(datset,
                       cf = NULL,
                       aff = NULL,
                       infraspp = NULL,
                       rename_cf = NULL,
                       rename_aff = NULL,
                       rename_infraspp = NULL,
                       shortaxlabel = TRUE,
                       multispp = TRUE) {

  if (multispp == TRUE) {
    # Putting back the names under cf. and aff.
    if (any(unlist(cf))) {
      names_temp_orig <- unique(unlist(rename_cf))
      names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
      names_temp <- unique(gsub("_cf_", "_cf", names_temp_orig))
      n = 0
      for (i in names_temp) {
        n = n + 1
        spp_labels <- lapply(datset, function(x) gsub(i, unique(names_temp_orig)[n], x[[1]]))

        for (j in seq_along(datset)) {
          datset[[j]][[1]] <- spp_labels[[j]]
        }
      }
    }
    if (any(unlist(aff))) {
      names_temp_orig <- unique(unlist(rename_aff))
      names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
      names_temp <- unique(gsub("_aff_", "_aff", names_temp_orig))
      n = 0
      for (i in names_temp) {
        n = n + 1
        spp_labels <- lapply(datset, function(x) gsub(i, unique(names_temp_orig)[n], x[[1]]))

        for (j in seq_along(datset)) {
          datset[[j]][[1]] <- spp_labels[[j]]
        }
      }
    }

    # Adjusting names with infraspecific taxa
    if (any(unlist(infraspp))) {
      names_temp_orig <- unique(unlist(rename_infraspp))
      names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
      names_temp <- unique(gsub("(_[^_]+)_", "\\1", names_temp_orig))
      n = 0
      for (i in names_temp) {
        n = n + 1
        spp_labels <- lapply(datset, function(x) gsub(i, unique(names_temp_orig)[n], x[[1]]))

        for (j in seq_along(datset)) {
          datset[[j]][[1]] <- spp_labels[[j]]
        }
      }
    }
  }

  if (multispp == FALSE) {
    # Putting back the names under cf. and aff.
    if (any(unlist(cf))) {
      names_temp_orig <- unique(unlist(rename_cf))
      if (shortaxlabel) {

        names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
      }
      names_temp <- gsub("_cf_", "_cf", names_temp_orig)
      n = 0
      for (i in names_temp) {
        n = n + 1
        spp_labels <- lapply(datset, function(x) gsub(i, names_temp_orig[n], x[[1]]))

        for (j in seq_along(datset)) {
          datset[[j]][[1]] <- spp_labels[[j]]
        }
      }
    }
    if (any(unlist(aff))) {
      names_temp_orig <- unique(unlist(rename_aff))
      if (shortaxlabel) {
        names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
      }
      names_temp <- gsub("_aff_", "_aff", names_temp_orig)
      n = 0
      for (i in names_temp) {
        n = n + 1
        spp_labels <- lapply(datset, function(x) gsub(i, names_temp_orig[n], x[[1]]))

        for (j in seq_along(datset)) {
          datset[[j]][[1]] <- spp_labels[[j]]
        }
      }
    }

    # Adjusting names with infraspecific taxa
    if (any(unlist(infraspp))) {
      names_temp_orig <- unique(unlist(rename_infraspp))
      if (shortaxlabel) {
        names_temp_orig <- gsub("(_[^_]+_[^_]+)_.*", "\\1", names_temp_orig)
      }
      names_temp <- gsub("(_[^_]+)_", "\\1", names_temp_orig)
      n = 0
      for (i in names_temp) {
        n = n + 1
        spp_labels <- lapply(datset, function(x) gsub(i, names_temp_orig[n], x[[1]]))

        for (j in seq_along(datset)) {
          datset[[j]][[1]] <- spp_labels[[j]]
        }
      }
    }

  }

  return(datset)
}


#-------------------------------------------------------------------------------
# Auxiliary function to get original names with cf., aff., or infraspecific taxa
# Used inside catmultGenes, catfullGenes and dropSeq functions

.namesTorename <- function(datset,
                           cf = NULL,
                           aff = NULL,
                           infraspp = NULL) {

  if (!is.data.frame(datset[[1]])) {
    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))) {
      # Finding species labels to rename
      rename_cf <- list()
      for (i in seq_along(datset)) {
        if (length(names(datset[[i]])[cf[[i]]]) > 0) {
          rename_cf[[i]] <- names(datset[[i]])[cf[[i]]]
          rename_cf <- rename_cf[!is.na(rename_cf)]
        }
      }
    }

    if (any(unlist(aff))) {
      # Finding species labels to rename
      rename_aff <- list()
      for (i in seq_along(datset)) {
        if (length(names(datset[[i]])[aff[[i]]]) > 0) {
          rename_aff[[i]] <- names(datset[[i]])[aff[[i]]]
          rename_aff <- rename_aff[!is.na(rename_aff)]
        }
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      # Finding species labels to rename
      rename_infraspp <- list()
      for (i in seq_along(datset)) {
        if (length(names(datset[[i]])[infraspp[[i]]]) > 0) {
          rename_infraspp[[i]] <- names(datset[[i]])[infraspp[[i]]]
          rename_infraspp <- rename_infraspp[!is.na(rename_infraspp)]
        }
      }
    }

  }

  if (is.data.frame(datset[[1]])) {
    # Adjusting names in dataset already passed through the catmultGenes

    # Adjusting species labels when they have cf or aff
    if (any(unlist(cf))){
      # Finding species labels to rename
      rename_cf <- list()
      for (i in seq_along(datset)) {
        if (length(datset[[i]][[1]][cf[[i]]]) > 0) {
          rename_cf[[i]] <- datset[[i]][[1]][cf[[i]]]
          rename_cf <- rename_cf[!is.na(rename_cf)]
        }
      }
    }

    if (any(unlist(aff))) {
      # Finding species labels to rename
      rename_aff <- list()
      for (i in seq_along(datset)) {
        if (length(datset[[i]][[1]][aff[[i]]]) > 0) {
          rename_aff[[i]] <- datset[[i]][[1]][aff[[i]]]
          rename_aff <- rename_aff[!is.na(rename_aff)]
        }
      }
    }

    # Adjusting species names with infraspecific taxa just for the cross-gene comparisons
    if (any(unlist(infraspp))) {
      # Finding species labels to rename
      rename_infraspp <- list()
      for (i in seq_along(datset)) {
        if (length(datset[[i]][[1]][infraspp[[i]]]) > 0) {
          rename_infraspp[[i]] <- datset[[i]][[1]][infraspp[[i]]]
          rename_infraspp <- rename_infraspp[!is.na(rename_infraspp)]
        }
      }
    }

  }

  orignames <- list(
    if (any(unlist(cf))) {
      rename_cf
    },
    if (any(unlist(aff))) {
      rename_aff
    },
    if (any(unlist(infraspp))) {
      rename_infraspp
    })

  names(orignames) <- c("rename_cf", "rename_aff", "rename_infraspp")

  return(orignames)
}


#-------------------------------------------------------------------------------
# Auxiliary function to insert a scale of character number for each individual
# alignment in the interleaved-formatted concatenated matrix, as well as for the
# entire non-interleaved concatenated matrix, and the associated gene names where
# they starts along the matrix.
# Used inside the function writeNexus

.charScale <- function(x,
                       numbchar = NULL,
                       interleave = NULL,
                       genenames = NULL) {

  if (interleave) {

    for (i in seq_along(x)) {

      # Finding the column with no brackets between accession numbers
      c <- which(grepl("[[]", x[[i]][["sequences"]]) == FALSE)[1]

      if (!is.na(c)) {
        x[[i]] <- rbind(x[[i]][rep(ifelse(c == 0, 1, c), 1), ], x[[i]])
        npad <- nchar(sub("\\s.*", "", x[[i]][1, ]))
        x[[i]][1, ] <- sub(".*?\\s", "", x[[i]][1, ])
        x[[i]][1, ] <- paste0(paste(rep(" ", npad), collapse = ""), x[[i]][1, ])
        x[[i]] <- rbind(x[[i]][rep(1, 1), ], x[[i]])
        x[[i]][1, ] <- gsub("^","[", x[[i]][1, ])
        x[[i]][1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[[i]][1, ])

        x[[i]][2, ] <- gsub("^","[", x[[i]][2, ])
        x[[i]][2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[[i]][2, ])

      } else {

        x[[i]] <- rbind(x[[i]][rep(1, 1), ], x[[i]])

        npad_acc <- gsub("\\s{2}.*", "", x[[i]][1, ])
        npad_acc <- gsub(".*\\s", "", npad_acc)

        x[[i]][1, ] <- gsub("[[]\\w+[]]",
                            paste(rep(" ", nchar(npad_acc)), collapse = ""),
                            x[[i]][1, ])

        npad <- nchar(sub("\\s.*", "", x[[i]][1, ]))
        x[[i]][1, ] <- sub(".*?\\s", "", x[[i]][1, ])
        x[[i]][1, ] <- paste0(paste(rep(" ", npad), collapse = ""), x[[i]][1, ])
        x[[i]] <- rbind(x[[i]][rep(1, 1), ], x[[i]])
        x[[i]][1, ] <- gsub("^","[", x[[i]][1, ])
        x[[i]][1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[[i]][1, ])

        x[[i]][2, ] <- gsub("^","[", x[[i]][2, ])
        x[[i]][2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[[i]][2, ])
      }
    }

    nbr <- list()
    nbrseq <- list()
    dotseq <- list()
    for (j in seq_along(numbchar)) {
      nbr[[j]] <- paste(seq(0, numbchar[[j]]+20, by = 10))
      tempA <- list()
      tempB <- list()
      for (i in seq_along(nbr[[j]])) {
        if (nbr[[j]][[i]] == 0) {
          tempA[[i]] <- paste(1, paste(rep(" ", 8-nchar(nbr[[j]][[i]])), collapse = ""), collapse = "")
          tempB[[i]] <- paste(".", paste(rep(" ", 8-nchar(nbr[[j]][[1]])), collapse = ""), collapse = "")
        } else {
          tempA[[i]] <- paste(nbr[[j]][[i]], paste(rep(" ", 9-nchar(nbr[[j]][[i]])), collapse = ""), collapse = "")
          tempB[[i]] <- paste(".", paste(rep(" ", 9-nchar(nbr[[j]][[1]])), collapse = ""), collapse = "")
        }
      }
      nbrseq[[j]] <- unlist(tempA)
      dotseq[[j]] <- unlist(tempB)

      x[[j]][1, ] <- paste0(x[[j]][1, ], paste0(paste(nbrseq[[j]], collapse = ""), "]"))
      x[[j]][2, ] <- paste0(x[[j]][2, ], paste0(paste(dotseq[[j]], collapse = ""), "]"))
    }

  } else {

    # Creating the scale of char length for the non-interleaved concatenated matrix

    # Finding the column with no brackets between accession numbers
    c <- which(grepl("[[]", x[["sequences"]]) == FALSE)[1]

    if (!is.na(c)) {

      x <- rbind(x[rep(ifelse(c == 0, 1, c), 1), ], x)
      npad <- nchar(sub("\\s.*", "", x[1, ]))
      x[1, ] <- sub(".*?\\s", "", x[1, ])
      x[1, ] <- paste0(paste(rep(" ", npad), collapse = ""), x[1, ])
      x <- rbind(x[rep(1, 1), ], x)
      x <- rbind(x[rep(1, 1), ], x)
      x[3, ] <- gsub("^","[", x[3, ])
      x[3, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[3, ])

      x[2, ] <- gsub("^","[", x[2, ])
      x[2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[2, ])

      x[1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[1, ])
      x[1, ] <- gsub("^"," ", x[1, ])

    } else {

      x <- rbind(x[rep(1, 1), ], x)

      npad_acc <- gsub("\\s{2}.*", "", x[1, ])
      npad_acc <- gsub(".*\\s", "", npad_acc)

      x[1, ] <- gsub("[[]\\w+[]]",
                     paste(rep(" ", nchar(npad_acc)), collapse = ""),
                     x[1, ])

      npad <- nchar(sub("\\s.*", "", x[1, ]))
      x[1, ] <- sub(".*?\\s", "", x[1, ])
      x[1, ] <- paste0(paste(rep(" ", npad), collapse = ""), x[1, ])
      x <- rbind(x[rep(1, 1), ], x)
      x <- rbind(x[rep(1, 1), ], x)
      x[3, ] <- gsub("^","[", x[3, ])
      x[3, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[3, ])

      x[2, ] <- gsub("^","[", x[2, ])
      x[2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[2, ])

      x[1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", x[1, ])
      x[1, ] <- gsub("^"," ", x[1, ])
    }

    s <- sum(unlist(numbchar))
    nbr <- paste(seq(0, s+20, by = 10))
    nbrseq <- list()
    dotseq <- list()
    for (i in seq_along(nbr)) {
      if (nbr[[i]] == 0) {
        nbrseq[[i]] <- paste(1, paste(rep(" ", 8-nchar(nbr[[i]])), collapse = ""), collapse = "")
        dotseq[[i]] <- paste(".", paste(rep(" ", 8-nchar(nbr[[1]])), collapse = ""), collapse = "")
      } else {
        nbrseq[[i]] <- paste(nbr[[i]], paste(rep(" ", 9-nchar(nbr[[i]])), collapse = ""), collapse = "")
        dotseq[[i]] <- paste(".", paste(rep(" ", 9-nchar(nbr[[1]])), collapse = ""), collapse = "")
      }
    }

    x[2, ] <- paste0(x[2, ], paste0(paste(unlist(nbrseq), collapse = ""), "]"))
    x[3, ] <- paste0(x[3, ], paste0(paste(unlist(dotseq), collapse = ""), "]"))

    g <- paste0("[", as.list(genenames), "]")

    geneseq <- list()
    for (i in seq_along(g)) {
      if (i == 1) {
        geneseq[[i]] <- paste(g[i], paste(rep(" ", numbchar[[i]]-nchar(g[i])-1),
                                          collapse = ""))
      } else {
        geneseq[[i]] <- paste(g[i], paste(rep(" ", numbchar[[i]]-nchar(g[i])-1),
                                          collapse = ""), collapse = "")
      }
    }

    x[1, ] <- paste0(x[1, ], paste(unlist(geneseq), collapse = ""))

  }

  return(x)
}

#-------------------------------------------------------------------------------
# Auxiliary function to replace terminal GAPs into missing character (?)
# Used inside the function nexusdframe writeNexus

.replace_terminal_gaps <- function (x) {
  tf <- grepl("^-", x$sequence)
  if (any(tf)) {
    temp <- gsub("[[:upper:]].*", "", x$sequence[tf])
    ll <- unlist(lapply(temp, function(y) length(unlist(strsplit(y, "")))))
    for (i in seq_along(temp)) {
      x$sequence[tf][i] <- gsub(paste0("^", temp[i]),
                                paste0(rep("?", ll[i]), collapse = ""),
                                x$sequence[tf][i])
    }
  }
  tf <- grepl("-$", x$sequence)
  if (any(tf)) {
    temp <- gsub(".*[[:upper:]]", "", x$sequence[tf])
    ll <- unlist(lapply(temp, function(y) length(unlist(strsplit(y, "")))))

    for (i in seq_along(temp)) {
      x$sequence[tf][i] <- gsub(paste0(temp[i], "$"),
                                paste0(rep("?", ll[i]), collapse = ""),
                                x$sequence[tf][i])
    }
  }
  return(x)
}

