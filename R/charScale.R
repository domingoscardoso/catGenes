# Auxiliary function to insert a scale of character number for each individual alignment

# Used inside the function writeNexus

# Author: Domingos Cardoso


.charScale <- function(datset,
                       numbchar = NULL,
                       interleave = NULL) {

  if (interleave) {

    for (i in seq_along(datset)) {

      # Finding the column with no brackets between accession numbers
      c <- which(grepl("[[]", datset[[i]][["sequences"]]) == FALSE)[1]

      if (!is.na(c)) {
        datset[[i]] <- rbind(datset[[i]][rep(ifelse(c == 0, 1, c), 1), ], datset[[i]])
        npad <- nchar(sub("\\s.*", "", datset[[i]][1, ]))
        datset[[i]][1, ] <- sub(".*?\\s", "", datset[[i]][1, ])
        datset[[i]][1, ] <- paste0(paste(rep(" ", npad), collapse = ""), datset[[i]][1, ])
        datset[[i]] <- rbind(datset[[i]][rep(1, 1), ], datset[[i]])
        datset[[i]][1, ] <- gsub("^","[", datset[[i]][1, ])
        datset[[i]][1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", datset[[i]][1, ])

        datset[[i]][2, ] <- gsub("^","[", datset[[i]][2, ])
        datset[[i]][2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", datset[[i]][2, ])

      } else {

        datset[[i]] <- rbind(datset[[i]][rep(1, 1), ], datset[[i]])

        npad_acc <- gsub("\\s{2}.*", "", datset[[i]][1, ])
        npad_acc <- gsub(".*\\s", "", npad_acc)

        datset[[i]][1, ] <- gsub("[[]\\w+[]]",
                                 paste(rep(" ", nchar(npad_acc)), collapse = ""),
                                 datset[[i]][1, ])


        npad <- nchar(sub("\\s.*", "", datset[[i]][1, ]))
        datset[[i]][1, ] <- sub(".*?\\s", "", datset[[i]][1, ])
        datset[[i]][1, ] <- paste0(paste(rep(" ", npad), collapse = ""), datset[[i]][1, ])
        datset[[i]] <- rbind(datset[[i]][rep(1, 1), ], datset[[i]])
        datset[[i]][1, ] <- gsub("^","[", datset[[i]][1, ])
        datset[[i]][1, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", datset[[i]][1, ])

        datset[[i]][2, ] <- gsub("^","[", datset[[i]][2, ])
        datset[[i]][2, ] <- gsub("[[:upper:]].*|[[:lower:]].*|[?]|[-]", "", datset[[i]][2, ])
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
          tempA[[i]] <- paste(nbr[[j]][[i]], paste(rep(" ", 8-nchar(nbr[[j]][[i]])), collapse = ""), collapse = "")
          tempB[[i]] <- paste(".", paste(rep(" ", 8-nchar(nbr[[j]][[1]])), collapse = ""), collapse = "")
        } else {
          tempA[[i]] <- paste(nbr[[j]][[i]], paste(rep(" ", 9-nchar(nbr[[j]][[i]])), collapse = ""), collapse = "")
          tempB[[i]] <- paste(".", paste(rep(" ", 9-nchar(nbr[[j]][[1]])), collapse = ""), collapse = "")
        }
      }
      nbrseq[[j]] <- unlist(tempA)
      dotseq[[j]] <- unlist(tempB)

      datset[[j]][1, ] <- paste0(datset[[j]][1, ], paste0(paste(nbrseq[[j]], collapse = ""), "]"))
      datset[[j]][2, ] <- paste0(datset[[j]][2, ], paste0(paste(dotseq[[j]], collapse = ""), "]"))
    }

    #} else {
    #number of chacter for non interleave
    #}


  }

  return(datset)
}
