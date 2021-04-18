# Auxiliary function to insert a scale of character number for each individual alignment

# Used inside the function writeNexus

# Author: Domingos Cardoso


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
          tempA[[i]] <- paste(nbr[[j]][[i]], paste(rep(" ", 8-nchar(nbr[[j]][[i]])), collapse = ""), collapse = "")
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
        nbrseq[[i]] <- paste(nbr[[i]], paste(rep(" ", 8-nchar(nbr[[i]])), collapse = ""), collapse = "")
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
