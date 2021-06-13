## code to prepare `DATASET` dataset goes here

library(ape)
# Load each example DNA alignments
# This will require the ape function read.nexus.data

# Loading Luetzelburgia example
genes <- list.files("data-raw/DNAlignments/Luetzelburgia")
Luetzelburgia <- list()
for (i in genes) {
  Luetzelburgia[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Luetzelburgia/",i))
}
names(Luetzelburgia) <- gsub(".*_(.+)[.].*", "\\1", names(Luetzelburgia))

# Adding dataset for tests
usethis::use_data(Luetzelburgia, overwrite = TRUE)


# Loading Ormosia example
genes <- list.files("data-raw/DNAlignments/Ormosia")
Ormosia <- list()
for (i in genes) {
  Ormosia[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Ormosia/",i))
}
names(Ormosia) <- gsub(".*_(.+)[.].*", "\\1", names(Ormosia))

# Adding dataset for tests
usethis::use_data(Ormosia, overwrite = TRUE)


# Loading Vataireoids example
genes <- list.files("data-raw/DNAlignments/Vataireoids")
Vataireoids <- list()
for (i in genes) {
  Vataireoids[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Vataireoids/",i))
}
names(Vataireoids) <- gsub("[.].*", "", names(Vataireoids))

# Adding dataset for tests
usethis::use_data(Vataireoids, overwrite = TRUE)


# Loading Gaya example
genes <- list.files("data-raw/DNAlignments/Gaya")
Gaya <- list()
for (i in genes) {
  Gaya[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Gaya/",i))
}
names(Gaya) <- gsub("[.].*", "", names(Gaya))

# Adding dataset for tests
usethis::use_data(Gaya, overwrite = TRUE)


# Loading Brongniartia example
genes <- list.files("data-raw/DNAlignments/Brongniartia")
Brongniartia <- list()
for (i in genes) {
  Brongniartia[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Brongniartia/",i))
}
names(Brongniartia) <- gsub("[.].*", "", names(Brongniartia))

# Adding dataset for tests
usethis::use_data(Brongniartia, overwrite = TRUE)


# Loading Cryptocarya example
genes <- list.files("data-raw/DNAlignments/Cryptocarya")
Cryptocarya <- list()
for (i in seq_along(genes)) {

  Cryptocarya[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Cryptocarya/", genes[i]))

}
names(Cryptocarya) <- gsub("[.].*", "", genes)

# Adding dataset for tests
usethis::use_data(Cryptocarya, overwrite = TRUE)


# Loading Vochysiaceae example
# See that the original DNA alignments are in Phylip format.
# So we will first import into R and them save in Nexus format
genes <- list.files("data-raw/DNAlignments/Vochysiaceae")
Vochysiaceae <- list()
for (i in genes) {
  chopper::alg2nex(paste0("data-raw/DNAlignments/Vochysiaceae/",i),
                   format = "sequential",
                   interleaved = FALSE, gap = "-",
                   missing = "?", partition.file = NULL)
}

f <- list.files("data-raw/DNAlignments/Vochysiaceae")
genes <- f[grepl("[.]nex", f)]
Vochysiaceae <- list()
for (i in genes) {
  Vochysiaceae[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Vochysiaceae/",i))

  #Vochysiaceae[[i]] <- phylotools::read.phylip(paste0("data-raw/DNAlignments/Vochysiaceae/",i))
  # renaming the terminals
  hitaxa <- gsub("(_[^_]+)_.*", "\\1", names(Vochysiaceae[[i]]))
  taxa_temp <- sub(".*?_", "",  names(Vochysiaceae[[i]]))
  taxa <- sub(".*?_", "",  taxa_temp)
  names(Vochysiaceae[[i]]) <- paste(taxa, hitaxa, sep="_")
  names(Vochysiaceae[[i]]) <- gsub("__", "_", names(Vochysiaceae[[i]]))
}
names(Vochysiaceae) <- gsub("[.].*", "", names(Vochysiaceae))

# Adding dataset for tests
usethis::use_data(Vochysiaceae, overwrite = TRUE)

