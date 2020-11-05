## code to prepare `DATASET` dataset goes here

library(ape)
# Load each example DNA alignments
# This will require the ape function read.nexus.data

# Loading Luetzelburgia example
genes <- list.files("data-raw/DNAlignments/Luetzelburgia")
Luetzelburgia <- list()
for (i in genes){
  Luetzelburgia[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Luetzelburgia/",i))
}
names(Luetzelburgia) <- gsub(".*_(.+)[.].*", "\\1", names(Luetzelburgia))

# Adding dataset for tests
usethis::use_data(Luetzelburgia, overwrite = TRUE)


# Loading Ormosia example
genes <- list.files("data-raw/DNAlignments/Ormosia")
Ormosia <- list()
for (i in genes){
  Ormosia[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Ormosia/",i))
}
names(Ormosia) <- gsub(".*_(.+)[.].*", "\\1", names(Ormosia))

# Adding dataset for tests
usethis::use_data(Ormosia, overwrite = TRUE)


# Loading Vataireoids example
genes <- list.files("data-raw/DNAlignments/Vataireoids")
Vataireoids <- list()
for (i in genes){
  Vataireoids[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Vataireoids/",i))
}
names(Vataireoids) <- gsub("[.].*", "", names(Vataireoids))

# Adding dataset for tests
usethis::use_data(Vataireoids, overwrite = TRUE)


# Loading Vataireoids example
genes <- list.files("data-raw/DNAlignments/Gaya")
Gaya <- list()
for (i in genes){
  Gaya[[i]] <- ape::read.nexus.data(paste0("data-raw/DNAlignments/Gaya/",i))
}
names(Gaya) <- gsub("[.].*", "", names(Gaya))

# Adding dataset for tests
usethis::use_data(Gaya, overwrite = TRUE)
