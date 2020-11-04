##############################################################################
###  Script to run gene concatenation functions from the catGenes package  ###
##############################################################################

# Authors: Domingos Cardoso
# Modified: 11 Oct 2020


################################################
###  IMPORTANT: before running the catGenes  ###
################################################

# It will works only for DNA aligments that is perfectly format like the following:
# Genus_species_everythingelse    #so... underscores separate the genus and species name from the collector+number+genbank accession
# G_species_everythingelse        #so... no problem if the genus name is abbreviated


################################################################
####  Reading the DNA alignments for concatenated analysis  ####
################################################################

library(ape)

# Install catGenes from GitHub
library(devtools)
# https://github.com/domingoscardoso/catGenes
devtools::install_github("domingoscardoso/catGenes")
library(catGenes)


################################################################################
# Load each gene dataset from the working directory
# This will require the ape function read.nexus.data

list.files("inst/DNAlignments")
ITS <- read.nexus.data("inst/DNAlignments/Gaya_ITS.nex")
petLpsbE <- read.nexus.data("inst/DNAlignments/Gaya_petL_psbE.nex")
rpl16 <- read.nexus.data("inst/DNAlignments/Gaya_rpl16.nex")


# Loading all DNA alignments from the working directory
dat_names <- list.files("inst/DNAlignments/Luetzelburgia")
Luetzelburgia <- list()
for (i in dat_names){
  Luetzelburgia[[i]] <- read.nexus.data(paste0("inst/DNAlignments/Luetzelburgia/",i))
}
names(Luetzelburgia) <- gsub("Luetz_", "", names(Luetzelburgia))
names(Luetzelburgia) <- gsub("[.].+", "", names(Luetzelburgia))
names(Luetzelburgia) <- gsub("_", "", names(Luetzelburgia))
# Adding dataset for tests
save(Luetzelburgia, file = "Luetzelburgia.rda")


# Loading all DNA alignments from the working directory
dat_names <- list.files("inst/DNAlignments/Gaya")
Gaya <- list()
for (i in dat_names){
  Gaya[[i]] <- read.nexus.data(paste0("inst/DNAlignments/Gaya/",i))
}
# Adding dataset for tests
save(Gaya, file = "Gaya.rda")


dat_names <- list.files("inst/DNAlignments/Vataireoids")
Vataireoids <- list()
for (i in dat_names){
  Vataireoids[[i]] <- read.nexus.data(paste0("inst/DNAlignments/Vataireoids/",i))
}


dat_names <- list.files("inst/DNAlignments/Ormosia")
Ormosia <- list()
for (i in dat_names){
  Ormosia[[i]] <- read.nexus.data(paste0("inst/DNAlignments/Ormosia/",i))
}
names(Ormosia) <- gsub("Ormosia_", "", names(Ormosia))
names(Ormosia) <- gsub("[.].+", "", names(Ormosia))
names(Ormosia) <- gsub("_", "", names(Ormosia))


# Adding dataset for tests
save(Ormosia, file = "Ormosia.rda")
################################################################################
################################################################################





################################################################################
# Load each gene dataset from the catGenes directory
# This will require the ape function read.nexus.data

list.files(system.file("DNAlignments/Gaya", package = "catGenes"))
ITS <- read.nexus.data(system.file("DNAlignments/Gaya", "ITS.nex", package = "catGenes"))
petLpsbE <- read.nexus.data(system.file("DNAlignments/Gaya", "petLpsbE.nex", package = "catGenes"))
rpl16 <- read.nexus.data(system.file("DNAlignments/Gaya", "rpl16.nex", package = "catGenes"))

# Loading all DNA alignments from within the catGenes directory
dat_names <- list.files(system.file("DNAlignments/Luetzelburgia", package = "catGenes"))
Luetzelburgia <- list()
for (i in dat_names[grepl("Luetz", dat_names)]){
  Luetzelburgia[[i]] <- read.nexus.data(system.file("DNAlignments/Luetzelburgia", i, package = "catGenes"))
}
names(Luetzelburgia) <- gsub("Luetz_", "", names(Luetzelburgia))
names(Luetzelburgia) <- gsub("[.].+", "", names(Luetzelburgia))
names(Luetzelburgia) <- gsub("_", "", names(Luetzelburgia))

################################################################################




##################################################################
#######################  Running catGenes  #######################
##################################################################

################################################################################
# You can run catGenes using an application by just running the function below

run_catGenes()

################################################################################



################################################################################
# You can run catGenes by running the specific functions in the code belo

###################################################################
###  "catfullGenes" function to match individual gene datasets  ###
###################################################################

# This newly developed "catfullGenes" function has the following arguments:
# "shortaxlabel": if TRUE it will keep the terminals/taxa with just the scientific name
# "missdata": if TRUE it will include missing data in the resulting dataset
# "outgroup": to define one or a list of taxon/terminal as outgroup; this argument is more useful when cocatenating genes that fully match and you know that one outgroup sequence is absent in one of the genes being compared. so, outgroup must be defined as because you do not want to miss the outgroup sequence
# "multiaccessions": if TRUE the analysis will consider that you have species represented
# with multiple accessions, which means you after the taxon names the collector+number must
# have to be prefectly formatted, because the functions will use de combination of
# taxon + collector + number as identifiers to compare among the input datasets.


################################################################################
################################################################################
# For combined analysis to get all individual dataset with full match...
# This means the function will find the intersection between the genes under comparison

# IMPORTANT: if you know that one outgroup is absent in any gene/DNa alignment...
# then outgroup must be defined as one or a list of taxon/terminal in your DNA alignments
# the resulting dataset of this analysis will always have the outgroups sorted alphabetically on top of the DNA aligments following by the ingroup sequences also sorted alphabetically

# We can use use a list of gene datasets as input
# Run this code when sequence of outgroup taxa is present in all DNA aligments
cat_datasets <- catfullGenes(ITS, petLpsbE, rpl16,
                             shortaxlabel = TRUE,
                             missdata = FALSE)

# When you do not have sequence of outgroup taxa for all genes/DNA alignments
cat_datasets <- catfullGenes(Ormosia,
                             shortaxlabel = TRUE,
                             missdata = FALSE,
                             outgroup = "Parasenegalia_muricata")


# When you have more that than one outgroup sequence...
# Define outgroup terminals as a list
outgrouptaxa <- c("Abutilon_costicalyx",
                  "Abutilon_itatiaie")
cat_datasets <- catfullGenes(DNAlignments,
                             shortaxlabel = TRUE,
                             missdata = FALSE,
                             multiaccessions = FALSE,
                             outgroup = outgrouptaxa)




################################################################################
################################################################################
# For combined analysis including "missing data" in the resulting dataset

# This means the species without a sequence or absent in one gene will be represented by missing data or "?"
# This is actually the case of incomplete taxon...
# See also the follwing papers:
# Wiens (2003) Systematic Biology 52(4):528-538; https://doi.org/10.1080/10635150390218330
# Wiens (2006) Journal of Biomedical Informatics 39(1):34-42; https://doi.org/10.1016/j.jbi.2005.04.001
# Wiens & Morrill (2011) Systematic Biology 60(5):719-731; https://doi.org/10.1093/sysbio/syr025

# You will see that it does not matter if you provide an outgroup or not... or if the outgroup is wrongly named...
# But note that the resulting dataset will have sequences sorted alphabetically by taxon/terminal names


cat_datasets1 <- catfullGenes(Vataireoids,
                             shortaxlabel = TRUE,
                             missdata = TRUE)

cat_datasets2 <- catfullGenes(Ormosia,
                             shortaxlabel = TRUE,
                             missdata = FALSE)


cat_datasets <- catmultGenes(Ormosia,
                             maxspp = TRUE,
                             shortaxlabel = TRUE,
                             missdata = TRUE)

# This is just a test with a wrongly defined outgroup taxon/terminal
cat_datasets <- catfullGenes(ITS, petLpsbE, rpl16,
                              shortaxlabel = TRUE,
                              missdata = TRUE,
                              multiaccessions = FALSE,
                              outgroup = "wrong_outgroup_JM1066")









#################################################################################
###  "writeNexus" and "writePhylip" functions for concatenating the datasets  ###
#################################################################################

# After running the "genefullcomp"...
# Run the output with the newly developed functions "writeNexus" or "writePhylip"
# These functions generate the concatenated dataset of all matched genes!

# If the argument bayesblock is TRUE...
# Then "writeNexus" will insert a preliminary Mr.Bayes command block
# Where the charset for each individual gene dataset is defined
# You will see several bracketed uppercase comments...
# where you have will have to edit the Mr.Bayes parameters accordingly


# Generating a nexus-formatted concatenated dataset without Mr.Bayes command block
writeNexus(cat_datasets1, file = "Vataireoids1.nex",
           bayesblock = TRUE, interleave = TRUE)

# Generating a nexus-formatted  concatenated dataset in Nexus format with Mr.Bayes command block
writeNexus(cat_datasets, file = "Gaya_missingdata_non_interleave.nex",
           bayesblock = TRUE, interleave = FALSE)

writeNexus(cat_datasets, file = "Luetz_missingdata_interleave.nex",
           bayesblock = TRUE, interleave = TRUE)


# Generating a concatenated dataset in Phylip format
writePhylip(cat_datasets, file = "Luetz_dataset.phy",
            catalignments = TRUE,
            partitionfile = TRUE)



################################################################################





################################################################################
# Dropping duplicated sequences

final_datasets <- dropSeq(cat_datasets)


# Generating a nexus-formatted concatenated dataset without Mr.Bayes command block
writeNexus(final_datasets, file = "Ormosia_missingdata_interleave_without_commands.nex",
           bayesblock = FALSE, interleave = TRUE)

# Generating a nexus-formatted  concatenated dataset in Nexus format with Mr.Bayes command block
writeNexus(final_datasets, file = "Ormosia_missingdata_non_interleave.nex",
           bayesblock = TRUE, interleave = FALSE)

writeNexus(final_datasets, file = "Ormosia_missingdata_interleave.nex",
           bayesblock = TRUE, interleave = TRUE)


# Generating a concatenated dataset in Phylip format
writePhylip(final_datasets, file = "Ormosia_dataset.phy",
            catalignments = TRUE,
            partitionfile = TRUE)



######################################################
###  To save individual datasets after gene match  ###
######################################################

# After running the gene match...
# If you want, you can also save the indivudal datasets using the following newly developed functions
# "fastadframe" - write each indivual gene (or dataframe formatted dataset) into fasta format
# "nexusdframe" - save each individual gene (or dataframe formatted dataset) into in nexus format


# Let's first extract each matched dataset to individual dataframes

ITS <- final_datasets[[1]]
matK <- final_datasets[[2]]
trnKintron <- final_datasets[[3]]
trnLF <- final_datasets[[4]]

##############################################################################
###  "nexusdframe" writes a dataframe formatted dataset into nexus format  ###
##############################################################################

nexusdframe(ITS, "ITS_Ormosia_comb.nex",
            dropmisseq = F)
nexusdframe(matK, "matK_Ormosia_comb.nex",
            dropmisseq = TRUE)
nexusdframe(trnKintron, "trnKintron_Ormosia_comb.nex",
            dropmisseq = TRUE)
nexusdframe(trnLF, "trnLF_Ormosia_comb.nex",
            dropmisseq = TRUE)


##############################################################################
###  "phylipdframe" writes a dataframe formatted dataset into phylip format  ###
##############################################################################

phylipdframe(ITS, "ITS_Ormosia_comb.phy",
            dropmisseq = TRUE)
phylipdframe(matK, "matK_Ormosia_comb.phy",
            dropmisseq = TRUE)
phylipdframe(trnKintron, "trnKintron_Ormosia_comb.phy",
            dropmisseq = TRUE)
phylipdframe(trnLF, "trnLF_Ormosia_comb.phy",
            dropmisseq = TRUE)



##############################################################################
###  "fastadframe" writes a dataframe formatted dataset into fasta format  ###
##############################################################################

fastadframe(ITS, "ITS_Ormosia_comb.fasta",
            dropmisseq = TRUE)
fastadframe(matK, "matK_Ormosia_comb.fasta",
            dropmisseq = TRUE)
fastadframe(trnKintron, "trnKintron_Ormosia_comb.fasta",
            dropmisseq = TRUE)
fastadframe(trnLF, "trnLF_Ormosia_comb.fasta",
            dropmisseq = TRUE)







################################################################################
# FORGET ABOUT THIS CODE BELOW

# How to build an R package
# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html


# To load all new functions in the R folder
devtools::document()

# Adding a dataset for tests
save(DNAlignments, file = "DNAlignments.rda")




################################################################################
################################################################################
library(read.gb)
# Under development
# I will be working on new R functions to retrieve specific DNA sequences from GenBank
GenBank_data <- read.gb("Data/Luetz_sequences.gb")
GenBank_data1 <- readChar("Data/Luetz_sequences.gb", file.info("Data/Luetz_sequences.gb")$size, nchars = 99999999)
################################################################################

