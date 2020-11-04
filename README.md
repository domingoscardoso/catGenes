
<!-- README.md is generated from README.Rmd. Please edit that file -->

# catGenes

<!-- badges: start -->

<!-- badges: end -->

The catGenes package is intended to help researchers in phylogenetics
and phylogenomics to build a fully concatenated or combined
(non-interleaved) dataset by automatically comparing individual DNA
alignments.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("domingoscardoso/catGenes")
```

## Standard format for input DNA alignments

#### *Naming the individual DNA alignment files*

Before loading the individual DNA alignments for concatenation, we
recommend that each file is simply named with the corresponding gene
names and do not use space or hyphen to separate the file name. For
example: **ITS.nex**, **rbcL.nex**, **psbAtrnH.nex**, **COX1.nex**, etc.
This is the best way to name the input files if you want to run
*catGenes* using the Shiny app. Also, it facilitates to load all DNA
alignments by just using a for loop that will create a list of the input
genes ready for the concatenation. Note that we load the DNA alignments
by using *ape*’s function **read.nexus.data**. See below an example
using the available DNA alignments of the Vataireoid legumes that are
stored within the *catGenes* directory:

``` r
library(catGenes)

# Loading all individual DNA alignments for the concatenation
genes <- list.files(system.file("DNAlignments/Vataireoids", 
                                package = "catGenes"))
Vataireoids <- list()
for (i in genes){
  Vataireoids[[i]] <- ape::read.nexus.data(system.file("DNAlignments/Vataireoids", i, 
                                                       package = "catGenes"))
}
names(Vataireoids) <- gsub("[.].*", "", names(Vataireoids))
```

The created object **Vataireoids** is a list of all input individual DNA
alignments ready for the concatenation by using the available *catGenes*
functions. As you can see below, the containing input DNA alignments
within the Vataireoids list are named exactly how the files were
originally named:

``` r
names(Vataireoids)
#> [1] "ETS"   "ITS"   "matK"  "psbA"  "rpS16" "trnDT" "trnL"  "trnQ"
```

Note that in your computer the individual DNA alignments could be loaded
into a single list by adjusting the abovementioned code as follows:

``` r
library(catGenes)

genes <- list.files("path_to_DNA_alignments_folder")
Vataireoids <- list()
for (i in genes){
  Vataireoids[[i]] <- ape::read.nexus.data(paste0("path_to_DNA_alignments_folder/", i))
}
names(Vataireoids) <- gsub("[.].*", "", names(Vataireoids))
```

If the input individual DNA alignments have any identifier associated to
the gene names, like **Vataireoids\_ITS.nex**,
**Vataireoids\_matK.nex**, **Vataireoids\_trnDT.nex**, etc., then you
can load all such files by either importing each DNA alignment
separately or by adjusting the abovementioned for loop as follows:

``` r
library(catGenes)

# Loading each DNA alignment from the working directory separately by using the
# ape's function read.nexus.data
ITS <- read.nexus.data("Vataireoids_ITS.nex")
matK <- read.nexus.data("Vataireoids_matK.nex")
trnDT <- read.nexus.data("Vataireoids_trnDT.nex")

# Loading just the Vataireoid DNA alignments even if there are other DNA alignments 
# in the same working directory
genes <- list.files("path_to_DNA_alignments_folder")
Vataireoids <- list()
for (i in genes[grepl("Vataireoids", genes)]){
  Vataireoids[[i]] <- ape::read.nexus.data(paste0("path_to_DNA_alignments_folder/", i))
}
# Deleting all characters between "Vataireoids_" and ".nex" so as to keep just the
# name of the genes in the vataireoids list of named DNA alignments
names(Vataireoids) <- gsub(".*_(.+)[.].*", "\\1", names(Vataireoids))
```

#### *Species labels and associated identifiers within the DNA alignments*

Before loading the DNA alignments to run catGenes functions, make sure
the DNA aligments have the species names format as following:

1)  Underscores separate the genus and species name from the
    collector+number+genbank accession

Genus\_species\_everythingelse

2)  so… no problem if the genus name is abbreviated

G\_species\_everythingelse

``` r
names(Luetzelburgia$ITS)
#>  [1] "Vataireopsis_araroba_Cardoso2175"      
#>  [2] "Vataireopsis_speciosa_Cardoso2208"     
#>  [3] "Vatairea_heteroptera_Cardoso2205"      
#>  [4] "Vatairea_macrocarpa_Cardoso2388"       
#>  [5] "Luetzelburgia_amazonica_Cardoso2920"   
#>  [6] "Luetzelburgia_amazonica_Maciel1313"    
#>  [7] "Luetzelburgia_andina_Cayola2358"       
#>  [8] "Luetzelburgia_andina_Michel4548"       
#>  [9] "Luetzelburgia_andradelimae_Cardoso1848"
#> [10] "Luetzelburgia_andradelimae_Cardoso2518"
#> [11] "Luetzelburgia_andradelimae_Cardoso3000"
#> [12] "Luetzelburgia_auriculata_Cardoso2372"  
#> [13] "Luetzelburgia_auriculata_Cardoso2385"  
#> [14] "Luetzelburgia_auriculata_Cardoso2617"  
#> [15] "Luetzelburgia_auriculata_Cardoso2620"  
#> [16] "Luetzelburgia_auriculata_Cardoso2667"  
#> [17] "Luetzelburgia_auriculata_Faria959"     
#> [18] "Luetzelburgia_auriculata_Chagas109"    
#> [19] "Luetzelburgia_bahiensis_Cardoso2690"   
#> [20] "Luetzelburgia_bahiensis_Cardoso2814"   
#> [21] "Luetzelburgia_bahiensis_Cardoso2942"   
#> [22] "Luetzelburgia_bahiensis_Cardoso2970"   
#> [23] "Luetzelburgia_guaissara_Cardoso2202"   
#> [24] "Luetzelburgia_guaissara_Cardoso2213"   
#> [25] "Luetzelburgia_guaissara_Reitz6384"     
#> [26] "Luetzelburgia_guianensis_Cardoso2853"  
#> [27] "Luetzelburgia_guianensis_Cardoso2857"  
#> [28] "Luetzelburgia_harleyi_Cardoso2076"     
#> [29] "Luetzelburgia_harleyi_Cardoso2186"     
#> [30] "Luetzelburgia_harleyi_Cardoso2187"     
#> [31] "Luetzelburgia_harleyi_Cardoso2190"     
#> [32] "Luetzelburgia_jacana_Banda386"         
#> [33] "Luetzelburgia_neurocarpa_Cardoso1853"  
#> [34] "Luetzelburgia_praecox_Cardoso2540"     
#> [35] "Luetzelburgia_praecox_Cardoso2545"     
#> [36] "Luetzelburgia_praecox_Ratter7975v"     
#> [37] "Luetzelburgia_purpurea_Cardoso1801"    
#> [38] "Luetzelburgia_purpurea_Cardoso2103"    
#> [39] "Luetzelburgia_sotoi_Soto1135"          
#> [40] "Luetzelburgia_sotoi_Soto1167"          
#> [41] "Luetzelburgia_sotoi_Soto1210"          
#> [42] "Luetzelburgia_trialata_Cardoso2196"    
#> [43] "Luetzelburgia_trialata_Moreno70"
```
