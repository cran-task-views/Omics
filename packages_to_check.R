# useful packages
library(tidytext)
library(tidyverse)
library(dplyr)
library(lubridate)
library(ggplot2)
library(SnowballC) # useful for stemming

# keywords for genetics and genomics data
omics_words <- c(
  "allele", "dna", "eqtl", "gene", "genes", "genetic", "genetics", "genetically", 
  "genome", "genomic", "genomics", "genotype", "genotypes", "gwas", "metabolomic",
  "metabolomics", "metagenomic", "metagenomics", "nucleotide", "omic", "omics",
  "protein", "proteomic", "proteomics",
  "rna", "snp", "snps", "transcriptomic", "transcriptomics", "variant", "linkage", "biomarkers")
# removed "cell" because was capturing e.g. {tidyr} and {googlesheets4}
# removed "sequencing" because was capturing e.g. {seriation}
# removed "evolution" because was capturing e.g. {BayesianTools}
# removed "phylogenetic", "phylogenetics" because of the Phylogenetics task view


# browse CRAN
cran_db <- tools::CRAN_package_db()

cran_tbl <- tibble::as_tibble(cran_db)#[-65])

cran_tbl_short <- cran_tbl %>% select(Package, Description, Published)

####
tidy_desc <- cran_tbl_short %>% unnest_tokens(word, Description)

data("stop_words")

cleaned_desc <- tidy_desc %>% anti_join(stop_words)

# Filtre sur les les packages non mis à jour depuis plus de 2 ou 3 ans 
omics_pkgs <- cleaned_desc %>% 
  group_by(Package) %>%
  filter(word %in% omics_words) %>% 
 # ungroup()
   filter(as.numeric(substr(Published, 1, 4)) >=2020) %>%
   ungroup()

n_distinct(omics_pkgs$Package) # 838

# Stemming step
omics_stem <- omics_pkgs %>%
  mutate(stem = wordStem(word)) 

omics_stem %>%
  count(stem, sort = TRUE)

liste_genet <- omics_stem %>% 
  filter(stem == "genet") 
n_distinct(liste_genet$Package) # 218

liste_others <- omics_pkgs %>% 
  filter(!(Package %in% liste_genet$Package)) # 620
n_distinct(liste_others$Package)

#write.table(unique(liste_genet$Package), file = "omics_genet_packages_2022-08-18.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(unique(liste_others$Package), file = "omics_others_packages_2022-08-25.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Avec info description du package
omics_others_desc <- cran_tbl_short[which(cran_tbl_short$Package %in% liste_others$Package),]
#write.table(omics_others_desc, file = "desc_omics_others_packages_2022-08-18.txt",
#            row.names = FALSE, col.names = FALSE, quote = FALSE)

