#This script uses the rep-seqs of the 18s dataset to assign taxonomy using the PR2 database and the SILVA database,
# It also saves the taxonomy results as RDS and CSV files for further analysis.

# Load Packages ---------------------------------
library(dada2)
library(Biostrings)
library(R.utils)


setwd(dir ="/scratch/tyjames_root/tyjames0/qmoon/SSBS/18S")


# Load FASTA and PR2 database ---------------------------------
fasta_file <- "dna-sequences_18s_SSBS.fasta"

unzip_ref_pr2 <- "../../AgeDiversityDistance/Moon18s/pr2_version_5.1.0_SSU_dada2.fasta"

# Read in the ASVs from the FASTA file
sequences_18s <- readDNAStringSet(fasta_file)

# Extract the headers (ASV identifiers) from the FASTA file
headers <- names(sequences_18s)

# PR2 database
pr2_db <- "pr2_version_5.1.0_SSU_dada2.fasta"
pr2_tax_levels <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")


# Assign PR2 taxonomy  ---------------------------------
taxa <- assignTaxonomy(sequences_18s,
                       unzip_ref_db,
                       multithread = TRUE,
                       minBoot = 80,           # Set the minimum bootstrap threshold
                       outputBootstraps = TRUE, # Return the bootstrap value  # Use the vector for taxonomic levels
                       taxLevels = tax_levels,  # Use the vector for taxonomic levels
                       verbose = TRUE)

# Convert the taxonomy result to a data frame
taxa_df <- as.data.frame(taxa)

# Add the headers (ASV identifiers) to the taxonomy data frame
taxa_df$ASV <- headers

# Move the row names (sequence identifiers) to a separate column
taxa_df$Sequence_ID <- rownames(taxa_df)


saveRDS(taxa_df, "taxa_df.rds")

# Optionally, save the results to a CSV
write.csv(taxa_df, "18s_taxonomic_pr2_bootstrap_SSBS.csv", row.names = FALSE)




# SILVA database ---------------------------------

silva_db <- "../../AgeDiversityDistance/Moon18s/SILVA_SSUfungi_nr99_v138_2_toGenus_trainset.fasta"



## --- Assign taxonomy with SILVA ---
taxa_silva <- assignTaxonomy(
  sequences_18s,
  silva_db,
  multithread = TRUE,
  minBoot = 80,
  outputBootstraps = TRUE,
  verbose = TRUE
)
taxa_silva_df <- as.data.frame(taxa_silva)
taxa_silva_df$ASV <- headers
taxa_silva_df$Sequence_ID <- rownames(taxa_silva_df)
saveRDS(taxa_silva_df, "taxa_df_silva.rds")
write.csv(taxa_silva_df, "18s_taxonomic_silva_bootstrap_SSBS.csv", row.names = FALSE)
