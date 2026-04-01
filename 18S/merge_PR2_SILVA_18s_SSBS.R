#This script loads the .csv output from the 18s results using the silva fungi and pr2 databases.
#I kept the pr2 assingments for all non-fungal OTUs and the silva assingments for all fungal OTUs.
#Be sure to manually check any OTUs that differ in taxonomy between two approaches.
# The output is a final 18s taxonomy with cleaned column names and missing taxonomy filled with "Unclassified".

#### Load Packages ####
library(dplyr)
library(stringr)

#### Load in PR2 and Silva Taxonomy ####
# Read both CSVs
setwd("/Users/quinnmoon/Downloads/PRB_SSBS/18S")

pr2 <- read.csv("18s_taxonomic_pr2_bootstrap_SSBS.csv")
silva <- read.csv("18s_taxonomic_silva_bootstrap_SSBS.csv")


#### Combine and clean OTU Tables ####
# 1. Identify OTUs classified as Fungi in pr2
fungal_otus <- pr2 %>%
  filter("tax.Kingdom/Subdivision" == "Fungi") %>%
  pull(ASV)  # Extract OTU identifiers

# 2. Remove fungal OTUs from pr2
pr2_trimmed <- pr2 %>%
  filter(!ASV %in% fungal_otus)

# 3. Keep only fungal OTUs in silva
silva_trimmed <- silva %>%
  filter(ASV %in% fungal_otus)

# Drop columns containing "boot" (bootstrap confidence scores)
pr2_trimmed <- pr2_trimmed %>%
  select(-contains("boot"))

silva_trimmed <- silva_trimmed %>%
  select(-contains("boot"))


# Check if both data frames have the same column names in the same order
identical(names(pr2_trimmed), names(silva_trimmed))

# Combine non-fungal PR2 OTUs with fungal SILVA OTUs
combined_taxonomy_18s <- bind_rows(pr2_trimmed, silva_trimmed)

# Rename identifiers and sequence columns for clarity
combined_taxonomy_18s <- combined_taxonomy_18s %>%
  rename(
    OTU_ID = ASV,
    Sequence = Sequence_ID
  )



# Remove "tax." prefix from taxonomy columns
combined_taxonomy_18s <- combined_taxonomy_18s %>%
  rename_with(~ sub("^tax\\.", "", .), starts_with("tax."))


# Replace taxonomy labels containing underscores with "Unclassified"
combined_taxonomy_18s <- combined_taxonomy_18s %>%
  mutate(across(-"Kingdom.Subdivision", ~ ifelse(str_detect(.x, "_"), "Unclassified", .x)))

View(combined_taxonomy_18s)
# Export final combined OTU taxonomy table to CSV
write.csv(combined_taxonomy_18s, "combined_taxonomy_18s_SSBS.csv", row.names = FALSE)

#Manually inspect all OTUs unclasifed above Class.
