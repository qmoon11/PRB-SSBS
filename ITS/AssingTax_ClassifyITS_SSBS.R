#This script uses the github package ClassifyITS to assign taxonomy to the ITS dataset using the UNITE database. 
#The ClassifyITS version (0.01) used is identical to the version currently under review at CRAN.


# Load Packages ---------------------------------
# Install ClassifyITS from GitHub
devtools::install_github("qmoon11/ClassifyITS")

# Load the package
library(ClassifyITS)


# Run ClassifyITS ---------------------------------

ITS_taxonomy <- ITS_assignment(
  blast_file = "megablast_ITS2.tsv",            # Path to BLAST .TSV file
  rep_fasta = "dna-sequences_ITS_SSBS.fasta"      # Path to FASTA file containg representative sequences used to generate the BLAST results
)

#Output csv for manual inspection
