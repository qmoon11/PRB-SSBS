# Load Packages ---------------------------------
library(dplyr)
library(ggplot2)
library(vegan)
library(metacoder)
library(phyloseq)
library(ggVennDiagram)
library(patchwork)
library(ggplotify)
library(tidyverse)
library(RColorBrewer)
library(Polychrome)
library(stringr)
library(readr)
# Read in OTU tables, taxonomy ---------------------------------
#ITS
ITS_otu_table <- read.delim("ITS/ITS2_OTU_table_SSBS.tsv", check.names = FALSE)
ITS_taxonomy <- read.csv("ITS/final_ITS_taxonomy.csv")

#update column 1 name to be OTU_ID
colnames(ITS_taxonomy)[1] <- "OTU_ID"

#drop the notes column
ITS_taxonomy <- ITS_taxonomy %>%
  select(-notes)


#18S
SSBS_18s_otu_table <- read.delim("18S/18s_OTU_table.tsv", check.names = FALSE)
SSBS_18s_taxonomy <- read.csv("18S/combined_taxonomy_18s_SSBS.csv")

#delete unneeded columns
SSBS_18s_taxonomy <- SSBS_18s_taxonomy %>%
  select(-Sequence)

SSBS_18s_taxonomy <- SSBS_18s_taxonomy %>%
  select(-Phylum.1)

#16S
SSBS_16s_otu_table <- read.delim("16S/ASV_table_16S_SSBS.tsv", check.names = FALSE)
SSBS_16s_taxonomy <- read.delim("16S/taxonomy_16S_SSBS.tsv", check.names = FALSE)

#drop the confidence column
SSBS_16s_taxonomy <- SSBS_16s_taxonomy %>%
  select(-Confidence)

# Parse 16S Taxonomy ---------------------------------
# split QIIME tax string into ranks
tax_mat_16s <- SSBS_16s_taxonomy %>%
  separate(
    Taxon,
    into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
    sep = ";\\s*",
    fill = "right",
    extra = "drop"
  ) %>%
  mutate(across(-OTU_ID, ~ sub("^[a-z]__", "", .x))) %>%  # remove k__/p__/...
  column_to_rownames("OTU_ID") %>%
  as.matrix()

#convert to phyloseq tax_table
TAX_16s <- phyloseq::tax_table(tax_mat_16s)






# Build VEGAN friendly format ---------------------------------
#ITS
# community matrix for vegan: rows = samples, cols = OTUs
ITS_vegan <- ITS_otu_table %>%
  tibble::column_to_rownames("OUT_ID") %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
  as.matrix() %>%
  t()

ITS_vegan <- rbind(ITS_vegan, Pooled = colSums(ITS_vegan))



#16S
# community matrix for vegan: rows = samples, cols = OTUs
vegan_16s <- SSBS_16s_otu_table %>%
  tibble::column_to_rownames("OTU_ID") %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
  as.matrix() %>%
  t()

vegan_16s <- rbind(vegan_16s, Pooled = colSums(vegan_16s))



#18S
# community matrix for vegan: rows = samples, cols = OTUs
vegan_18s <- SSBS_18s_otu_table %>%
  tibble::column_to_rownames("OUT_ID") %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
  as.matrix() %>%
  t()

vegan_18s <- rbind(vegan_18s, Pooled = colSums(vegan_18s))



# Contaminant Check ---------------------------------
#18s look for animals, land plants. I inspected manually and found one OTU that is likely a contaminant: b124aef652447037e2e2ca37f44dc248, which is classified as Metazoa in PR2.
#After BLASTing the sequence, it is 100% match to many different mammals, so definitly a contaminant (but unclear exactly what it is). It is also very low abundance, so I will remove it from the dataset.
# % of total reads contributed by OTU 6
otuMetazoa_pct <- 100 * vegan_18s["Pooled", "b124aef652447037e2e2ca37f44dc248"] / sum(vegan_18s["Pooled", ])
#otuMetazoa_pct is 0.09%
otualg_pct <- 100 * vegan_18s["Pooled", "f858b11b4a48fee92b30a46f05828d6e"] / sum(vegan_18s["Pooled", ])

#drop observed Metazoa OTU
vegan_18s <- vegan_18s[, colnames(vegan_18s) != "b124aef652447037e2e2ca37f44dc248", drop = FALSE]
# sum(vegan_18s["Pooled", ])  # total reads after removing contaminant OTU. 24252 compared to 24,274 before. 

#ITS
#I will look for obligate human associated taxa using ClassifyITS (see BLAST script)
#no contaminants found in ITS dataset, so no need to remove any OTUs
#find any OTUs that are not kingdom=fungi in ITS_taxonomy. Green algae(600dfeb158b1481fb4a1de17e0b5e4f4, f8b025f755232b59ce31914e56591686) and plant (a16d03d3f78958a5f3b301a7b17709d1)
plant_ids <- c(
  "600dfeb158b1481fb4a1de17e0b5e4f4",
  "f8b025f755232b59ce31914e56591686",
  "a16d03d3f78958a5f3b301a7b17709d1"
)

# total % of pooled reads from these OTUs
plant_pct_total <- 100 * sum(ITS_vegan["Pooled", plant_ids]) / sum(ITS_vegan["Pooled", ])
plant_pct_total

plant_pct_each <- 100 * ITS_vegan["Pooled", plant_ids] / sum(ITS_vegan["Pooled", ])
plant_pct_each

#remove plant IDs from dataset
ITS_vegan <- ITS_vegan[, !colnames(ITS_vegan) %in% plant_ids, drop = FALSE]

sum(ITS_vegan["Pooled", ]) #4733 reads


#16S (look for contaminant genera list from Ruff et al 2024)
tax_df_16s <- as.data.frame(tax_mat_16s, stringsAsFactors = FALSE)

contam_16S <- read_lines("16S/contam_list_16S") |>
  str_trim() |>
  (\(x) x[x != "" & !str_starts(x, "#")])()

genus_col <- "Genus"

hits_exact <- tax_df_16s |>
  filter(tolower(.data[[genus_col]]) %in% tolower(contam_16S))

View(hits_exact)

#Erysipelothrix, Rhodococcus,and Romboutsia found. However I will leave them in as literature search shows all can be found in soil and water.



# Observed vs Chao1 richness (after remove contaminants)---------------------------------
#18S
t(estimateR(vegan_18s)[c("S.obs","S.chao1","se.chao1"), ])
#results 18S: Observed richness=34, Chao1 richness=34, SE=0

#16S
t(estimateR(vegan_16s)[c("S.obs","S.chao1","se.chao1"), ])
#results 16s: Observed richness=3994, Chao1 richness=4050, SE=14.04

#ITS
t(estimateR(ITS_vegan)[c("S.obs","S.chao1","se.chao1"), ])
#results ITS: Observed richness=42, Chao1 richness=42, SE=0


# Venn Diagram of replicates, build function---------------------------------
venn2_otu <- function(df,
                      col1 = NULL, col2 = NULL,
                      id_col = "OTU_ID",
                      set_names = c("Replicate 1", "Replicate 2"),
                      rotate_side_by_side = TRUE,
                      outline_color = "black",
                      outline_size = 0.7,
                      set_label_size = 4,
                      region_label_size = 3.5) {
  stopifnot(length(set_names) == 2)
  
  # pick columns
  if (is.null(col1) || is.null(col2)) {
    # default: first two non-id columns
    cols <- setdiff(names(df), id_col)
    stopifnot(length(cols) >= 2)
    col1 <- cols[1]; col2 <- cols[2]
  }
  
  # ids
  ids <- if (id_col %in% names(df)) df[[id_col]] else rownames(df)
  if (is.null(ids)) stop("No OTU IDs found: provide id_col or rownames(df).")
  
  # numeric counts + NA->0
  c1 <- as.numeric(df[[col1]]); c1[is.na(c1)] <- 0
  c2 <- as.numeric(df[[col2]]); c2[is.na(c2)] <- 0
  
  keep <- (c1 + c2) > 0
  ids <- ids[keep]; c1 <- c1[keep]; c2 <- c2[keep]
  
  # sets
  sets <- setNames(list(ids[c1 > 0], ids[c2 > 0]), set_names)
  
  # partitions
  only1  <- c1 > 0 & c2 == 0
  shared <- c1 > 0 & c2 > 0
  only2  <- c2 > 0 & c1 == 0
  
  reads_only1  <- sum(c1[only1])
  reads_shared <- sum(c1[shared]) + sum(c2[shared])
  reads_only2  <- sum(c2[only2])
  total_reads  <- reads_only1 + reads_shared + reads_only2
  
  lab_only1  <- sprintf("%d ASVs\n(%.1f%% pooled reads)", sum(only1),  100 * reads_only1  / total_reads)
  lab_shared <- sprintf("%d ASVs\n(%.1f%% pooled reads)", sum(shared), 100 * reads_shared / total_reads)
  lab_only2  <- sprintf("%d ASVs\n(%.1f%% pooled reads)", sum(only2),  100 * reads_only2  / total_reads)
  
  # ggVennDiagram plot data
  venn <- ggVennDiagram::Venn(sets)
  pd <- ggVennDiagram::process_data(venn)
  
  if (rotate_side_by_side) {
    swap_xy <- function(d) { tmp <- d$X; d$X <- d$Y; d$Y <- tmp; d }
    pd$setEdge     <- swap_xy(pd$setEdge)
    pd$setLabel    <- swap_xy(pd$setLabel)
    pd$regionEdge  <- swap_xy(pd$regionEdge)
    pd$regionLabel <- swap_xy(pd$regionLabel)
  }
  
  shared_name <- paste(sort(set_names), collapse = "/")
  
  region_df <- dplyr::as_tibble(pd$regionLabel) %>%
    dplyr::mutate(label = dplyr::case_when(
      name == set_names[1] ~ lab_only1,
      name == set_names[2] ~ lab_only2,
      name == shared_name  ~ lab_shared,
      TRUE ~ name
    ))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = pd$setEdge,
                       ggplot2::aes(X, Y, group = id),
                       color = outline_color, linewidth = outline_size) +
    ggplot2::geom_text(data = pd$setLabel,
                       ggplot2::aes(X, Y, label = name),
                       size = set_label_size) +
    ggplot2::geom_text(data = region_df,
                       ggplot2::aes(X, Y, label = label),
                       size = region_label_size) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  
  summary <- data.frame(
    region = c(paste0(set_names[1], "_only"), "shared", paste0(set_names[2], "_only")),
    n_asvs = c(sum(only1), sum(shared), sum(only2)),
    reads = c(reads_only1, reads_shared, reads_only2),
    pct_pooled_reads = c(100 * reads_only1/total_reads,
                         100 * reads_shared/total_reads,
                         100 * reads_only2/total_reads)
  )
  
  list(plot = p, summary = summary, sets = sets)
}

# Venn 16S ---------------------------------
venn_16s <- venn2_otu(SSBS_16s_otu_table,
                 col1 = "Shadmani_1_16s", col2 = "Shadmani_2_16s",
                 id_col = "OTU_ID",
                 set_names = c("Replicate 1", "Replicate 2"))
venn_16s$plot


# Venn 18S ---------------------------------
#drop otub124aef652447037e2e2ca37f44dc248 from SSBS_18s_otu_table

SSBS_18s_otu_table <- SSBS_18s_otu_table |>
  dplyr::filter(OUT_ID != "b124aef652447037e2e2ca37f44dc248")

venn_18s <- venn2_otu(SSBS_18s_otu_table,
                      col1 = "Shadmani_1_18S", col2 = "Shadmani_2_18S",
                      id_col = "OTU_ID",
                      set_names = c("Replicate 1", "Replicate 2"))
venn_18s$plot


# Venn ITS ---------------------------------
venn_ITS <- venn2_otu(ITS_otu_table,
                      col1 = "Shadmani_1_ITS", col2 = "Shadmani_2_ITS",
                      id_col = "OTU_ID",
                      set_names = c("Replicate 1", "Replicate 2"))
venn_ITS$plot



# Rarefaction curve function ---------------------------------
rarecurve_2rep_pooled <- function(df,
                                  col1, col2,
                                  id_col = "OTU_ID",
                                  rep_names = c("Replicate 1", "Replicate 2"),
                                  step = 100,
                                  cols = c("#0073C2FF", "red", "black"),
                                  lwd = 2,
                                  xlab = "Reads",
                                  ylab = "ASVs",
                                  legend_pos = "bottomright",
                                  label = FALSE,
                                  ...) {
  stopifnot(length(rep_names) == 2, length(cols) == 3)
  
  ids <- if (id_col %in% names(df)) df[[id_col]] else rownames(df)
  if (is.null(ids)) stop("No OTU IDs found: provide id_col or rownames(df).")
  
  mat <- df[, c(col1, col2), drop = FALSE]
  rownames(mat) <- ids
  mat[] <- lapply(mat, as.numeric)
  mat[is.na(mat)] <- 0
  
  comm <- t(as.matrix(mat))              # rows = replicates, cols = OTUs
  rownames(comm) <- rep_names
  comm <- rbind(comm, Pooled = comm[rep_names[1], ] + comm[rep_names[2], ])
  
  vegan::rarecurve(
    comm,
    step = step,
    sample = min(rowSums(comm)),
    col = cols,
    lwd = lwd,
    xlab = xlab,
    ylab = ylab,
    label = label,
    ...
  )
  
  legend(legend_pos, legend = rownames(comm), col = cols, lwd = lwd, bty = "n")
  
  invisible(comm)  # returns the vegan community matrix invisibly
}
# Rarefaction curve 16s---------------------------------
rarecurve_plot_16s <- as.ggplot(function() {
  rarecurve_2rep_pooled(
    SSBS_16s_otu_table,
    col1 = "Shadmani_1_16s",
    col2 = "Shadmani_2_16s",
    id_col = "OTU_ID",
    ylab = "16S ASVs"
  )
})

# Rarefaction curve 18S ---------------------------------
rarecurve_plot_18s <- as.ggplot(function() {
  rarecurve_2rep_pooled(
    SSBS_18s_otu_table,
    col1 = "Shadmani_1_18S",
    col2 = "Shadmani_2_18S",
    id_col = "OTU_ID",
    ylab = "18S OTUs"
  )
})


# Rarefaction curve ITS ---------------------------------
rarecurve_plot_ITS <- as.ggplot(function() {
  rarecurve_2rep_pooled(
    ITS_otu_table,,
    col1 = "Shadmani_1_ITS",
    col2 = "Shadmani_2_ITS",
    id_col = "OTU_ID",
    step = 10,
    ylab = "ITS OTUs"
  )
})

rarecurve_plot_ITS


# Combine venn and rarefication plots  ---------------------------------
combined_venn_rare <- ((venn_16s$plot | rarecurve_plot_16s) /
                         (venn_18s$plot | rarecurve_plot_18s) /
                         (venn_ITS$plot | rarecurve_plot_ITS)) +
  plot_annotation(tag_levels = "A")

# PDF
ggsave("figures/combined_venn_rare.pdf", combined_venn_rare, width = 12, height = 14)

# JPG
ggsave("figures/combined_venn_rare.jpg", combined_venn_rare, width = 12, height = 14, dpi = 300)

# Load in phyloseq ---------------------------------
#16S
#already made TAX_16s above

# OTU table -> matrix (taxa/OTUs as rows)
otu_mat_16s <- SSBS_16s_otu_table %>%
  column_to_rownames("OTU_ID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# add pooled sample column
otu_mat_16s <- cbind(otu_mat_16s, Pooled = rowSums(otu_mat_16s))
#Chnage Shadmani_1_16s and Shadmani_2_16s to Replicate 1 and Replicate 2
colnames(otu_mat_16s) <- gsub("Shadmani_1_16s", "Replicate 1", colnames(otu_mat_16s))
colnames(otu_mat_16s) <- gsub("Shadmani_2_16s", "Replicate 2", colnames(otu_mat_16s))

# taxonomy matrix (already made with rownames = OTU_ID)
ps_16s <- phyloseq(
  otu_table(otu_mat_16s, taxa_are_rows = TRUE),
  tax_table(TAX_16s)
)

#18S
# OTU table -> matrix (taxa/OTUs as rows)
otu_mat_18s <- SSBS_18s_otu_table %>%
  column_to_rownames("OUT_ID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
# add pooled sample column
otu_mat_18s <- cbind(otu_mat_18s, Pooled = rowSums(otu_mat_18s))
#Change Shadmani_1_18S and Shadmani_2_18S to Replicate 1 and Replicate 2
colnames(otu_mat_18s) <- gsub("Shadmani_1_18S", "Replicate 1", colnames(otu_mat_18s))
colnames(otu_mat_18s) <- gsub("Shadmani_2_18S", "Replicate 2", colnames(otu_mat_18s))


## Taxonomy table -> matrix with rownames matching OTU IDs
tax_mat_18s <- SSBS_18s_taxonomy %>%
  column_to_rownames("OTU_ID") %>%
  mutate(across(everything(), as.character)) %>%
  as.matrix()

# keep only OTUs present in the OTU table, and in the same order
tax_mat_18s <- tax_mat_18s[rownames(otu_mat_18s), , drop = FALSE]

## Build phyloseq object
ps_18s <- phyloseq(
  otu_table(otu_mat_18s, taxa_are_rows = TRUE),
  tax_table(tax_mat_18s)
)


# ITS
# remove plant from OTU table (ID column is OUT_ID)
ITS_otu_table <- ITS_otu_table %>%
  dplyr::filter(!OUT_ID %in% plant_ids)

# remove from taxonomy (ID column is OTU_ID)
ITS_taxonomy <- ITS_taxonomy %>%
  dplyr::filter(!OTU_ID %in% plant_ids)

# OTU table -> matrix (taxa/OTUs as rows)
otu_mat_ITS <- ITS_otu_table %>%
  column_to_rownames("OUT_ID") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
# add pooled sample column
otu_mat_ITS <- cbind(otu_mat_ITS, Pooled = rowSums(otu_mat_ITS))
#Change Shadmani_1_ITS and Shadmani_2_ITS to Replicate 1 and Replicate 2
colnames(otu_mat_ITS) <- gsub("Shadmani_1_ITS", "Replicate 1", colnames(otu_mat_ITS))
colnames(otu_mat_ITS) <- gsub("Shadmani_2_ITS", "Replicate 2", colnames(otu_mat_ITS)) 

## Taxonomy table -> matrix with rownames matching OTU IDs
# assumes your taxonomy table has an OTU ID column called OUT_ID (adjust if different)
tax_mat_ITS <- ITS_taxonomy %>%
  column_to_rownames("OTU_ID") %>%
  mutate(across(everything(), as.character)) %>%
  as.matrix()


## Build phyloseq object
ps_ITS <- phyloseq(
  otu_table(otu_mat_ITS, taxa_are_rows = TRUE),
  tax_table(tax_mat_ITS)
)




# Stacked Bar Function ---------------------------------
plot_rank_bar_other <- function(ps, rank, threshold = NULL, xlab = NULL,
                                legend_title_size = 14,
                                legend_text_size  = 12,
                                y_axis_title_size = 14,
                                y_axis_text_size  = 12,
                                legend_ncol = 1,
                                ylab = "Relative Abundance") {
  ps_r <- phyloseq::tax_glom(ps, taxrank = rank, NArm = FALSE)
  ps_r <- phyloseq::transform_sample_counts(ps_r, function(x) x / sum(x))
  df <- phyloseq::psmelt(ps_r)
  
  df[[rank]] <- as.character(df[[rank]])
  df[[rank]][is.na(df[[rank]]) | df[[rank]] == ""] <- "Unclassified"
  
  # only lump to "Other" if threshold is provided
  if (!is.null(threshold)) {
    df[[rank]][df$Abundance < threshold] <- "Other"
    df[[rank]] <- factor(df[[rank]], levels = c(setdiff(unique(df[[rank]]), "Other"), "Other"))
  } else {
    df[[rank]] <- factor(df[[rank]])
  }
  
  lvls <- levels(df[[rank]])
  non_other <- setdiff(lvls, "Other")
  cols_non_other <- Polychrome::palette36.colors(length(non_other))
  pal <- c(setNames(cols_non_other, non_other), Other = "black")
  
  ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +  # keep all legend levels
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol)) +
    ggplot2::labs(x = xlab, y = ylab, fill = rank) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.title.y = ggplot2::element_text(size = y_axis_title_size, color = "black"),
      axis.text.y  = ggplot2::element_text(size = y_axis_text_size, color = "black"),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text  = ggplot2::element_text(size = legend_text_size)
    )
}

# Stacked Bar 16s ---------------------------------
#all three
p_phylum_16s <- plot_rank_bar_other(
  ps_16s,
  rank = "Phylum",
  threshold = 0.01,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 2
)


p_class_16s <- plot_rank_bar_other(
  ps_16s,
  rank = "Class",
  threshold = 0.01,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 2
)


#combie p_phylum_16s and p_class_16s into one plot
combined_16s <- combined_16s & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_16s
ggsave("figures/combined_16s_phylum_class.pdf", combined_16s, width = 10, height = 12)
#save as jpg
ggsave("figures/combined_16s_phylum_class.jpg", combined_16s, width = 10, height = 12, dpi = 300)

#pooled
ps_16s_pooled <- prune_samples(sample_names(ps_16s) == "Pooled", ps_16s)

p_class_16s_pooled <- plot_rank_bar_other(
  ps_16s_pooled,
  rank = "Class",
  threshold = 0.01,
  xlab = NULL,
  ylab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),  # remove tick labels
    axis.ticks.x = ggplot2::element_blank(),  # remove tick marks
    axis.title.x = ggplot2::element_blank()   # remove x-axis title
  )

p_phylum_16s_pooled <- plot_rank_bar_other(
  ps_16s_pooled,
  rank = "Phylum",
  threshold = 0.01,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),  # remove tick labels
    axis.ticks.x = ggplot2::element_blank(),  # remove tick marks
    axis.title.x = ggplot2::element_blank()   # remove x-axis title
  )


# Stacked Bar 18S ---------------------------------
#all three

# Use Kingdom.Subdivision if it's fungi; otherwise use Division
tax_df <- as.data.frame(tax_table(ps_18s), stringsAsFactors = FALSE)

tax_df$Kingdom_or_Division <- ifelse(
  !is.na(tax_df$Kingdom.Subdivision) &
    grepl("Fungi", tax_df$Kingdom.Subdivision, ignore.case = TRUE),
  paste0(tax_df$Kingdom.Subdivision, "*"),  # fungi at Kingdom.Subdivision -> keep it (+ *)
  tax_df$Division                           # otherwise -> Division
)

tax_df$Kingdom_or_Division[
  is.na(tax_df$Kingdom_or_Division) | tax_df$Kingdom_or_Division == ""
] <- "Unclassified"

tax_table(ps_18s) <- tax_table(as.matrix(tax_df))

p_kingdom_18s <- plot_rank_bar_other(
  ps_18s,
  rank = "Kingdom_or_Division",
  threshold = NULL,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::labs(fill = "Kingdom (fungi)\nDivision (non-fungi)")


p_kingdom_18s



tax_df <- as.data.frame(tax_table(ps_18s), stringsAsFactors = FALSE)

# Class for all; add * if fungi at Kingdom.Subdivision
tax_df$Class_star <- tax_df$Class
tax_df$Class_star[
  !is.na(tax_df$Kingdom.Subdivision) &
    grepl("Fungi", tax_df$Kingdom.Subdivision, ignore.case = TRUE) &
    !is.na(tax_df$Class) & tax_df$Class != ""
] <- paste0(tax_df$Class[
  !is.na(tax_df$Kingdom.Subdivision) &
    grepl("Fungi", tax_df$Kingdom.Subdivision, ignore.case = TRUE) &
    !is.na(tax_df$Class) & tax_df$Class != ""
], "*")

tax_df$Class_star[is.na(tax_df$Class_star) | tax_df$Class_star == ""] <- "Unclassified"

tax_table(ps_18s) <- tax_table(as.matrix(tax_df))

p_class_18s <- plot_rank_bar_other(
  ps_18s,
  rank = "Class_star",
  threshold = NULL,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::labs(fill = "Class")

p_class_18s

#combie p_phylum_18s and p_class_18s into one plot
combined_18s <- (p_kingdom_18s / p_class_18s) + plot_annotation(tag_levels = "A")

combined_18s <- combined_18s & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_18s
ggsave("figures/combined_18s_phylum_class.pdf", combined_18s, width = 10, height = 12)
#save as jpg
ggsave("figures/combined_18s_phylum_class.jpg", combined_18s, width = 10, height = 12, dpi = 300)

#pooled
## --- keep the same star/grouping logic on ps_18s ---
tax_df <- as.data.frame(tax_table(ps_18s), stringsAsFactors = FALSE)

# Kingdom.Subdivision* for fungi, Division for non-fungi
tax_df$Kingdom_or_Division <- ifelse(
  !is.na(tax_df$Kingdom.Subdivision) &
    grepl("Fungi", tax_df$Kingdom.Subdivision, ignore.case = TRUE),
  paste0(tax_df$Kingdom.Subdivision, "*"),
  tax_df$Division
)
tax_df$Kingdom_or_Division[is.na(tax_df$Kingdom_or_Division) | tax_df$Kingdom_or_Division == ""] <- "Unclassified"

# Class for all, but add * if fungi at Kingdom.Subdivision
tax_df$Class_star <- tax_df$Class
is_fungi <- !is.na(tax_df$Kingdom.Subdivision) & grepl("Fungi", tax_df$Kingdom.Subdivision, ignore.case = TRUE)
has_class <- !is.na(tax_df$Class) & tax_df$Class != ""
tax_df$Class_star[is_fungi & has_class] <- paste0(tax_df$Class[is_fungi & has_class], "*")
tax_df$Class_star[is.na(tax_df$Class_star) | tax_df$Class_star == ""] <- "Unclassified"

tax_table(ps_18s) <- tax_table(as.matrix(tax_df))

## --- pooled plots (18S) ---
ps_18s_pooled <- prune_samples(sample_names(ps_18s) == "Pooled", ps_18s)

p_kingdom_18s_pooled <- plot_rank_bar_other(
  ps_18s_pooled,
  rank = "Kingdom_or_Division",
  threshold = NULL,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::labs(fill = "Kingdom (fungi)\nDivision (non-fungi)") +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

p_class_18s_pooled <- plot_rank_bar_other(
  ps_18s_pooled,
  rank = "Class_star",
  threshold = NULL,
  xlab = NULL,
  ylab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::labs(fill = "Class") +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

p_kingdom_18s_pooled
p_class_18s_pooled

# Stacked Bar ITS ---------------------------------
#all three
p_phylum_ITS <- plot_rank_bar_other(
  ps_ITS,
  rank = "phylum",
  threshold = NULL,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
)

p_phylum_ITS

p_class_ITS <- plot_rank_bar_other(
  ps_ITS,
  rank = "class",
  threshold = NULL,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
)
p_class_ITS 


#combie p_phylum_ITS and p_class_ITS into one plot
combined_ITS <- (p_phylum_ITS / p_class_ITS) + plot_annotation(tag_levels = "A")

combined_ITS <- combined_ITS & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_ITS
ggsave("figures/combined_ITS_phylum_class.pdf", combined_ITS, width = 10, height = 12)
#save as jpg
ggsave("figures/combined_ITS_phylum_class.jpg", combined_ITS, width = 10, height = 12, dpi = 300)





#pooled

ps_ITS_pooled <- prune_samples(sample_names(ps_ITS) == "Pooled", ps_ITS)


p_class_ITS_pooled <- plot_rank_bar_other(
  ps_ITS_pooled,
  rank = "class",
  threshold = NULL,
  xlab = NULL,
  ylab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::labs(fill = "Class") +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),  # remove tick labels
    axis.ticks.x = ggplot2::element_blank(),  # remove tick marks
    axis.title.x = ggplot2::element_blank()   # remove x-axis title
  )

p_phylum_ITS_pooled <- plot_rank_bar_other(
  ps_ITS_pooled,
  rank = "phylum",
  threshold = NULL,
  xlab = NULL,
  legend_title_size = 16,
  legend_text_size  = 16,
  y_axis_title_size = 16,
  y_axis_text_size  = 16,
  legend_ncol = 1
) +
  ggplot2::labs(fill = "Phylum") +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank()
  )

p_class_ITS_pooled
p_phylum_ITS_pooled





# Combine stacked bars pooled samples ---------------------------------
combined_pooled <- (p_phylum_16s_pooled | p_class_16s_pooled) /
  (p_kingdom_18s_pooled | p_class_18s_pooled) /
  (p_phylum_ITS_pooled | p_class_ITS_pooled) +
  plot_annotation(tag_levels = "A")
combined_pooled <- combined_pooled & theme(plot.tag = element_text(size = 18, face = "bold"))

combined_pooled

ggsave("figures/combined_pooled_phylum_class.pdf", combined_pooled, width = 12, height = 14)
ggsave("figures/combined_pooled_phylum_class.jpg", combined_pooled, width = 12, height = 14, dpi = 300)

# Metacoder ---------------------------------



# Rank Abundance ---------------------------------

otu_table <- ITS_vegan

highlight_ids <- c(
  "f4ca0e638d89d15fdf2efa177740c9ed",
  "9d66dd8f2085b9f1745d039dbfb65b30",
  "feafda5ba95fcf36d141ed4058f141fc"
)

highlight_labels <- c(
  "f4ca0e638d89d15fdf2efa177740c9ed" = "Scedosporium",
  "9d66dd8f2085b9f1745d039dbfb65b30" = "Exophiala",
  "feafda5ba95fcf36d141ed4058f141fc" = "Paracremonium"
)

pooled_counts <- colSums(otu_table, na.rm = TRUE)
pooled_rel <- pooled_counts / sum(pooled_counts)

rank_abundance <- data.frame(
  otu_id = names(pooled_rel),
  rel_abundance = as.numeric(pooled_rel),
  stringsAsFactors = FALSE
)

rank_abundance <- rank_abundance[rank_abundance$rel_abundance > 0, , drop = FALSE]
rank_abundance <- rank_abundance[order(rank_abundance$rel_abundance, decreasing = TRUE), , drop = FALSE]
rank_abundance$rank <- seq_len(nrow(rank_abundance))

rank_abundance$highlight <- rank_abundance$otu_id %in% highlight_ids
rank_abundance$taxon_label <- unname(highlight_labels[rank_abundance$otu_id])

rank_abundance <- ggplot(rank_abundance, aes(rank, rel_abundance)) +
  geom_line(color = "black", linewidth = 0.3) +
  geom_point(
    aes(fill = highlight),
    shape = 21, color = "black", stroke = 0.8, size = 2.3
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(`TRUE` = "black", `FALSE` = "white"),
    breaks = c(TRUE, FALSE),
    labels = c("Isolated", "Not isolated")
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "OTU rank (pooled)", y = "Relative abundance (pooled)") +
  theme_classic() +
  geom_text(
    data = subset(rank_abundance, highlight),
    aes(label = taxon_label),
    nudge_y = 0.01,
    check_overlap = TRUE,
    size = 3
  )

#ggsave rank abundance plot
ggsave("figures/rank_abundance_ITS.pdf", rank_abundance, width = 6, height = 4)
ggsave("figures/rank_abundance_ITS.jpg", rank_abundance, width = 6, height = 4, dpi = 300)
