#Generate KEGG gene-to-KO files

# ============================================================
# CREATE TWO-COLUMN KEGG MAPPER INPUT FILES
#
# Input files:
#   query.emapper.annotations_Acr
#   query.emapper.annotations_Exo
#   query.emapper.annotations_Pen
#
# Output format:
#
#   gene1    K02874
#   gene4    K00416
#
# The separator is a tab. Output files have no header.
#
# Actual protein identifiers from the eggNOG query column are
# retained. If one protein has multiple KO assignments, each
# assignment is written on a separate row.
# ============================================================


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

annotation_files <- c(
  Acr = "query.emapper.annotations_Acr",
  Exo = "query.emapper.annotations_Exo",
  Pen = "query.emapper.annotations_Pen",
  Para = "query.emapper.annotations_Para"
)


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

output_directory <- "KEGG_Mapper_inputs"

dir.create(
  output_directory,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# SAFETY CHECKS
# ============================================================

missing_files <- annotation_files[
  !file.exists(annotation_files)
]

if (length(missing_files) > 0) {
  stop(
    "The following files were not found:\n",
    paste(
      names(missing_files),
      missing_files,
      sep = ": ",
      collapse = "\n"
    )
  )
}


# ============================================================
# READ ONE EGGNOG-MAPPER ANNOTATION FILE
# ============================================================

read_emapper_file <- function(file_path) {
  
  lines <- readLines(
    file_path,
    warn = FALSE
  )
  
  header_index <- grep(
    "^#query\\t",
    lines
  )[1]
  
  if (is.na(header_index)) {
    stop(
      "Could not find the #query header in: ",
      file_path
    )
  }
  
  header_names <- strsplit(
    sub(
      "^#",
      "",
      lines[header_index]
    ),
    split = "\t",
    fixed = TRUE
  )[[1]]
  
  # Retain annotation rows after the header, excluding comments.
  data_lines <- lines[
    seq.int(
      header_index + 1,
      length(lines)
    )
  ]
  
  data_lines <- data_lines[
    nzchar(trimws(data_lines)) &
      !grepl("^#", data_lines)
  ]
  
  if (length(data_lines) == 0) {
    stop(
      "No annotation rows were found in: ",
      file_path
    )
  }
  
  temporary_file <- tempfile(
    fileext = ".tsv"
  )
  
  on.exit(
    unlink(temporary_file),
    add = TRUE
  )
  
  writeLines(
    data_lines,
    temporary_file
  )
  
  annotations <- read.delim(
    temporary_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (ncol(annotations) < length(header_names)) {
    stop(
      "The annotation rows contain fewer columns than the header in: ",
      file_path
    )
  }
  
  if (ncol(annotations) > length(header_names)) {
    extra_names <- paste0(
      "extra_column_",
      seq_len(
        ncol(annotations) - length(header_names)
      )
    )
    
    names(annotations) <- c(
      header_names,
      extra_names
    )
  } else {
    names(annotations) <- header_names
  }
  
  annotations
}


# ============================================================
# EXTRACT GENE-TO-KO ASSIGNMENTS
# ============================================================

extract_gene_to_ko <- function(
    file_path,
    genome_name
) {
  
  annotations <- read_emapper_file(
    file_path
  )
  
  if (!"query" %in% names(annotations)) {
    stop(
      "No query column was found in: ",
      file_path
    )
  }
  
  if (!"KEGG_ko" %in% names(annotations)) {
    stop(
      "No KEGG_ko column was found in: ",
      file_path,
      "\nAvailable columns:\n",
      paste(
        names(annotations),
        collapse = ", "
      )
    )
  }
  
  gene_to_ko_rows <- vector(
    "list",
    nrow(annotations)
  )
  
  for (row_index in seq_len(nrow(annotations))) {
    
    gene_id <- as.character(
      annotations$query[row_index]
    )
    
    ko_text <- as.character(
      annotations$KEGG_ko[row_index]
    )
    
    if (
      is.na(ko_text) ||
      ko_text == "" ||
      ko_text == "-"
    ) {
      next
    }
    
    # Extract all valid K numbers from the annotation.
    ko_positions <- gregexpr(
      "K[0-9]{5}",
      ko_text,
      perl = TRUE
    )
    
    ko_ids <- regmatches(
      ko_text,
      ko_positions
    )[[1]]
    
    if (
      length(ko_ids) == 0 ||
      identical(ko_ids, character(0))
    ) {
      next
    }
    
    ko_ids <- sort(
      unique(ko_ids)
    )
    
    gene_to_ko_rows[[row_index]] <- data.frame(
      genome = genome_name,
      gene = gene_id,
      KO = ko_ids,
      stringsAsFactors = FALSE
    )
  }
  
  gene_to_ko_rows <- gene_to_ko_rows[
    !vapply(
      gene_to_ko_rows,
      is.null,
      logical(1)
    )
  ]
  
  if (length(gene_to_ko_rows) == 0) {
    warning(
      genome_name,
      ": no valid KEGG KO assignments were found."
    )
    
    return(
      data.frame(
        genome = character(),
        gene = character(),
        KO = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  gene_to_ko <- do.call(
    rbind,
    gene_to_ko_rows
  )
  
  gene_to_ko <- unique(
    gene_to_ko
  )
  
  gene_to_ko <- gene_to_ko[
    order(
      gene_to_ko$gene,
      gene_to_ko$KO
    ),
    ,
    drop = FALSE
  ]
  
  rownames(gene_to_ko) <- NULL
  
  gene_to_ko
}


# ============================================================
# PROCESS ALL GENOMES
# ============================================================

all_results <- list()
summary_results <- list()

for (genome_name in names(annotation_files)) {
  
  cat(
    "\nProcessing ",
    genome_name,
    "...\n",
    sep = ""
  )
  
  gene_to_ko <- extract_gene_to_ko(
    file_path = annotation_files[[genome_name]],
    genome_name = genome_name
  )
  
  # ----------------------------------------------------------
  # Write the required two-column file:
  #
  # gene1<TAB>K02874
  #
  # No header and no row names.
  # ----------------------------------------------------------
  
  two_column_file <- file.path(
    output_directory,
    paste0(
      genome_name,
      "_gene_to_KO_for_KEGG_Mapper.txt"
    )
  )
  
  write.table(
    gene_to_ko[
      ,
      c(
        "gene",
        "KO"
      ),
      drop = FALSE
    ],
    file = two_column_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    na = ""
  )
  
  # Also write a one-column KO file as an optional alternative.
  ko_only_file <- file.path(
    output_directory,
    paste0(
      genome_name,
      "_KOs_only.txt"
    )
  )
  
  writeLines(
    sort(
      unique(
        gene_to_ko$KO
      )
    ),
    ko_only_file
  )
  
  all_results[[genome_name]] <- gene_to_ko
  
  summary_results[[genome_name]] <- data.frame(
    genome = genome_name,
    annotation_file = annotation_files[[genome_name]],
    n_gene_KO_assignments = nrow(gene_to_ko),
    n_genes_with_KO = length(
      unique(
        gene_to_ko$gene
      )
    ),
    n_unique_KOs = length(
      unique(
        gene_to_ko$KO
      )
    ),
    output_file = two_column_file,
    stringsAsFactors = FALSE
  )
  
  cat(
    "  Gene–KO assignments: ",
    nrow(gene_to_ko),
    "\n",
    "  Genes with KO: ",
    length(unique(gene_to_ko$gene)),
    "\n",
    "  Unique KOs: ",
    length(unique(gene_to_ko$KO)),
    "\n",
    "  Output: ",
    two_column_file,
    "\n",
    sep = ""
  )
}


# ============================================================
# WRITE COMBINED TABLES
# ============================================================

all_gene_to_ko <- do.call(
  rbind,
  all_results
)

rownames(all_gene_to_ko) <- NULL

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(summary_table) <- NULL


write.table(
  all_gene_to_ko,
  file = file.path(
    output_directory,
    "all_genomes_gene_to_KO.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


write.table(
  summary_table,
  file = file.path(
    output_directory,
    "all_genomes_KO_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# ============================================================
# FINAL SUMMARY
# ============================================================

cat(
  "\n============================================================\n",
  "KEGG MAPPER FILES CREATED\n",
  "============================================================\n",
  sep = ""
)

print(
  summary_table,
  row.names = FALSE
)

cat(
  "\nUpload these two-column files to KEGG Mapper:\n"
)

for (genome_name in names(annotation_files)) {
  cat(
    "  ",
    file.path(
      output_directory,
      paste0(
        genome_name,
        "_gene_to_KO_for_KEGG_Mapper.txt"
      )
    ),
    "\n",
    sep = ""
  )
}

cat(
  "\nKEGG Mapper Reconstruct Pathway:\n",
  "https://www.kegg.jp/kegg/mapper/reconstruct.html\n",
  sep = ""
)

