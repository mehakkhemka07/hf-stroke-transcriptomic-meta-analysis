library(DESeq2)

# ---------------------------------------------------------------
# Reusable DESeq2 runner
# ---------------------------------------------------------------
run_deseq2 <- function(counts_path, metadata_path, treatment_label, reference_label,
                        annot, output_path) {

  # Load count matrix
  counts <- read.table(counts_path, header = TRUE, sep = "\t", row.names = 1)

  # Load sample metadata
  sample_meta <- read.csv(metadata_path, header = TRUE)

  # Build metadata with GSM as rownames
  metadata <- data.frame(
    row.names = sample_meta$GSM,
    condition = factor(sample_meta$condition)
  )

  # Set healthy (or other reference) as reference level
  metadata$condition <- relevel(metadata$condition, ref = reference_label)

  # Check GSM IDs match count matrix columns
  stopifnot(all(rownames(metadata) %in% colnames(counts)))

  # Subset and reorder counts columns to match metadata
  counts <- counts[, rownames(metadata)]

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~ condition
  )

  # Filter low-count genes
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  # Run DESeq2
  dds <- DESeq(dds)

  # Get results: treatment vs reference
  res <- results(dds, contrast = c("condition", treatment_label, reference_label), alpha = 0.05)

  # Convert to dataframe
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)
  res_df <- res_df[, c("GeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

  # Merge gene symbol annotation
  annot_sub <- annot[, c("GeneID", "Symbol")]
  res_df$GeneID <- as.integer(res_df$GeneID)
  annot_sub$GeneID <- as.integer(annot_sub$GeneID)
  res_df <- merge(res_df, annot_sub, by = "GeneID", all.x = TRUE)

  # Reorder columns
  res_df <- res_df[, c("GeneID", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

  # Save
  write.csv(res_df, output_path, row.names = FALSE)
  message("Done: ", output_path)

  return(res_df)
}

# ---------------------------------------------------------------
# Gene annotation table (shared across all datasets)
# GRCh38.p13 NCBI annotation file (bundled with the GEO series'
# raw count matrices) — maps Entrez GeneID -> gene Symbol.
# Expects columns: GeneID, Symbol
# ---------------------------------------------------------------
annot <- read.csv("data/Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t")

# ---------------------------------------------------------------
# Dataset configuration — split by disease group so condition
# labels stay uniform within each group (treatment/ref should
# match the exact wording used in each dataset's condition column).
# Fill in metadata file paths after confirming filenames.
# ---------------------------------------------------------------

STROKE_TREATMENT <- "ischemic stroke"
STROKE_REFERENCE <- "healthy"

stroke_datasets <- list(
  list(id = "GSE122709",
       counts = "data/GSE122709_all.counts",
       meta = "data/GSE122709_sample_metadata.csv"),

  list(id = "GSE140275",
       counts = "data/GSE140275_raw_counts_GRCh38.p13_NCBI.tsv",
       meta = "data/GSE140275_sample_metadata.csv"),

  list(id = "GSE56267",
       counts = "data/GSE56267_raw_counts_GRCh38.p13_NCBI.tsv",
       meta = "data/GSE56267_sample_metadata.csv"),

  list(id = "GSE178764",
       counts = "data/GSE178764_raw_counts_GRCh38.p13_NCBI.tsv",
       meta = "data/GSE178764_sample_metadata.csv")
)

HF_TREATMENT <- "heart failure"   # TODO: confirm exact wording used in metadata
HF_REFERENCE <- "healthy"         # TODO: confirm exact wording used in metadata

heart_failure_datasets <- list(
  list(id = "GSE116250",
       counts = "data/GSE116250_raw_counts_GRCh38.p13_NCBI.tsv",
       meta = "data/GSE116250_sample_metadata.csv"),

  list(id = "GSE120852",
       counts = "data/GSE120852_raw_counts_GRCh38.p13_NCBI.tsv",
       meta = "data/GSE120852_sample_metadata.csv"),

  list(id = "GSE135055",
       counts = "data/GSE135055_raw_counts_GRCh38.p13_NCBI.tsv",
       meta = "data/GSE135055_sample_metadata.csv")
)

# ---------------------------------------------------------------
# Run DESeq2 — stroke group (uniform treatment/reference labels)
# ---------------------------------------------------------------
dir.create("02_deseq2/results", recursive = TRUE, showWarnings = FALSE)

stroke_results <- list()
for (d in stroke_datasets) {
  stroke_results[[d$id]] <- run_deseq2(
    counts_path      = d$counts,
    metadata_path    = d$meta,
    treatment_label  = STROKE_TREATMENT,
    reference_label  = STROKE_REFERENCE,
    annot            = annot,
    output_path      = paste0("02_deseq2/results/", d$id, "_deseq2.csv")
  )
}

# ---------------------------------------------------------------
# Run DESeq2 — heart failure group (uniform treatment/reference labels)
# ---------------------------------------------------------------
hf_results <- list()
for (d in heart_failure_datasets) {
  hf_results[[d$id]] <- run_deseq2(
    counts_path      = d$counts,
    metadata_path    = d$meta,
    treatment_label  = HF_TREATMENT,
    reference_label  = HF_REFERENCE,
    annot            = annot,
    output_path      = paste0("02_deseq2/results/", d$id, "_deseq2.csv")
  )
}

all_results <- c(stroke_results, hf_results)
