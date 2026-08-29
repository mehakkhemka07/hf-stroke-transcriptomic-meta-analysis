library(GEOquery)

# ---------------------------------------------------------------
# Reusable metadata downloader
# ---------------------------------------------------------------
download_geo_metadata <- function(geo_id, out_dir_base = "01_datasets/metadata", data_dir = "data") {

  out_dir <- file.path(out_dir_base, paste0(geo_id, "_metadata"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  message(">>> Fetching metadata for ", geo_id, " ...")
  gse <- getGEO(geo_id, GSEMatrix = FALSE, destdir = out_dir)

  # ---- Series-level info ----
  cat("\n========== SERIES INFO:", geo_id, "==========\n")
  cat("Title         :", Meta(gse)$title, "\n")
  cat("Summary       :", Meta(gse)$summary, "\n")
  cat("Overall design:", Meta(gse)$overall_design, "\n")
  cat("Organism      :", Meta(gse)$sample_organism, "\n")
  cat("Platform      :", Meta(gse)$platform_id, "\n")
  cat("Submission    :", Meta(gse)$submission_date, "\n")
  cat("Pubmed ID     :", Meta(gse)$pubmed_id, "\n")
  cat("n Samples     :", length(GSMList(gse)), "\n")

  # ---- Sample-level metadata ----
  gsm_list <- GSMList(gse)

  flatten <- function(x) paste(x, collapse = "; ")

  sample_meta <- lapply(names(gsm_list), function(gsm_name) {
    gsm  <- gsm_list[[gsm_name]]
    meta <- Meta(gsm)

    data.frame(
      GSM              = gsm_name,
      title            = flatten(meta$title),
      source_name      = flatten(meta$source_name_ch1),
      organism         = flatten(meta$organism_ch1),
      characteristics  = flatten(meta$characteristics_ch1),
      molecule         = flatten(meta$molecule_ch1),
      extract_protocol = flatten(meta$extract_protocol_ch1),
      data_processing  = flatten(meta$data_processing),
      platform         = flatten(meta$platform_id),
      stringsAsFactors = FALSE
    )
  })

  meta_df <- do.call(rbind, sample_meta)

  # ---- Print summary to console ----
  cat("\n========== UNIQUE CONDITIONS:", geo_id, "==========\n")
  cat("Source names:\n");       print(unique(meta_df$source_name))
  cat("\nCharacteristics:\n");  print(unique(meta_df$characteristics))
  cat("\nOrganisms:\n");        print(unique(meta_df$organism))
  cat("\nPlatforms:\n");        print(unique(meta_df$platform))
  cat("\nMolecule types:\n");   print(unique(meta_df$molecule))

  # ---- Save to CSV ----
  # Raw/verbose copy (kept alongside downloaded series info)
  out_file <- file.path(out_dir, "sample_metadata.csv")
  write.csv(meta_df, out_file, row.names = FALSE)

  # Flat copy used by run_deseq2.R, named "<GEO_ID>_sample_metadata.csv" in data/
  deseq_file <- file.path(data_dir, paste0(geo_id, "_sample_metadata.csv"))
  write.csv(meta_df, deseq_file, row.names = FALSE)

  message(">>> Metadata saved to: ", out_file, " and ", deseq_file, "\n")

  return(meta_df)
}

# ---------------------------------------------------------------
# Run for all 7 datasets used in this project
# ---------------------------------------------------------------
geo_ids <- c(
  "GSE122709",  # Stroke
  "GSE140275",  # Stroke
  "GSE56267",   # Stroke
  "GSE178764",  # Stroke
  "GSE116250",  # Heart Failure
  "GSE120852",  # Heart Failure
  "GSE135055"   # Heart Failure
)

all_metadata <- lapply(geo_ids, download_geo_metadata)
names(all_metadata) <- geo_ids
