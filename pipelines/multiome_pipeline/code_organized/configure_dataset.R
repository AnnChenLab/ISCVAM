# ============================================================================
# DATASET CONFIGURATION SCRIPT
# ============================================================================
# Use this script to easily run the pipeline with different datasets
# Simply source this file before running the main pipeline

# === Brain 3k Dataset ===
# dataset_name <- "multiome_brain_3k"
# species <- "Human"

dataset_name <- "multiome_human_kidney"
species <- "Human"
# ============================================================================
# AUTO-GENERATED PATHS (DO NOT MODIFY)
# ============================================================================

# For single sample: data path to raw_data folder
data.path <- paste0("example_data/", dataset_name, "/raw_data")

# For multi-sample: root directory containing sample subfolders (GSM...)
# Used by run_pipeline_multisamples.R
root_dir <- paste0("example_data/", dataset_name, "/")

results_dir <- paste0("results/results_", dataset_name)
filename <- paste0(dataset_name, "_analysis.h5")

# Load SingleR references
if (file.exists("code_organized/ref_data/singler.refs.rdata")) {
  load("code_organized/ref_data/singler.refs.rdata")
  cat("✅ SingleR references loaded\n")
} else {
  cat("⚠️  WARNING: SingleR references not found at code_organized/ref_data/singler.refs.rdata\n")
}

cat("=== DATASET CONFIGURATION LOADED ===\n")
cat("Dataset:", dataset_name, "\n")
cat("Species:", species, "\n") 
cat("Data path:", data.path, "\n")
cat("Results directory:", results_dir, "\n")
cat("Output filename:", filename, "\n")
cat("=====================================\n\n")

# Check if dataset directory exists
if (!dir.exists(data.path)) {
  cat("⚠️  WARNING: Dataset directory not found:", data.path, "\n")
  cat("   Please make sure you have the dataset in the correct location.\n\n")
} else {
  cat("✅ Dataset directory found:", data.path, "\n\n")
}

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("✅ Created results directory:", results_dir, "\n\n")
} else {
  cat("✅ Results directory exists:", results_dir, "\n\n")
}

cat("Ready to run pipeline!\n")
cat("source('code_organized/core/run_pipeline.R')\n")