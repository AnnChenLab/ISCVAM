# Complete setup script for multiome pipeline
# This script installs R packages AND sets up MACS3 environment
# Run this script from the main project directory

cat("=== MULTIOME PIPELINE SETUP SCRIPT ===\n")
cat("Setting up R packages and MACS3 environment...\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Step 1: Install R packages
cat("STEP 1: Installing R packages...\n")
cat("=====================================\n")

# Install CRAN packages
cran_packages <- c(
  "dplyr",
  "future",
  "rhdf5"
)

cat("Installing CRAN packages...\n")
install.packages(cran_packages, dependencies = TRUE)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "Seurat",
  "Signac", 
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "SingleR"
)

cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_packages, dependencies = TRUE, update = FALSE)

# Check R package installation
cat("Checking R package installation...\n")
all_packages <- c(cran_packages, bioc_packages)
missing_packages <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) == 0) {
  cat("✅ All R packages installed successfully!\n\n")
} else {
  cat("❌ The following R packages failed to install:\n")
  cat(paste(missing_packages, collapse = ", "), "\n\n")
}

# Step 2: Setup MACS3 environment
cat("STEP 2: Setting up MACS3 environment...\n")
cat("=======================================\n")

# Check if conda is available
conda_available <- system("which conda", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0

if (!conda_available) {
  cat("❌ Conda not found. Please install Miniconda or Anaconda first.\n")
  cat("Download from: https://docs.conda.io/en/latest/miniconda.html\n")
} else {
  cat("✅ Conda found. Setting up MACS3 environment...\n")
  
  # Create conda environment and install MACS3
  cat("Creating conda environment 'macs3' with Python 3.9...\n")
  system("conda create -n macs3 python=3.9 -y", ignore.stdout = FALSE)
  
  cat("Installing MACS3 in the macs3 environment...\n")
  system("conda run -n macs3 pip install MACS3", ignore.stdout = FALSE)
  
  # Test MACS3 installation
  cat("Testing MACS3 installation...\n")
  macs3_test <- system("conda run -n macs3 macs3 --version", ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  if (macs3_test == 0) {
    cat("✅ MACS3 installed successfully!\n")
    
    # Get MACS3 path
    macs3_path <- system("conda run -n macs3 which macs3", intern = TRUE)
    cat("MACS3 path:", macs3_path, "\n")
    
    # Update the multiome_functions_pipeline.R file
    cat("Updating multiome_functions_pipeline.R with MACS3 path...\n")
    
    # Read the current file
    multiome_file <- "code_organized/core/multiome_functions_pipeline.R"
    if (file.exists(multiome_file)) {
      content <- readLines(multiome_file)
      
      # Find and replace the call.peaks_from_signac function
      start_idx <- grep("call.peaks_from_signac.*function", content)
      end_idx <- grep("^}", content[start_idx:length(content)])[1] + start_idx - 1
      
      if (length(start_idx) > 0 && length(end_idx) > 0) {
        # Create new function with correct MACS3 path
        new_function <- c(
          "call.peaks_from_signac<- function(seurat.raw, data.path){",
          "  DefaultAssay(seurat.raw) <- \"ATAC\"",
          paste0("  seurat.raw <- CallPeaks(seurat.raw, "),
          "                          outdir = data.path,",
          paste0("                          macs2.path = \"", macs3_path, "\")"),
          "  return(seurat.raw)",
          "}"
        )
        
        # Replace the function in the file
        new_content <- c(content[1:(start_idx-1)], new_function, content[(end_idx+1):length(content)])
        writeLines(new_content, multiome_file)
        cat("✅ Updated multiome_functions_pipeline.R with MACS3 path\n")
      }
    }
    
  } else {
    cat("❌ MACS3 installation failed\n")
  }
}

# Step 3: Final verification
cat("\nSTEP 3: Final verification...\n")
cat("=============================\n")

# Test R package loading
cat("Testing R package loading...\n")
test_packages <- c("Seurat", "Signac")
for (pkg in test_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat("✅", pkg, "loaded successfully\n")
  }, error = function(e) {
    cat("❌ Error loading", pkg, ":", e$message, "\n")
  })
}

# Summary
cat("\n=== SETUP SUMMARY ===\n")
if (length(missing_packages) == 0) {
  cat("✅ R packages: All installed\n")
} else {
  cat("❌ R packages: Some failed\n")
}

if (conda_available && macs3_test == 0) {
  cat("✅ MACS3 environment: Ready\n")
} else {
  cat("❌ MACS3 environment: Not ready\n")
}

cat("\n🎉 Setup complete! You can now run the multiome pipeline.\n")
cat("1. Configure dataset: source('code_organized/configure_dataset.R')\n")
cat("2. Run pipeline:      source('code_organized/core/run_pipeline_multisamples.R')\n")