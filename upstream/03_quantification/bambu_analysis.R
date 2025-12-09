#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)
library(bambu)
library(Rsamtools)

# Define and parse command-line arguments
option_list <- list(
  make_option("--bam_files", type = "character", help = "Comma-separated list of BAM file paths"),
  make_option("--output_dir", type = "character", help = "Directory to save results"),
  make_option("--gtf", type = "character", help = "GTF annotation file"),
  make_option("--fasta", type = "character", help = "Genome FASTA file"),
  make_option("--threads", type = "integer", default = 4, help = "Number of threads [default: 4]")
)

args <- parse_args(OptionParser(option_list = option_list))

# Prepare annotation file to avoid re-processing
prepared_gtf <- prepareAnnotations(args$gtf)

# Split BAM file paths by comma and create BamFileList
bam_files <- unlist(strsplit(args$bam_files, ","))
stopifnot(all(file.exists(bam_files)))  # Ensure all BAM files exist

bamFileList <- Rsamtools::BamFileList(bam_files)

cat("Processing BAM files:\n", paste(bam_files, collapse = "\n"), "\n")

# Run bambu with error handling
tryCatch({
  se.multiSample <- bambu(
    reads = bamFileList,
    annotations = prepared_gtf,
    genome = args$fasta
  )

  # Save results
  save.dir <- args$output_dir
  writeBambuOutput(se.multiSample, path = save.dir, prefix = "Quartet_")

}, error = function(e) {
  cat("Error during bambu execution:\n", e$message, "\n")
  quit(status = 1)
})

cat("Bambu analysis completed successfully.\n")
