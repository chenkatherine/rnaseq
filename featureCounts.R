#' Generates count matrices from .bam files

#' Load libraries
#' "Rsubread" for featureCounts
#' "openxlsx" to write results to disk
library(Rsubread)
library(openxlsx)

#' Sets input and output directories and annotation file
bam_directory <- "BAM_DIRECTORY"
annotation_file <- "protein_coding_exons.gtf"
output_directory <- "PROCODING_COUNT_MATRICES"

#' Creates output directory if it doesn't already exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

#' Gets all BAM files in the specified BAM directory
bam_files <- list.files(path = bam_directory, 
                        pattern = "\\.bam$", 
                        full.names = TRUE)

#' Stops if no BAM files were found in the specified directory
if (length(bam_files) == 0) {
  stop("No BAM files found in the specified directory")
}

#' Processes each BAM file
for (bam_file in bam_files) {
  
  #' Extracts basename
  base_name <- tools::file_path_sans_ext(basename(bam_file))
  
  #' Prints progress
  cat("Processing:", base_name, "\n")
  
  #' Runs featureCounts on current BAM file
  fc <- featureCounts(files = bam_file,
                      annot.ext = annotation_file,
                      isGTFAnnotationFile = TRUE,
                      isPairedEnd = TRUE)
  
  #' Creates count matrix output 
  output_file <- file.path(output_directory, paste0(base_name, ".csv"))
  
  #' Writes count matrix to csv
  write.table(x = data.frame(fc$annotation[, c("GeneID")],
                             fc$counts,
                             stringsAsFactors = FALSE),
              file = output_file,
              quote = FALSE,
              sep = ",",
              row.names = FALSE)
  
  #' Completion message
  cat("Completed:", base_name, "-> saved as", output_file, "\n")
}

cat("All files processed successfully!\n")