#' Pseudocode for aligning raw reads to a specified GTF file
#' Uses Rsubread

#' Loads Rsubread
library(Rsubread)

#' Performs alignment
align(index="my_index",readfile1="sample1_R1.fastq",readfile2="sample1_R2.fastq",type="dna",
      output_file="sample1.bam", nthreads=24)
