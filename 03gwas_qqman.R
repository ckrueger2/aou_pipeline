#!/usr/bin/R

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")

library(data.table)
library(qqman)
library(argparse)

#set up argparse
parser <- ArgumentParser()
parser$add_argument("--phecode", help="all of us phenotype ID")
parser$add_argument("--pop", help="all of us population ID")
parser$add_argument("--snp_count", type = "integer", help="number of SNPs in summary statistics table")

args <- parser$parse_args()

#find bucket
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')

#build file name
name_of_file_in_bucket <- paste0(args$pop, "_formatted_filtered_", args$phecode,".tsv")

#copy file from the bucket to the current workspace
system(paste0("gsutil cp ", my_bucket, "/data/", name_of_file_in_bucket, " ."), intern=T)

#load the file into a dataframe
hail_table  <- fread(name_of_file_in_bucket, header=TRUE)

#what to save file as
destination_filename <- paste0(args$pop, "_", args$phecode,"_gwas_manhattan.png")

#manhattan plot title name
title <- paste0(args$pop, " ", args$phecode, " GWAS Manhattan Plot")

#png destination
png(filename = destination_filename, width = 1200, height = 800, res = 100)

#create the Manhattan plot; only use Bonferroni line if SNP count is provided
if (!is.null(args$snp_count)) {
  #calc Bonferroni correction value
  b_value <- (0.05/args$snp_count)
  log_b_value <- -log10(b_value)
  cat("Bonferroni corrected P-value: ", b_value, "\n")
  
  #manhattan plot with suggestive line
  manhattan(hail_table,
            main = title,
            chr = "CHR", 
            bp = "POS", 
            p = "Pvalue", 
            snp = "SNP", 
            col = c("lightblue3", "lightblue4"),
            suggestiveline = log_b_value)
} else {
  #manhattan plot without suggestive line
  cat("SNP count not provided, plotting without Bonferroni correction line\n")
  manhattan(hail_table,
            main = title,
            chr = "CHR", 
            bp = "POS", 
            p = "Pvalue", 
            snp = "SNP", 
            col = c("lightblue3", "lightblue4"),
            suggestiveline = FALSE)
}
dev.off()

cat("GWAS Manhattan plot complete\n")