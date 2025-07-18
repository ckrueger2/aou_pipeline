#!/usr/bin/R

#using biomart to find genomic coordinates of genes from spredixcan output
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

#data table and dplyr for data and table manipulation
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")

library(biomaRt)
library(data.table)
library(dplyr)
library(qqman)
library(argparse)

#set up argparse
parser <- ArgumentParser()
parser$add_argument("--phecode", help="all of us phenotype ID")
parser$add_argument("--pop", help="all of us population ID")
parser$add_argument("--ref", help="predixcan tissue")

args <- parser$parse_args()

#get the bucket name
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')

#pull file from bucket
name_of_file_in_bucket <- paste0(args$pop, "_predixcan_output_", args$phecode, "_", args$ref, ".csv")
read_in_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_file_in_bucket, " .")

#copy file from the bucket to the current workspace
system(read_in_command, intern=TRUE)
df <- fread(name_of_file_in_bucket, sep = ",", header=TRUE)

#accessing the Ensembl biomart database for 'genes', specifically the human genes version 113 (can be changed or left out)
biomart_access <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#removes the decimal after the gene, helps with merging and searching in general 
#this doesnt rename the gene, instead creates a new column so we dont lose data
df$gene_id <- sub("\\..*", "", df$gene)

#head to see the data and make sure it was read correctly 
cat("Preview of S-PrediXcan output \n")
head(df, 10)

#use biomart to pull chromosomal location information 
#we need the gene_id to find it, the chromosomal name (chromosome), the start position and end position (we will only use start) and what the gene is called in the database
gene_coords <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = df$gene_id,
  mart = biomart_access
)

cat("Preview of Biomart Query:\n")
head(gene_coords)

# attributes are the columns you want to retrieve/add from the database to your data
# filters are for ensuring that youre using the same data base as was used for the previous tools (spredixcan)
# values are what the database will use to search for your data, in our case the "cleaned up" gene_ids 
# mart is just the accessing that we created in the previous cell

#merge our TWAS table and biomart results 
merged_df <- left_join(df, gene_coords, by=c("gene_id" ="ensembl_gene_id"))
head(merged_df)

#merging the original "df", with new "gene_coords", gene_id in df and ensembl_gene_id in gene_coords into a new dataframe merged_df

#make P column is numeric and handle NA p-values and 0s
merged_df$P <- as.numeric(merged_df$pvalue)
zeros <- which(merged_df$P == 0)
#filter only for required fields
merged_df <- merged_df %>% filter(!is.na(chromosome_name) & !is.na(start_position) & !is.na(P))
if(length(zeros) > 0) {
  merged_df$P[merged_df$P == 0] <- 1e-300
}

#make CHR rows numeric and remove rows with NA
merged_df$CHR <- as.numeric(merged_df$chromosome_name)
merged_df <- merged_df[!is.na(merged_df$CHR), ]

#double check data
cat("Preview of BioMart merged table \n")
head(merged_df, 10)

#finding sample size to calculate threshold
sample_size <- nrow(df)

#calculate the new bonferroni and threshold based on sample size 
bonferroni_threshold <- 0.05 / sample_size
new_suggestive_threshold <- -log10(bonferroni_threshold)
cat("Bonferroni corrected P-value: ", bonferroni_threshold , "\n")

#manhattan plot title name
title <- paste0(args$pop, " ", args$phecode, " TWAS Manhattan Plot")

#name of saved file
destination_filename <- paste0(args$pop, "_", args$phecode,"_twas_manhattan_", args$ref, ".png")

#use qqman to plot the chromosome, location, snp, and pvalue into manhattan plot
png(filename = destination_filename, width = 1200, height = 800, res = 100)
manhattan(merged_df, 
          chr = "chromosome_name", 
          bp = "start_position", 
          snp = "gene_name", 
          p = "pvalue",
          main = title,
          col = c("lightblue3", "lightblue4"),
          suggestiveline = new_suggestive_threshold)
dev.off()

cat("TWAS Manhattan plot complete\n")