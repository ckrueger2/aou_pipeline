#!/usr/bin/R

cat("Beginning table formatting...")

#if needed, install packages
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")

#load packages
library(data.table)
library(dplyr)
library(tidyr)
library(argparse)

#set up argparse
parser <- ArgumentParser()
parser$add_argument("--phecode", help="all of us phenotype ID")
parser$add_argument("--pop", help="all of us population ID")

args <- parser$parse_args()

#find bucket
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')

#PERFORM COMMAND LINE FORMATTING FOR rsID VCF FILE
#download reference genome
file_name <- "All_20180418.vcf.gz"
#check if file exists in the home directory
if (!file.exists(file_name)) {
  #if file doesn't exist, download it
  command <- paste0("wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/", file_name)
  system(command)
} else {
  #if file exists, print a message
  cat("File", file_name, "already exists. Skipping download.\n")
}

cat("Formatting reference files...")

#create file of chr and pos columns only to use for filtering
command2 <- paste0("gsutil cat ", my_bucket, "/data/", args$pop, "_filtered_", args$phecode, ".tsv | awk 'NR > 1 {print $8, $9}' > /tmp/subset_", args$phecode, ".tsv")
system(command2)

#remove chr prefix
command3 <- paste0("sed -e 's/chr//' -e 's/^X /23 /' /tmp/subset_", args$phecode, ".tsv > /tmp/nochr", args$phecode, ".tsv")
system(command3)

#filter large file, eliminating SNPs not present in sumstats file
command4 <- paste0("zcat All_20180418.vcf.gz | awk 'NR==FNR {a[$1\" \"$2]=1; next} !/^#/ && ($1\" \"$2) in a' /tmp/nochr", args$phecode, ".tsv - > /tmp/filtered_20180418.vcf")
system(command4)

#remove metadata rows
command5 <- paste0("awk '!/^##/' /tmp/filtered_20180418.vcf > /tmp/", args$phecode, "ref.vcf")
system(command5)

#copy to bucket
command6 <- paste0("gsutil cp /tmp/", args$phecode, "ref.vcf ", my_bucket, "/data/")
system(command6)

#check bucket for vcf file
check_result <- system(paste0("gsutil ls ", my_bucket, "/data/ | grep ", args$phecode, "ref.vcf"), ignore.stderr = TRUE)

if (check_result != 0) {
  stop(paste0("ERROR: File '", args$phecode, "ref.vcf' was not found in bucket ", my_bucket, "/data/"))
} else {
  cat("Reference VCF file successfully transferred to bucket.\n")
}

#PERFORM COMMAND LINE FORMATTING FOR S-PREDIXCAN FILE
#upload GTEx SNP file to workspace bucket
command7 <- paste0("gsutil -m cp -v ~/aou_predixcan/predixcan_models_varids-effallele_mesa.txt.gz ", my_bucket, "/data/")
system(command7, intern=TRUE)
#command7.5 <- paste0("gsutil -m cp -v ~/aou_predixcan/", args$pop, "_mesa_snps.txt.gz ", my_bucket, "/data/")
#system(command7.5, intern=TRUE)

#unzip files
command8 <- paste0("gsutil cat ", my_bucket, "/data/predixcan_models_varids-effallele_mesa.txt.gz | gunzip > /tmp/predixcan_models_varids-effallele_mesa.txt")
system(command8)
#command8.5 <- paste0("gsutil cat ", my_bucket, "/data/", args$pop, "_mesa_snps.txt.gz | gunzip > /tmp/", args$pop, "_mesa_snps.txt")
#system(command8.5)

#format reference file
system("awk -F'[,:]' 'NR>1 {print $1\":\"$2}' /tmp/predixcan_models_varids-effallele_mesa.txt > /tmp/chrpos_allele_table.tsv", intern=TRUE)

#make temp files
command9 <- paste0("gsutil cp ", my_bucket, "/data/", args$pop, "_full_", args$phecode,".tsv /tmp/")
system(command9)

#filter SNPs
command10 <- paste0("awk 'NR==FNR{a[$1];next} $1 in a' /tmp/chrpos_allele_table.tsv /tmp/", args$pop, "_full_", args$phecode, ".tsv > /tmp/gtex_", args$phecode, ".tsv")
system(command10)
#command10.5 <- paste0("awk 'NR==FNR{if(NR>1){split($0,arr,\",\"); a[arr[2]\":\"arr[3]]}; next} FNR>1 && $1 in a' /tmp/", args$pop, "_mesa_snps.txt /tmp/", args$pop, "_full_", args$phecode, ".tsv > /tmp/mesa_", args$phecode, ".tsv")
#system(command10.5)

#save to bucket
command11 <- paste0("gsutil cp /tmp/gtex_", args$phecode, ".tsv ", my_bucket, "/data/", args$pop, "_gtex_", args$phecode,".tsv")
system(command11)
#command11.5 <- paste0("gsutil cp /tmp/mesa_", args$phecode, ".tsv ", my_bucket, "/data/", args$pop, "_mesa_", args$phecode,".tsv")
#system(command11.5)

#check bucket
check_result2 <- system(paste0("gsutil ls ", my_bucket, "/data/ | grep ", args$pop, "_gtex_", args$phecode, ".tsv"), ignore.stderr = TRUE)
#check_result3 <- system(paste0("gsutil ls ", my_bucket, "/data/ | grep ", args$pop, "_mesa_", args$phecode, ".tsv"), ignore.stderr = TRUE)

if (check_result2 != 0) {
  stop(paste0("ERROR: File '", args$pop, "_gtex_", args$phecode, ".tsv' was not found in bucket ", my_bucket, "/data/"))
} else {
  cat("Reference GTEx file successfully transferred to bucket.\n")
}

#if (check_result3 != 0) {
#  stop(paste0("ERROR: File '", args$pop, "_mesa_", args$phecode, ".tsv' was not found in bucket ", my_bucket, "/data/"))
#} else {
#  cat("Reference mesa file successfully transferred to bucket.\n")
#}

#FORMAT TABLES
#read in gtex filtered table
name_of_gtex_file <- paste0(args$pop, "_gtex_", args$phecode, ".tsv")
gtex_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_gtex_file, " .")

system(gtex_command, intern=TRUE)

gtex_table <- fread(name_of_gtex_file, header=FALSE, sep="\t")
colnames(gtex_table) <- c("locus","alleles","BETA","SE","Het_Q","Pvalue","Pvalue_log10","CHR","POS","rank","Pvalue_expected","Pvalue_expected_log10")

#check table
cat("GTEx filtered table preview:\n")
head(gtex_table)

#read in mesa filtered table
#name_of_mesa_file <- paste0(args$pop, "_mesa_", args$phecode, ".tsv")
#mesa_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_mesa_file, " .")

#system(mesa_command, intern=TRUE)

#mesa_table <- fread(name_of_mesa_file, header=FALSE, sep="\t")
#colnames(mesa_table) <- c("locus","alleles","BETA","SE","Het_Q","Pvalue","Pvalue_log10","CHR","POS","rank","Pvalue_expected","Pvalue_expected_log10")

#check table
#cat("MESA filtered table preview:\n")
#head(mesa_table)

#read in pvalue filtered table
name_of_filtered_file <- paste0(args$pop, "_filtered_", args$phecode, ".tsv")
filtered_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_filtered_file, " .")

system(filtered_command, intern=TRUE)

filtered_table <- fread(name_of_filtered_file, header=TRUE)

#check table
cat("pvalue filtered table preview:\n")
head(filtered_table)

#GTEX TABLE
#reformat locus column to chr_pos_ref_alt_b38
gtex_table$locus_formatted <- gsub(":", "_", gtex_table$locus) #colon to underscore
gtex_table$alleles_formatted <- gsub('\\["', "", gtex_table$alleles)  #remove opening [
gtex_table$alleles_formatted <- gsub('"\\]', "", gtex_table$alleles_formatted)  #remove closing ]
gtex_table$alleles_formatted <- gsub('","', "_", gtex_table$alleles_formatted)  #comma to underscore

#split allele column
gtex_table <- gtex_table %>%
  separate(alleles_formatted, into = c("REF", "ALT"), sep = "_", remove=F)

#combine strings
gtex_table$SNP <- paste0(gtex_table$locus_formatted, "_", gtex_table$alleles_formatted, "_b38")
gtex_table$ID <- paste0(gtex_table$locus, ":", gtex_table$REF, ":", gtex_table$ALT)

#remove intermediate columns
gtex_table$locus_formatted <- NULL
gtex_table$alleles_formatted <- NULL

#edit sex chromosomes
gtex_table$CHR <- gsub("X", "23", gtex_table$CHR)
gtex_table$CHR <- gsub("Y", "24", gtex_table$CHR)

#MESA TABLE
#reformat locus column to chr_pos_ref_alt_b38
#mesa_table$locus_formatted <- gsub(":", "_", mesa_table$locus) #colon to underscore
#mesa_table$alleles_formatted <- gsub('\\["', "", mesa_table$alleles)  #remove opening [
#mesa_table$alleles_formatted <- gsub('"\\]', "", mesa_table$alleles_formatted)  #remove closing ]
#mesa_table$alleles_formatted <- gsub('","', "_", mesa_table$alleles_formatted)  #comma to underscore

#split allele column
#mesa_table <- mesa_table %>%
#  separate(alleles_formatted, into = c("REF", "ALT"), sep = "_", remove=F)

#combine strings
#mesa_table$SNP <- paste0(mesa_table$locus_formatted, "_", mesa_table$alleles_formatted, "_b38")
#mesa_table$ID <- paste0(mesa_table$locus, ":", mesa_table$REF, ":", mesa_table$ALT)

#remove intermediate columns
#mesa_table$locus_formatted <- NULL
#mesa_table$alleles_formatted <- NULL

#edit sex chromosomes
#mesa_table$CHR <- gsub("X", "23", mesa_table$CHR)
#mesa_table$CHR <- gsub("Y", "24", mesa_table$CHR)

#FILTERED TABLE
filtered_table$locus_formatted <- gsub(":", "_", filtered_table$locus) #colon to underscore
filtered_table$alleles_formatted <- gsub('\\["', "", filtered_table$alleles)  #remove opening [
filtered_table$alleles_formatted <- gsub('"\\]', "", filtered_table$alleles_formatted)  #remove closing ]
filtered_table$alleles_formatted <- gsub('","', "_", filtered_table$alleles_formatted)  #comma to underscore

filtered_table <- filtered_table %>%
  separate(alleles_formatted, into = c("REF", "ALT"), sep = "_", remove=F)

#filtered_table$SNP <- paste0(filtered_table$locus_formatted, "_", filtered_table$alleles_formatted, "_b38")
#filtered_table$ID <- paste0(filtered_table$locus, ":", filtered_table$REF, ":", filtered_table$ALT)
filtered_table$SNP <- paste0(filtered_table$locus_formatted, "_b38")
filtered_table$ID <- paste0(filtered_table$locus)

filtered_table$locus_formatted <- NULL
filtered_table$alleles_formatted <- NULL

filtered_table$CHR <- gsub("X", "23", filtered_table$CHR)
filtered_table$CHR <- gsub("Y", "24", filtered_table$CHR)

#MERGE rsIDs TO S-PREDIXCAN TABLE
#read in rsID reference file
name_of_vcf <- paste0(args$phecode, "ref.vcf")
reference_command <- paste0("gsutil cp ", my_bucket, "/data/", name_of_vcf, " .")

system(reference_command, intern=T)

reference_data <- fread(name_of_vcf, header = FALSE, sep='\t')
reference_data <- reference_data[,1:3]
colnames(reference_data) <- c("CHR", "POS", "rsID")

#format data for matching
filtered_table$CHR <- as.character(filtered_table$CHR)
filtered_table$POS <- as.character(filtered_table$POS)

reference_data$CHR <- paste0("chr", reference_data$CHR)
reference_data$CHR <- as.character(reference_data$CHR)
reference_data$POS <- as.character(reference_data$POS)

#check tables
cat("rsID reference table preview:\n")
head(reference_data)

#merge files
merged_table <- merge(filtered_table, reference_data[, c("CHR", "POS", "rsID")], by = c("CHR", "POS"), all.x = TRUE)
head(merged_table)

#remove un-needed columns
filtered_merged_table <- merged_table[, c(1, 2, 13, 14, 15, 17, 5, 6, 8, 3)]

#check table
cat("rsID merged table preview:\n")
head(filtered_merged_table)

#FINAL FORMATTING
#format chromosomes
filtered_merged_table$CHR <- gsub("chr", "", filtered_merged_table$CHR)
gtex_table$CHR <- gsub("chr", "", gtex_table$CHR)
gtex_table$CHR <- gsub("X", "23", gtex_table$CHR)
gtex_table$CHR <- gsub("Y", "24", gtex_table$CHR)
#mesa_table$CHR <- gsub("chr", "", mesa_table$CHR)
#mesa_table$CHR <- gsub("X", "23", mesa_table$CHR)
#mesa_table$CHR <- gsub("Y", "24", mesa_table$CHR)

#make numeric
filtered_merged_table$CHR <- as.numeric(filtered_merged_table$CHR)
gtex_table$CHR <- as.numeric(gtex_table$CHR)
#mesa_table$CHR <- as.numeric(mesa_table$CHR)
filtered_merged_table$POS <- as.numeric(filtered_merged_table$POS)

#read in mesa rsIDs
#mesa_data <- fread("/tmp/", args$pop, "_mesa_snps.txt", header=TRUE, sep=",")

#this merge is finicky, re-format everything
#mesa_data$chr <- paste0("chr", gsub("^chr", "", mesa_data$chr))
#mesa_table$CHR <- paste0("chr", gsub("^chr", "", mesa_table$CHR))
#mesa_data$pos <- as.integer(as.character(mesa_data$pos))
#mesa_table$POS <- as.integer(as.character(mesa_table$POS))
#mesa_data$chr <- trimws(mesa_data$chr)
#mesa_data$rsid <- trimws(mesa_data$rsid)
#mesa_table$CHR <- trimws(mesa_table$CHR)

#create unique key for matching 
#mesa_data$merge_key <- paste(mesa_data$chr, mesa_data$pos, sep="_")
#mesa_table$merge_key <- paste(mesa_table$CHR, mesa_table$POS, sep="_")

#merge mesa table and data
#merged_mesa_table <- merge(mesa_table, mesa_data[, c("merge_key", "rsid")], by = "merge_key", all.x = TRUE)

#clean up the merge key column
#merged_mesa_table$merge_key <- NULL

#sort by chr, pos
filtered_merged_table <- filtered_merged_table %>%
  arrange(CHR, POS)
gtex_table <- gtex_table %>%
  arrange(CHR, POS)
#merged_mesa_table <- merged_mesa_table %>%
#  arrange(CHR, POS)

#rename header
gtex_table$"#CHROM" <- gtex_table$CHR
gtex_table$CHR <- NULL
#merged_mesa_table$"#CHROM" <- merged_mesa_table$CHR
#merged_mesa_table$CHR <- NULL

#select columns
gtex_table <- gtex_table %>%
  select(locus, alleles, ID, REF, ALT, "#CHROM", BETA, SE, Pvalue, SNP)
#merged_mesa_table <- merged_mesa_table %>%
#  select(locus, alleles, ID, REF, ALT, "#CHROM", BETA, SE, Pvalue, SNP, rsid, ID)

#edit X chromosome SNP file for gtex table
gtex_table$SNP <- gsub("^chrX_", "X_", gtex_table$SNP)

#edit locus column for pvalue filtered table
filtered_merged_table$locus <- gsub("^chrX:", "X:", filtered_merged_table$locus)

#check tables
cat("Final pvalue filtered table:\n")
head(filtered_merged_table)

cat("Final GTEx filtered table:\n")
head(gtex_table)

#cat("Final MESA filtered table:\n")
#head(merged_mesa_table)

#write gtex table
gtex_destination_filename <- paste0(args$pop, "_formatted_gtex_", args$phecode,".tsv")

#store the dataframe in current workspace
write.table(gtex_table, gtex_destination_filename, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#copy the file from current workspace to the bucket
system(paste0("gsutil cp ./", gtex_destination_filename, " ", my_bucket, "/data/"), intern=TRUE)

#write mesa table
#mesa_destination_filename <- paste0(args$pop, "_formatted_mesa_", args$phecode,".tsv")

#store the dataframe in current workspace
#write.table(merged_mesa_table, mesa_destination_filename, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#copy the file from current workspace to the bucket
#system(paste0("gsutil cp ./", mesa_destination_filename, " ", my_bucket, "/data/"), intern=TRUE)

#write p-filtered table
filtered_destination_filename <- paste0(args$pop, "_formatted_filtered_", args$phecode,".tsv")

#store the dataframe in current workspace
write.table(filtered_merged_table, filtered_destination_filename, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#copy the file from current workspace to the bucket
system(paste0("gsutil cp ./", filtered_destination_filename, " ", my_bucket, "/data/"), intern=TRUE)

#CHECK IF FILES ARE IN THE BUCKET
#GTEx file
check_gtex <- system(paste0("gsutil ls ", my_bucket, "/data/ | grep ", gtex_destination_filename), ignore.stderr = TRUE)

if (check_gtex != 0) {
  stop(paste0("ERROR: File '", gtex_destination_filename, "' was not found in bucket ", my_bucket, "/data/"))
} else {
  cat("GTEx formatted file successfully saved to bucket.\n")
}

#mesa file
#check_mesa <- system(paste0("gsutil ls ", my_bucket, "/data/ | grep ", mesa_destination_filename), ignore.stderr = TRUE)

#if (check_mesa != 0) {
#  stop(paste0("ERROR: File '", mesa_destination_filename, "' was not found in bucket ", my_bucket, "/data/"))
#} else {
#  cat("MESA formatted file successfully saved to bucket.\n")
#}

#filtered file
check_filtered <- system(paste0("gsutil ls ", my_bucket, "/data/ | grep ", filtered_destination_filename), ignore.stderr = TRUE)

if (check_filtered != 0) {
  stop(paste0("ERROR: File '", filtered_destination_filename, "' was not found in bucket ", my_bucket, "/data/"))
} else {
  cat("Filtered pvalue formatted file successfully saved to bucket.\n")
}

#clean up tmp files
system(paste0("rm -f /tmp/subset_", args$phecode, ".tsv /tmp/nochr", args$phecode, ".tsv /tmp/filtered_20180418.vcf /tmp/", args$phecode, "ref.vcf /tmp/predixcan_models_varids-effallele_mesa.txt /tmp/chrpos_allele_table.tsv /tmp/", args$pop, "_full_", args$phecode, ".tsv /tmp/gtex_", args$phecode, ".tsv"), intern=TRUE)
