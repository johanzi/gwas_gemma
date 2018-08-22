# Format df 

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("4 arguments need to be provided (accession order file, accessions details, methylation data, output file name).n", call.=FALSE)
} else if (length(args) == 4) {
  df_order <- args[1]
  df_accession <- args[2]
  df_methylation <- args[3]
  name_output <- args[4]
}
  
# Create the final file for GWAS, the samples need to be aligned as precised in the vcf file
df_order <- read.table(df_order, header=TRUE)

df_accession <- read.table(df_accession, header=TRUE)

df_methylation <- read.table(df_methylation, header=TRUE)

# First merge
df_order_accession <- merge(df_order, df_accession, by="seq_ID", sort=FALSE)

# Second merge
df_meth_accession <- merge(df_order_accession, df_methylation, by="sample", sort=FALSE)

# Export file for GWAS (no header, only phenotype)
write.table(df_meth_accession$percent_methylation, file=name_output, row.names=FALSE, col.names=FALSE, quote=FALSE)
