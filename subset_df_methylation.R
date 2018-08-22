# Load the methylation data in format rda (R object)

# Function can be called in Rscript with command line argument
# Rscript subset_df_methylation.R CpG df_mean_genes.rda CpG_filtered_test.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("3 arguments need to be provided (methylation context, name of rda object, name output text file).n", call.=FALSE)
} else if (length(args) == 3) {
  context <- args[1]
  rda_object <- args[2]
  name_output <- args[3]
}

df <- load(rda_object)
df <- get(df)

# Define context of interest
df_context <- df[(df$context == context),]


# Export the new df
write.table(df_context, file=name_output, row.names=FALSE, quote=FALSE, sep="\t")

