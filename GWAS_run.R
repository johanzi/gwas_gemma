GWAS_run <- function(output_gemma, threshold_pvalue="0", highlighted_SNP=""){
  
  # Highlighted_SNP allows to display in green the SNP of interested on the Manahattan plot
  # It can be 1 SNP (e.g. highlighted_SNP="4:10420088") or several SNPs, passed as a vector
  # (e.g. highlighted_SNP=c("4:10420088","5:112000"). No SNP highlighted by default
  
  # Import GEMMA output file
  gwas.results <- read.delim(path.file, sep="\t")
  
  # Plot QQ plot (need to precise the package as lattice has a similar function
  #qqman::qq(gwas.results$P, main=file.name)
  
  # One can select SNPs above the Bonferroni corrected p-value threshold
  # by using the argument "bonferroni"
  if(threshold_pvalue == "bonferroni"){
    # Calculate Bonferroni threshold with risk 5%
    ## Get total number of SNPs
    nb_snps <- dim(gwas.results)[[1]]
    
    ## Calculate Bonferroni corrected P-value threshold
    bonferroni_threshold <- 0.05/nb_snps
    
    threshold_pvalue <- bonferroni_threshold
  } else {
    # In case the variable was entered as string and is not "bonferroni"
    # convert to numeric. Set to 0 by default if user does not want any threshold
    threshold_pvalue <-  as.numeric(threshold_pvalue)
  }
  
  # Get positions of the chromosome with SNPs having a -log(P) > 5
  gwas_significant <- subset(gwas.results, P < threshold_pvalue)
  
  # Default p-value threshold line commonly used in GWAS -> -log10(5e-8) => red line. 
  # Set genomewideline to False has it makes little sense for Arabidopsis genome
  
  # suggestive line = Bonferroni corrected P-value threshold => blue line
  
  # Plot manhattan plot
  manhattan(gwas.results, highlight=highlighted_SNP, main=file.name, suggestiveline = -log10(threshold_pvalue), genomewideline = FALSE)
  
  #Check if dataframe is not empty (no SNPs above threshold value
  if(dim(gwas_significant)[[1]] != 0){ 
    return(gwas_significant)
  }
}


