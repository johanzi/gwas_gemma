GWAS pipeline with GEMMA
=====

Simplified pipeline that requires from the user a VCF file with the samples of interest and a phenotype.

# Softwares needed
* Python 2.7
* R (>3.3.0) with the libraries indicated in the scripts
* Unix/Mac OS X environment for manipulation in bash
* p-link (v1.07)
* vcftools (v0.1.14)
* bcftools (v1.2)
* gemma (v0.94)

*NB: The versions indicated were the ones used and tested. One can first try already installed versions before installing the recommended version.*

# Files needed
* VCF file with the accessions/samples
* text file containing the order of the accessions in the VCF file (order_accession.txt)
* file containing the phenotypes of interest with a column containing the name of the accession
as described in the vcf file. Alternatively, one can process directly a tsv file containing
the phenotype for interest with one value per row, with the same order than for the VCF file (no header)



## VCF file preprocessing
Consider a VCF file containing 100 *Arabidopsis thaliana*, but only 80 accessions should be used in the GWAS. The VCF file must be first subsetted to these 80 accessions before being used in gemma. To do this, [vcftools](https://vcftools.github.io/man_latest.html) can be used.

```

# Subset the vcf file (list file contains the ID of each accession on different rows)
vcftools --keep list_accessions_to_keep.txt --gzvcf file.vcf.gz --recode --out subset_80

# The output file will be
subset_80.recode.vcf

# Compress and tabix the file
bgzip subset_80.recode.vcf && tabix subset_80.recode.vcf.gz

```

Get list of accessions in vcf file:

```
$ bcftool query -l  subset_80.recode.vcf.gz > order_accession.txt
$ cat order_accession.txt
1001
1002
1003

```

Note that the vcf file should also be filtered to remove singletons and low quality calls.

```
# Get singleton positions
vcftools --singletons --gvcf subset_80.recode.vcf.gz

# The file out.singletons should be generated. Then exclude these positions from the vcf file
vcftools --gvcf subset_80.recode.vcf --exclude-positions out.singletons \
 --recode --recode-INFO-all --out subset_80_without_singletons

# File generated subset_80_without_singletons.recode.vcf
# Compress and tabix file
bgzip subset_80_without_singletons.recode.vcf && tabix subset_80_without_singletons.recode.vcf.gz 

```

## Phenotype file

Example of content for phenotype file:

```
$ cat phenotype.tsv 
12.3
13.4
15.3
```
The value 12.3, 13.4, 15.3, ... being the height of the accessions 1001, 1002, and 1003, ..., respectively.

# General pipeline
1. Generate a dataframe with all the variables and order the accessions so that they match VCF sample order
2. Generate a file for one phenotype from the previously generated dataframe (test_export_df.txt)
3. Process the phenotype file with plink and gemma (run_gwas_gemma.sh). This step should be run directly into the directory containing the vcf file and the phenotype file.
4. Import in R the output phenotype.assoc.clean.txt file to visualize GWAS results (script gemma_analysis.Rmd)

* The part 1, 2, and 3 are done interactively in R (need to be adjusted according to the dataframe used)
* The part 3 is done in bash through the run_gwas_gemma.sh script. The only variables being the input vcf file used (change path in the script file) and the phenotype file given as first argument in command  line:
```
bash run_gwas_gemma.sh phenotype.tsv vcf_file.vcf
```

* The part 4 is done interactively in R

## Use a covariate
A strong peak can hide other peaks. In this case, the potential causative SNP in the peak can be used as covariate and the GWAS can be run again to assess what is the weight of the other SNPs when the covariate SNP weight is removed from the analysis.
To do so one needs to generate a file with 2 columns, the first containings 1s and the second a code for the SNP to use as covariate. For example, if the SNP can be coded as 1 and the reference allele as 0, therefore, a set of 4 accessions were only the 2 first accessions have the SNP would yield a covariate file looking like this:

```
$ cat covariate_file.txt
1	1
1	1
1	0
1	0
```

This file can then be used in gemma such as:

```
gemma -bfile ${prefix_vcf} -k ${current_path}/output/${prefix_vcf}.cXX.txt -lmm 2 -o ${prefix_gwas} -c covariate_cmt2.txt
```

TODO: implement covariate analysis in the script. Currently, one needs to proceed manually from the script and do the polishing step.





# Pipeline with ChromatinJ
If GWAS is to be performed in ChromatinJ output:
1. Generate a dataframe with all the variables and order the accessions so that they match VCF sample order with the R script process_chromatinJ_output.R (see example output file test_export_df.txt in directory 'examples')
2. Generate a file for one phenotype from the previously generated dataframe (test_export_df.txt) with the R script generate_phenotype_file.Rmd. The output is a file with a single column (no header) containing the phenotype of interest, scaled or not. See example files unscaled "Area_nucleus.tsv"
3. Process the phenotype file with plink and gemma (run_gwas_gemma.sh). This step should be run directly into the directory containing the vcf file and the phenotype file.
4. Import in R the output phenotype.assoc.clean.txt file to visualize GWAS results (script gemma_analysis.Rmd)

## User specific modifications

### process_chromatinJ_output.Rmd
-Indicate path of the ChromatinJ output results files
-Indicate path to file containing the order of the accessions from the VCF file

### generate_phenotype_file.Rmd (for the output of chromatinJ pipeline)
-Change the path of the file to import and indicate which phenotype is wanted and which name for this phenotype should be given

### run_gwas_gemma.sh
-Provide the two arguments (phenotype file and VCF file). Note that the two input files should be located in the same directory.

### gemma_analysis.Rmd
-Change dir_file and file.name variables

# Authors
* **Johan Zicola** - [johanzi](https://github.com/johanzi)

# License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details



