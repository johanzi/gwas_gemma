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
Consider a VCF file containing 100 *Arabidopsis thaliana*, but the phenotype of only 80 accessions is available. The VCF file must then be first subset to these 80 accessions before being used in gemma. To do this, [vcftools](https://vcftools.github.io/man_latest.html) can be used. The VCF file input for GWAS should not contains indels, singletons, and keep only biallelic positions (non-alternative position can be removed to reduce file size). Also, the calls should be filtered by their quality for every individual, for instance a quality of minimum 25 (GQ>=25) and a coverage of minimum 3 reads (DP>=3) to keep the genotype for a certain SNP and individual (GT field in VCF).

```

# Subset the vcf file (list file contains the ID of each accession on different rows)
vcftools --keep list_accessions_to_keep.txt --gzvcf file.vcf.gz --recode recode-INFO-all --out subset_80

# The output file will be
subset_80.recode.vcf

# VCF file can be heavy and usually most of lines do not contain information (no ALT allele), remove then
# now to reduce size of the file and go faster with next step
# Keep only positions with an alternative allele (--min-ac) and only biallelic positions (--max-alleles) (1 REF + 1 ALT)
bcftools view --min-ac=1 --max-alleles 2  subset_80.recode.vcf > subset_80_biallelic_only_alt.recode.vcf



# Remove indels and hide GT of individuals with low quality call (DP>=3 and GQ>=25)
vcftools --vcf subset_80_biallelic_only_alt.recode.vcf --remove-indels --minDP 3 --minGQ 25 \
			 --recode --recode-INFO-ALL --out subset_80_biallelic_only_alt_no_indels

# Remove singletons
## Generates out.singletons (positions of all singletons)
vcftools --singletons --vcf subset_80_biallelic_only_alt_no_indels.recode.vcf

## Exclude singleton positions 
vcftools --vcf subset_80_biallelic_only_alt_no_indels.recode.vcf \
			--exclude-positions out.singletons --recode --recode-INFO-all \
			--out subset_80_biallelic_only_alt_no_indels_no_singletons 

# Compress and tabix the file
bgzip  subset_80_biallelic_only_alt_no_indels_no_singletons.recode.vcf && \
			tabix subset_80_biallelic_only_alt_no_indels_no_singletons.recode.vcf.gz 

# Note: If you have chromosomes or organelle genomes to be excluded (mitochondria, chloroplasts, ...), 
# you should remove them as they increase the number of SNPs to be tested and therefore 
# decrease the threshold of significance (if Bonferroni correction is used for instance). 
# In this case I want to keep only the 5 chromosomes of Arabidopsis thaliana
vcftools --gzvcf subset_80_biallelic_only_alt_no_indels_no_singletons.recode.vcf.gz \
			--chr Chr1 --chr Chr2 --chr Chr3 --chr Chr4 --chr Chr5 \
			--recode --out  subset_80_biallelic_only_alt_no_indels_no_singletons_only_chr

```

Get list of accessions in vcf file:

```
$ bcftool query -l  subset_80.recode.vcf.gz > order_accession.txt
$ cat order_accession.txt
1001
1002
1003

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
4. Import in R the output phenotype.assoc.clean.txt file to visualize GWAS results, see [gemma_analysis.Rmd](gemma_analysis.Rmd)

* The part 1, 2, and 4 are done interactively in R (need to be adjusted according to the dataframe used)
* The part 3 is done in bash through the run_gwas_gemma.sh script. The only variables being the input vcf file used (change path in the script file) and the phenotype file given as first argument in command  line:
* The part 4 is done interactively in R


```
bash run_gwas_gemma.sh phenotype.tsv vcf_file.vcf.gz
```

The output file is created in the subdirectory `output/`, which is automatically created by gemma from the directory containing the input vcf file. The name of the output file is `phenotype.assoc.clean.txt`. Note that the file `phenotype.assoc.txt` is the direct output of gemma but cannot be loaded by the `qqman` package in R due to a wrong organization of the column.


In addition, 1 log files are generated

One log file named `phenotype.log` is generated in the same directory `output/` and contains the different parameters of the run. Example of log output:

```
LOG FILE

Command: run_gwas_gemma.sh CHH_genes_cluster6_subset_64.tsv subset_64_accessions_only_alt_wo_singletons_biallelic_only_wo_indels_minDP3_minGQ25.recode.vcf.gz
Phenotype file: CHH_genes_cluster6_subset_64.tsv
VCF file: subset_64_accessions_only_alt_wo_singletons_biallelic_only_wo_indels_minDP3_minGQ25.recode.vcf.gz
Output file: /srv/biodata/dep_coupland/grp_hancock/johan/GWAS/dna_methylation/output/CHH_genes_cluster6_subset_64.assoc.clean.txt
Run finished on Fri Nov 16 13:16:58 CET 2018
Total time of the run: 13 seconds

Log output from GEMMA:
##
## GEMMA Version = 0.94
##
## Command Line Input = -bfile subset_64_accessions_only_alt_wo_singletons_biallelic_only_wo_indels_minDP3_minGQ25 -k /srv/biodata/dep_coupland/grp_hancock/johan/GWAS/dna_methylation/output/subset_64_accessions_only_alt_wo_singletons_biallelic_only_wo_indels_minDP3_minGQ25.cXX.txt -lmm 2 -o CHH_genes_cluster6_subset_64
##
## Summary Statistics:
## number of total individuals = 64
## number of analyzed individuals = 64
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs = 1566021
## number of analyzed SNPs = 17175
## REMLE log-likelihood in the null model = -114.18
## MLE log-likelihood in the null model = -114.937
## pve estimate in the null model = 0.848988
## se(pve) in the null model = 0.0639243
## vg estimate in the null model = 10.6163
## ve estimate in the null model = 0.73981
## beta estimate in the null model =   5.19004
## se(beta) =   0.107515
##
## Computation Time:
## total computation time = 0.204333 min
## computation time break down:
##      time on eigen-decomposition = 0 min
##      time on calculating UtX = 0.0015 min
##      time on optimization = 0.0626667 min
##

```


More details about the run is also displayed as standard output when running `run_gwas_gemma.sh` and can be redirected to a log file if needed.
```
bash run_gwas_gemma.sh phenotype.tsv vcf_file.vcf.gz > log.txt
```


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

This file can then be used as third argument in run_gwas_gemma.sh such as:

```
bash run_gwas_gemma.sh phenotype.tsv vcf_file.vcf.gz covariate_file.txt
```

TODO: implement covariate analysis in the script. Currently, one needs to proceed manually from the script and do the polishing step. Done but need to be tested.


# Analysis in R

The output file `phenotype.assoc.clean.txt` generated by gemma can be imported in R and plotted using qqman package. See the example in [gemma_analysis.Rmd](gemma_analysis.Rmd). 
An example of GWAS plot would look like this:

![](images/example_manhattan_plot.png)

In this case, no SNP has a p-value below the threshold of -log10(10E-5) (indicated by the blue line). To know more about the qqman package

* https://www.biorxiv.org/content/early/2014/05/14/005165.full.pdf+html
* https://cran.r-project.org/web/packages/qqman/qqman.pdf
* https://github.com/stephenturner/qqman



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


# Pipeline with ChromatinJ
If GWAS is to be performed in ChromatinJ output:
1. Generate a dataframe with all the variables and order the accessions so that they match VCF sample order with the R script process_chromatinJ_output.R (see example output file test_export_df.txt in directory 'examples')
2. Generate a file for one phenotype from the previously generated dataframe (test_export_df.txt) with the R script generate_phenotype_file.Rmd. The output is a file with a single column (no header) containing the phenotype of interest, scaled or not. See example files unscaled "Area_nucleus.tsv"
3. Process the phenotype file with plink and gemma (run_gwas_gemma.sh). This step should be run directly into the directory containing the vcf file and the phenotype file.
4. Import in R the output phenotype.assoc.clean.txt file to visualize GWAS results (script gemma_analysis.Rmd)


# Authors
* **Johan Zicola** - [johanzi](https://github.com/johanzi)

# License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details



