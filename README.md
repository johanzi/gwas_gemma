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


**For instance:**

Get list of accessions in vcf file:
```
$ bcftool query -l myvcffile.vcf > order_accession.txt
$ cat order_accession.txt
1001
1002
1003
```

Content of the phenotype file:
```
$ cat phenotype.tsv 
12.3
13.4
15.3
```
The value 12.3, 13.4, 15.3 being the height of the accessions 1001, 1002, and 1003, respectively.

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



