# J-statistic

## Scope

J-statistic is an easy-to-use statistical pipeline designed to infer associations between different mutation types, including but not limited to single nucleotide polymorphisms (SNP) and copy number variations (CNV) observed from the microarray platform. This R package applies the J-statistic to test the following null hypotheses:

1. SNP differences outside of CNV regions follow complete spatial randomness.
2. SNP differences outside of CNV regions have similar properties in clearly defined CNV nearby regions as for regions further away.
3. SNP differences outside of CNV regions have similar properties everywhere on the chromosome.

Thus, the J-statistic approach aims to infer whether the existence of CNVs could influence the spacing of SNP differences outside of CNVs or whether there tend to be more or less SNP differences in regions near the CNVs compared to those farther away. 

## Installation

J-statistic is freely available on [GitHub](https://github.com/HillLab/J-statistic). Installation requires [git](https://git-scm.com/) and [git lfs](https://git-lfs.github.com/) installed. 

Install J-statistic to your working directory using the following command in a terminal.

```sh
git clone https://github.com/HillLab/J-statistic
```

## Requirements

J-statistic has the following dependencies:
- R version 4.1.2
- [Devtools 2.4.3](https://cran.r-project.org/web/packages/devtools/index.html)
- [Intervals 0.15.2](https://cran.r-project.org/web/packages/intervals/index.html)
- [Optparse 1.7.1](https://cran.r-project.org/web/packages/optparse/index.html)
- [Rowr](https://github.com/cvarrichio/rowr)

## Usage

Summary of the different arguments that J-statistic.R uses as input. The script will output summary statistics for each chromosome of each individual microarray, along with Rainbow plots, Rainfall plots, and J function plots.

Short-Form Argument Name| Long-Form Argument Name| Argument Type | Argument Description | Argument Range 
--- | --- | --- | --- | ---
-s | --snp | Character | Absolute file path of SNP calls (CSV file) |
-c | --cnv | Character | Absolute file path of CNV calls (CSV file) |  
-o | --output | Character | Absolute file path of output directory |  
-n | --nrun | Integer | Number of bootstrap simulations |  Default:10 ; Recommended Range: 10-100
-h | --het | Integer | Minimum number of heterozygous SNPs |  Default:10 ; Recommended Range: 10-100
-x | --seed | Integer | Number of bootstrap simulations |  Default:12345

## Tutorial

The `J_statistic.R` converts the input SNP and CNV data files into a single processed file, then tests for statistical association between SNPs and CNVs and output the statistical results as a summary Excel file with associated plots. 

Using the example input data found at `./example/input`, `J_statistic.R` run using the following command in terminal will produce the output found at `./example/output':

`cd [Directory where J-statistic.R is located]`

`Rscript J_statistic.R --snp ./example/input/example_SNP.csv --cnv ./example/input/example_CNV.csv --output ./example/output --nrun 10 --het 10 --seed 1`

## Structure of J-statistic package
<pre>
├── README.md                                                                     // README. 
├── LICENSE                                                                       // Copy of the Creative Commons Attribution 4.0 License (CC BY). 
├── requirements.txt                                                              // List of package dependencies.   
├── data                                                                          // Directory containing 2 reference datasets. 
│   └── GGG.csv                                                                   // Reference dataset. 
│   └── SNP cluster detection - SNP CNV association test - FUNCTION v4.RData      // Reference dataset. 
├── example                                                                       // Directory containing example input and output files. 
│   └── input                                                                     // Directory containing example input SNP and CNV files. 
│   │   └── example_CNV.csv                                                       // Example SNP file used as input into J-statistic script. 
│   │   └── example_SNP.csv                                                       // Example CNV file used as input into J-statistic script. 
│   └── output                                                                    // Directory containing example output summary statistics and Rainfall, Rainbow, and J function plots. 
│   │   └── SNP_mDIV_A1.SNP09_319_111109                                          // Directory containing summary statistics and plots for one microarray. 
│   │   └── SNP_mDIV_A1.SNP09_319_111109.csv                                      // Processed data combining SNP and CNV input files. 
├── J_statistic.R                                                                 // J-statistic R script. 
</pre>

## Citing J-statistic

Placeholder for citation.

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
