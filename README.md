# J-statistic

## Scope

J-statistic is an easy-to-use statistical pipeline designed to infer associations between different mutation types, including but not limited to single nucleotide polymorphisms (SNP) and copy number variations (CNV) observed from the microarray platform. This R package applies the J-statistic to test the following null hypotheses:

1. SNP differences outside of CNV regions follow complete spatial randomness.
2. SNP differences outside of CNV regions have similar properties in clearly defined CNV nearby regions as for regions further away.
3. SNP differences outside of CNV regions have similar properties everywhere on the chromosome.

Thus, the J-statistic approach aims to infer whether the existence of CNVs could influence the spacing of SNP differences outside of CNVs or whether there tend to be more or less SNP differences in regions near the CNVs compared to those farther away. J-statistic was originally tested use SNP calls from [placeholder] and CNV calls from [PennCNV](http://penncnv.openbioinformatics.org/en/latest/). 

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

The `J_statistic.R` will automatically attempt to install dependencies not found in the current working environment. 

To install Devtools, Intervals, and Optparse from CRAN manually, use the following command in an interactive R environment: 

```sh
install.packages(c("intervals", "devtools", "optparse"))  
```

To install Rowr from GitHub manually, use the following command in an interactive R environment:

```sh
library(devtools)
install_github("cvarrichio/rowr")
```

## Parameters

Summary of the different arguments that `J-statistic.R` uses as input. The script will output summary statistics for each chromosome of each individual microarray, along with Rainbow plots, Rainfall plots, and J function plots. Using the long-form arguments when running the `J-statistic.R` in a terminal is highly recommended. 

Short-Form Argument Name| Long-Form Argument Name| Argument Type | Argument Description | Argument Range 
--- | --- | --- | --- | ---
-s | --snp | Character | Absolute file path of SNP calls (CSV file) | User-specified file path
-c | --cnv | Character | Absolute file path of CNV calls (CSV file) | User-specified file path
-o | --output | Character | Absolute file path of output directory | User-specified file path
-r | --nrun | Integer | Number of bootstrap simulations |  Default:10 ; Recommended Range: 10-100
-h | --het | Integer | Minimum number of heterozygous SNPs |  Default:10 ; Recommended Range: 10-100
-x | --seed | Integer | Number of bootstrap simulations |  Default:12345
-m | --assocation_max_distance | Integer | Maximum distance between SNPs and CNVs to test association between SNPs and CNVs. | Default:10000000; Recommended Range: 1000000-50000000
-i | --association_interval_distance | Integer | Interval step size to test association between SNPs and CNVs. | Default: 5000; Recommended Range: 1000-10000
-n | --cluster_max_distance | Integer | Maximum distance between SNPs to test for existence of SNP clusters. | Default:100000; Recommended Range: 10000-500000
-j | --cluster_interval_distance | Integer | Interval step size to test for existence of SNP clusters. | Default: 5000; Recommended Range: 1000-10000
-a | --alpha | Float | Alpha value for significance threshold of statical tests. | Default: 0.05; Recommended Range: 0.01-0.10

The `--max_distance` and `--interval_distance` parameters set the maximum distance and step size to test for existence of SNP clustering and SNP-CNV association. Shown below is a visualization of different `--max_distance` and `--interval_distance` parameters. 

# <img src="https://github.com/HillLab/J-statistic/blob/main/J_statistic_Grid_Points_Figure.jpg" alt="Grid Points Figure" width="1000" style="float: right;"/>


## Input

The expected input into the `J-statistic.R` script is one SNP file in CSV format and one CNV file in CSV format. Examples of a correctly formatted SNP and CNV input file can be found in `./example/input`. See Section `Example Dataset and Tutorial` to run this example data. The required format, column names, and data fields for custom SNP and CNV input is described below.

### SNP Input File
The SNP input file requires a unique identifier for each SNP (ID), the chromosome where the SNP exists (SNP.Chromosome), the base position on the chromosome where the SNP exists (Position), and a binary matrix for each sample (e.g. SNP_mDIV_A1 and SNP_mDIV_C4). An unlimited number of samples can be included in one SNP input file, where SNPs associated with each sample is represented by one unique column in the SNP input file.

The binary matrix for each sample column (e.g. SNP_mDIV_A1 and SNP_mDIV_C4) uses the integer value of `1` to show that a SNP exists in the sample and the integer value of `0` to show that a SNP does not exist in the sample. In the example SNP input file below, the JAX00000002 SNP only exists in SNP_mDIV_A1, the JAX00000003 SNP only exists in SNP_mDIV_C4, and the JAX00240566 SNP does not exist in either SNP_mDIV_A1 or SNP_mDIV_C4.  

ID | SNP.Chromosome | Position | SNP_mDIV_A1 | SNP_mDIV_C4
--- | --- | --- | --- | ---
JAX00000002 | 1 | 30460970 | 1 | 0
JAX00000003 | 1 | 30461840 | 0 | 1
JAX00240566 | 1 | 30491060 | 0 | 0
... | ... | ... | ... | ...

### CNV Input File
The CNV input file requires the chromosome where the CNV exists (CNV.Chromosome), the start base position of the CNV segment (Start), the end base position of the CNV segment (End), and the sample name (Sample). An unlimited number of samples can be included in one CNV input file, where CNVs associated with each sample is represented by one or more rows in the CNV input file.

In the example CNV input file below, there exists a CNV on Chromosome 1 from 30688847-30934770 base position observed in Sample SNP_mDIV_A1. 

CNV.Chromosome | Start | End | Sample
--- | --- | --- | --- 
1 | 30688847 | 30934770 | SNP_mDIV_A1 
1 | 132977602 | 133459657 | SNP_mDIV_A1 
1 | 132977602 | 133459657 | SNP_mDIV_C4 
... | ... | ... | ...
 
## Output

Using SNP and CNV data from two samples (Sample 1 and Sample 2) as input, the expected structure of the `J-statistic.R` script output in the user-specified directory is shown.
<pre>
├── Output Directory                                                              // User-specified output directory for all test results.
│   └── summary_output_statistics.csv                                             // Summary file of statistical results for all chromosomes of all samples.
│   └── Sample_1.csv                                                              // Sample 1 merged table of input SNP and CNV data. 
│   └── Sample_2.csv                                                              // Sample 1 merged table of input SNP and CNV data. 
│   └── Sample_1                                                                  // Sample 1 directory of all chromosome-specific statistical results and plots. 
│   │   └── Sample_1_Chr_1                                                        // Sample 1 sub-directory of Chromosome 1-specific statistical results and plots.
│   │   │   └── Sample_1_Chr_1_Clustering.Rdata                                   // Sample 1 Chromosome 1: Rdata file of KS test for existence of SNP clusters.
│   │   │   └── Sample_1_Chr_1_Association_Test.Rdata                             // Sample 1 Chromosome 1: Rdata file of J-statistic test for SNP-CNV association.
│   │   │   └── Sample_1_Chr_1_J_statistic.pdf                                    // Sample 1 Chromosome 1: J-statistic plot.
│   │   │   └── Sample_1_Chr_1_Rainfall.pdf                                       // Sample 1 Chromosome 1: Rainfall plot.
│   │   │   └── Sample_1_Chr_1_Rainbow.pdf                                        // Sample 1 Chromosome 1: Rainbow plot.
│   │   │   └── Sample_1_Chr_1_stats.txt                                          // Sample 1 Chromosome 1: SNP cluster and SNP-CNV association test results.
│   │   └── Sample_1_Chr_2                                                        // Sample 1 sub-directory of Chromosome 2-specific statistical results and plots.
│   │   └── Sample_1_Chr_3                                                        // Sample 1 sub-directory of Chromosome 3-specific statistical results and plots.
│   │   └── ...
│   └── Sample_2                                                                  // Sample 2 directory of all chromosome-specific statistical results and plots. 
│   │   └── Sample_2_Chr_1                                                        // Sample 2 sub-directory of Chromosome 1-specific statistical results and plots.
│   │   └── Sample_2_Chr_2                                                        // Sample 2 sub-directory of Chromosome 2-specific statistical results and plots. 
│   │   └── Sample_2_Chr_3                                                        // Sample 2 sub-directory of Chromosome 3-specific statistical results and plots.
│   │   └── ...
</pre>

## Example Dataset and Tutorial

The `J_statistic.R` converts the input SNP and CNV data files into a single processed file, then tests for existence of SNP clusters and association between SNPs and CNVs, as well as outputs the statistical results for each sample as a summary Excel file with associated plots. 

The example input data found at `./example/input` includes one SNP data file (`./example/input/example_SNP.csv`) containing SNPs and one CNV data file (`./example/input/example_CNV.csv`) containing CNVs for 2 unique samples. 

`J_statistic.R` run using the following command in terminal will produce the output found at `./example/output`:

```sh
cd [Directory where J-statistic.R is located]

Rscript J_statistic.R --snp ./example/input/example_SNP.csv --cnv ./example/input/example_CNV.csv --output ./example/output
```

## Use Cases

## Structure of J-statistic package
<pre>
├── README.md                                                                     // README. 
├── LICENSE                                                                       // Copy of the Creative Commons Attribution 4.0 License (CC BY). 
├── requirements.txt                                                              // List of package dependencies.   
├── example                                                                       // Directory containing example input and output files. 
│   └── input                                                                     // Directory containing example input SNP and CNV files. 
│   │   └── example_CNV.csv                                                       // Example SNP file used as input into J-statistic script. 
│   │   └── example_SNP.csv                                                       // Example CNV file used as input into J-statistic script. 
│   └── output                                                                    // Directory containing example output summary statistics and Rainfall, Rainbow, and J-statistic plots. 
│   │   └── SNP_mDIV_A1.SNP09_319_111109                                          // Directory containing summary statistics and plots for SNP_mDIV_A1.SNP09_319_111109 sample. 
│   │   └── SNP_mDIV_A1.SNP09_319_111109.csv                                      // Processed data combining SNP and CNV input files. 
│   │   └── SNP_mDIV_C4_SNP09_300_102709                                          // Directory containing summary statistics and plots for SNP_mDIV_C4_SNP09_300_102709 sample. 
│   │   └── SNP_mDIV_C4_SNP09_300_102709.csv                                      // Processed data combining SNP and CNV input files. 
├── J_statistic.R                                                                 // J-statistic R script. 
├── reference_functions.R                                                         // Supporting functions needed for J-statistic R script. 
</pre>

## Citing J-statistic

Placeholder for citation.

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
