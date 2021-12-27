# J-statistic

## Scope

J-statistic is an easy-to-use statistical pipeline designed to infer associations between different mutation types, including but not limited to single nucleotide polymorphisms (SNP) and copy number variations (CNV) observed from the microarray platform. This R package applies the J-statistic to test the following null hypotheses:

1. SNP differences outside of CNV regions follow complete spatial randomness.
2. SNP differences outside of CNV regions have similar properties in clearly defined CNV nearby regions as for regions further away.
3. SNP differences outside of CNV regions have similar properties everywhere on the chromosome.

Thus, the J-statistic approach aims to infer whether the existence of CNVs could influence the spacing of SNP differences outside of CNVs or whether there tend to be more or less SNP differences in regions near the CNVs compared to those farther away. J-statistic was originally tested using microarray probe data.

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
-s | --snp_file | Character | Absolute file path of SNP calls (CSV file) | User-specified file path
-c | --cnv_file | Character | Absolute file path of CNV calls (CSV file) | User-specified file path
-o | --output_dir | Character | Absolute file path of output directory | User-specified file path
-r | --nrun | Integer | Number of bootstrap simulations |  Default:1000 ; Recommended Range: 1000-10000
-h | --min_snp | Integer | Minimum number of SNPs |  Default:10 ; Recommended Range: 10-100
-z | --seed | Integer | Number of bootstrap simulations |  Default:12345
-m | --assocation_max_distance | Integer | Maximum distance between SNPs and CNVs to test association between SNPs and CNVs. | Default:10000000; Recommended Range: 1000000-50000000
-i | --association_interval_distance | Integer | Interval step size to test association between SNPs and CNVs. | Default: 5000; Recommended Range: 1000-10000
-n | --cluster_max_distance | Integer | Maximum distance between SNPs to test for existence of SNP clusters. | Default:100000; Recommended Range: 10000-500000
-j | --cluster_interval_distance | Integer | Interval step size to test for existence of SNP clusters. | Default: 5000; Recommended Range: 1000-10000
-a | --alpha | Float | Alpha value for significance threshold of statical tests. | Default: 0.05; Recommended Range: 0.01-0.10
-w | --wgs_file | Character | Absolute file path of file with start and stop coordinates for all segments in whole genome/exome (CSV file) | User-specified file path (Mandatory argument for WGS/WES input data; Null argument for microarray probe input data)
-x | --wgs_nsample | Integer | Number of randomly sampled locations in whole genome/exome segments for mutation position null distribution estimate. | Default: 500000; Recommended Range: 100000-1000000 for WGS and 10000-100000 for WES

The `--max_distance` and `--interval_distance` parameters set the maximum distance and step size to test for existence of SNP clustering and SNP-CNV association. Shown below is a visualization of different `--max_distance` and `--interval_distance` parameters. 

# <img src="https://github.com/HillLab/J-statistic/blob/main/example/J_statistic_Grid_Points_Figure.jpg" alt="Grid Points Figure" width="1000" style="float: right;"/>

The `--wgs_file` and `--wgs_nsample` parameters must only be set when using whole genome or whole exome sequencing input data. See `./example/input/example_exome.csv` and `./example/input/example_genome.csv` for an example of the format and required fields of the `--wgs_file` file needed.

## Input

The expected input into the `J-statistic.R` script is one SNP file in CSV format and one CNV file in CSV format. Examples of a correctly formatted SNP and CNV input file can be found in `./example/input`. See Section `Example Dataset and Tutorial` to run this example data. The required format, column names, and data fields for custom SNP and CNV input is described below.

### SNP Input File 

##### Parameter: `--snp_file`

The SNP input file requires the chromosome where the SNP exists (SNP.Chromosome), the base position on the chromosome where the SNP exists (Position), and a binary matrix for each sample (e.g. SNP_mDIV_A1 and SNP_mDIV_C4). An unlimited number of samples can be included in one SNP input file, where SNPs associated with each sample is represented by one unique column in the SNP input file. 

The binary matrix for each sample column (e.g. SNP_mDIV_A1 and SNP_mDIV_C4) uses the integer value of `1` to show that a SNP exists in the sample and the integer value of `0` to show that a SNP does not exist in the sample. Some SNPs may not exist (denoted by `0`) in any of the samples. In the example SNP input file below, the Chr. 1 (30460970) SNP only exists in SNP_mDIV_A1, the Chr. 1 (30461840) SNP only exists in SNP_mDIV_C4, and the Chr. 1 (30491060) SNP does not exist in either SNP_mDIV_A1 or SNP_mDIV_C4.  

SNP.Chromosome | Position | SNP_mDIV_A1 | SNP_mDIV_C4
--- | --- | --- | ---
1 | 30460970 | 1 | 0
1 | 30461840 | 0 | 1
1 | 30491060 | 0 | 0
... | ... | ... | ... | ...

### CNV Input File

##### Parameter: `--cnv_file`

The CNV input file requires the chromosome where the CNV exists (CNV.Chromosome), the start base position of the CNV segment (Start), the end base position of the CNV segment (End), and the sample name (Sample). An unlimited number of samples can be included in one CNV input file, where CNVs associated with each sample is represented by one or more rows in the CNV input file.

In the example CNV input file below, one CNV on Chromosome 1 from 30688847-30934770 base position is observed in Sample SNP_mDIV_A1. 

CNV.Chromosome | Start | End | Sample
--- | --- | --- | --- 
1 | 30688847 | 30934770 | SNP_mDIV_A1 
1 | 132977602 | 133459657 | SNP_mDIV_A1 
1 | 132977602 | 133459657 | SNP_mDIV_C4 
... | ... | ... | ...
 
### Optional: WGS/WES Input Data

When using WGS/WES data as input into the `J-statistic.R` script, three changes to the parameters and input data must be made. 

##### Parameter: `--snp_file`

For WGS/WES input data, the SNP input file must contain only single base mutations that exist (denoted by `1`) in at least one of the samples, unlike the SNP input data for microarray probes. Thus, each mutation (row) in the SNP input file must have at least one `1` value in at least one of the samples (column), see below for an example SNP input file. 

SNP.Chromosome | Position | SNP_mDIV_A1 | SNP_mDIV_C4
--- | --- | --- | ---
1 | 30460970 | 1 | 0
1 | 30461840 | 0 | 1
1 | 30491060 | 1 | 1
... | ... | ... | ... | ...

##### Parameter: `--wgs_file `

For WGS/WES input data, a file that describes the start (Start) and stop (Stop) base pair locations as well as chromosome (Chromosome) of each possible segment where single base mutations may be observed is required. See `./example/input/example_exome.csv` and `./example/input/example_genome.csv` for an example of the format and required fields of the `--wgs_file` file needed. In the example WGS/WES input file below, 

Chromosome | Start | End
--- | --- | ---
1 | 69090 | 70008
1 | 450739 | 451678
1 | 685715 | 686654
... | ... | ...

##### Parameter: `--wgs_nsample`

For WGS/WES input data, possible base pair locations in whole genome/exome segments are randomly sampled to generate a null distribution of single base mutation positions to test for SNP cluster existence and SNP-CNV association. Generally, a minimum random sample of 100000 base pair locations for whole genome input data or a minimum random sample of 10000 base pair locations for whole exome input data is recommended.

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

## Example Use Cases

Example Use Case 1. Test for SNP-CNV association using smaller interval sizes than default. Low number of bootstrap simulations (`--nrun 10`) to speed up run for proof of concept. 

```sh
Rscript J_statistic.R --snp_file ./example/input/example_SNP.csv --cnv_file ./example/input/example_CNV.csv --output_dir ./example/example_use_case_1 --association_max_distance 10000000 --association_interval_distance 1000 --nrun 10
```

Example Use Case 2. Test for SNP cluster existence using smaller maximum distances between SNPs and interval sizes than default. Low number of bootstrap simulations (`--nrun 10`) to speed up run for proof of concept. 

```sh
Rscript J_statistic.R --snp_file ./example/input/example_SNP.csv --cnv_file ./example/input/example_CNV.csv --output_dir ./example/example_use_case_2 --cluster_max_distance 50000 —-cluster_interval_distance 1000 --nrun 10
```

Example Use Case 3. Test for SNP cluster existence and SNP-CNV association using a larger alpha value than default. Low number of bootstrap simulations (`--nrun 10`) to speed up run for proof of concept. 

```sh
Rscript J_statistic.R --snp_file ./example/input/example_SNP.csv --cnv_file ./example/input/example_CNV.csv --output_dir ./example/example_use_case_3 --alpha 0.10 --nrun 10
```

Example Use Case 4. Test for SNP cluster existence and SNP-CNV association using WGS single base mutation data instead of SNP data from microarray platforms (default). 

```sh
Rscript J_statistic.R --snp_file ./example/input/example_wgs.csv --cnv_file ./example/input/example_CNV.csv --output_dir ./example/example_use_case_4 --wgs_file ./example/input/example_genome.csv --wgs_nsample 500000
```

The output data for Example Use Case 1-3 can be found open-access on Zenodo: <a href="https://doi.org/10.5281/zenodo.5804599"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5804599.svg" alt="DOI"></a>

## Structure of J-statistic package
<pre>
├── README.md                                                                     // README. 
├── LICENSE                                                                       // Copy of the Creative Commons Attribution 4.0 License (CC BY). 
├── requirements.txt                                                              // List of package dependencies.   
├── example                                                                       // Directory containing example input and output files. 
│   └── input                                                                     // Directory containing example input SNP and CNV files. 
│   │   └── example_CNV.csv                                                       // Example SNP file used as input into J-statistic script. 
│   │   └── example_SNP.csv                                                       // Example CNV file used as input into J-statistic script. 
│   │   └── example_exome.csv                                                     // Chromosome and base pair start/end locations of all human exome segments. 
│   │   └── example_genome.csv                                                    // Chromosome and base pair start/end locations of all human chromosomes. 
│   └── output                                                                    // Directory containing example output summary statistics and Rainfall, Rainbow, and J-statistic plots. 
│   │   └── SNP_mDIV_A1.SNP09_319_111109                                          // Directory containing summary statistics and plots for SNP_mDIV_A1.SNP09_319_111109 sample. 
│   │   └── SNP_mDIV_A1.SNP09_319_111109.csv                                      // Processed data combining SNP and CNV input files. 
│   │   └── SNP_mDIV_C4_SNP09_300_102709                                          // Directory containing summary statistics and plots for SNP_mDIV_C4_SNP09_300_102709 sample. 
│   │   └── SNP_mDIV_C4_SNP09_300_102709.csv                                      // Processed data combining SNP and CNV input files. 
├── J_statistic.R                                                                 // J-statistic R script. 
├── reference_functions.Rdata                                                     // Supporting functions needed for J-statistic R script. 
</pre>

## Citing J-statistic

Placeholder for citation.

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
