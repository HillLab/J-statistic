# J-statistic

## Scope

## Installation

J-statistic is freely available on [GitHub](https://github.com/HillLab/J-statistic). Installation requires [git](https://git-scm.com/) and [git lfs](https://git-lfs.github.com/) installed. 

Install J-statistic to your working directory using the follow command in a terminal.

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

## Tutorial

## Structure of J-statistic package

├── README.md                                                                     // README. \
├── LICENSE                                                                       // Copy of the Creative Commons Attribution 4.0 License (CC BY). \
├── requirements.txt                                                              // List of package dependencies.   \
├── data                                                                          // Directory containing 2 reference datasets. \
│   └── GGG.csv                                                                   // Reference dataset. 
│   └── SNP cluster detection - SNP CNV association test - FUNCTION v4.RData      // Reference dataset. 
├── example                                                                       // Directory containing example input and output files. 
│   └── input                                                                     // Directory containing example input SNP and CNV files. 
│   │.  └── example_CNV.csv                                                       // Example SNP file used as input into J-statistic script. 
│   │.  └── example_SNP.csv                                                       // Example CNV file used as input into J-statistic script. 
│   └── output                                                                    // Directory containing example output summary statistics and Rainfall, Rainbow, and J-statistic plots. 
│   │.  └── SNP_mDIV_A1.SNP09_319_111109                                          // Directory containing summary statistics and plots for one microarray. 
│   │.  └── SNP_mDIV_A1.SNP09_319_111109.csv                                      // Processed data combining SNP and CNV input files. 
├── J_statistic.R                                                                 // J-statistic R script. 

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
