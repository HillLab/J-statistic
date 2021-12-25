###################
#Load dependencies# 
###################

# Install libraries from CRAN and GitHub if not already installed
cran_packages <- c("intervals", "devtools", "optparse")
install.packages(setdiff(cran_packages, rownames(installed.packages())))  

library(devtools)
library(intervals)
library(optparse)

github_packages <- setdiff("rowr", rownames(installed.packages()))
if (length(github_packages) > 0) {
  install_github("cvarrichio/rowr")
}

library(rowr)

# Load reference functions
load(file = "reference_functions.RData")

######################################
#Stage 1: Setup user-input parameters#
######################################

option_list = list(
  make_option(c("-s", "--snp"), type="character", default=NULL, 
              help="Absolute file path of SNP calls (CSV file).", metavar="SNP_CALL_FILE_PATH"),
  make_option(c("-c", "--cnv"), type="character", default=NULL, 
              help="Absolute file path of CNV calls (CSV file).", metavar="CNV_CALL_FILE_PATH"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Absolute file path of output directory.", metavar="OUTPUT_DIRECTORY_FILE_PATH"),
  make_option(c("-r", "--n_run"), type="integer", default=10, 
              help="Number of bootstrap simulations/", metavar="NUMBER_OF_RUNS"),
  make_option(c("-h", "--min_snp"), type="integer", default=1000, 
              help="Minimum number of SNPs.", metavar="MIN_NUMBER_OF_SNPS"),
  make_option(c("-x", "--seed"), type="integer", default=12345, 
              help="Random seed.", metavar="SEED"),
  make_option(c("-m", "--association_max_distance"), type="integer", default=10000000, 
            help="Maximum distance to test association between SNPs and CNVs.", metavar="SNP-CNV_ASSOCIATION_MAX_DISTANCE"),
  make_option(c("-i", "--association_interval_distance"), type="integer", default=5000, 
            help="Interval step size to test association between SNPs and CNVs.", metavar="SNP-CNV_ASSOCIATION_INTERVAL_DISTANCE"),
  make_option(c("-n", "--cluster_max_distance"), type="integer", default=100000, 
              help="Maximum distance between SNPs to test for existence of SNP clusters.", metavar="SNP_CLUSTER_MAX_DISTANCE"),
  make_option(c("-j", "--cluster_interval_distance"), type="integer", default=5000, 
              help="Interval distance for SNP clustering test.", metavar="SNP_CLUSTER_INTERVAL_DISTANCE"),
  make_option(c("-a", "--alpha"), type="numeric", default=0.05, 
              help="Alpha value for significance threshold.", metavar="ALPHA")
)

# Parse user-specified parameters in terminal as vector
opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE);
opt = parse_args(opt_parser)

# CHECK: Mandatory SNP file path argument
if(is.null(opt$snp)) { 
  stop("Missing: Absolute file path of SNP calls (CSV file).")
}

# CHECK: Mandatory CNV file path argument
if(is.null(opt$cnv)) { 
  stop("Missing: Absolute file path of CNV calls (CSV file).")
}

# CHECK: Mandatory Output Directory argument
if(is.null(opt$output)) { 
  stop("Missing: Absolute file path of output directory.")
}

cat("Stage 1: Setup user-input parameters.\n")

# Set up all user-specified parameters
calls = read.csv(opt$snp,check.names=FALSE)
cnvsData = read.csv(opt$cnv,check.names=FALSE)
output_directory = opt$output
numRuns = opt$n_run
heterozygousCallCutoff = opt$min_snp
seed = opt$seed
max_distance = opt$association_max_distance
interval_distance = opt$association_interval_distance
cluster_max_distance = opt$cluster_max_distance
cluster_interval_distance = opt$cluster_interval_distance
alpha = opt$alpha

########################################################
#Stage 2: Pre-process input data for J-statistic script#
########################################################

cat("Stage 2: Pre-process input data for J-statistic script.\n")

# Create output folder
dir.create(output_directory, showWarnings = FALSE) 

# Initialize Stage 2 Progress Bar
stage_2_pb <- txtProgressBar(min = 0, max = length(4:ncol(calls)), style = 3)
stage_2_pb_counter = 0

# Output one processed file (contains both SNPs and CNVs) for each sample
for (i in 4:ncol(calls)){

  snpData <- calls[,c(1,2,3,i)]
  
  # Remove all SNPs found in X, Y, and M chromosomes from downstream statistical analysis
  snpData <- snpData[grep("X|Y|M",snpData$'SNP.Chromosome', invert=TRUE),]
  sampleName <- colnames(snpData)[4]

  # Remove all CNVs found in X, Y, and M chromosomes from downstream statistical analysis
  cnvData <- cnvsData[grep(sampleName,cnvsData$'Sample'),]
  cnvData <- cnvData[grep("X|Y|M",cnvData$'CNV.Chromosome', invert=TRUE),]
  
  # Skip output of processed file if sample has no CNVs
  if (nrow(cnvData)==0){
    write.csv(NA, file.path(output_directory, paste0(sampleName,"-NULL", ".csv")), row.names=FALSE, na='' )
    next
  }
  
  # Sort SNP dataframe by chromosome number and mutation position
  snpData <- snpData[order(as.numeric(as.character(snpData$'SNP.Chromosome')),snpData$'Position'),]
  
  # Sort CNV dataframe by chromosome number and segment start position
  cnvData <- cnvData[order(cnvData$'CNV.Chromosome',cnvData$'Start'),]
  
  # Merge SNP and CNV dataframe into a single dataframe
  snpData$nullID <- seq.int(nrow(snpData))
  cnvData$nullID <- seq.int(nrow(cnvData))
  outputtable<-merge(snpData,cnvData, by ="nullID", all.x=T, check.names=F)
  outputtable$'nullID'<-NULL

  # Output merged dataframe as CSV file
  write.csv(outputtable, file.path(output_directory, paste0(sampleName,".csv")), row.names=FALSE, na='' )
  
  # Update Stage 2 progress bar 
  stage_2_pb_counter = stage_2_pb_counter + 1
  setTxtProgressBar(stage_2_pb, stage_2_pb_counter)
    
}
close(stage_2_pb)

############################################
#Stage 3: J-statistic calculation and plots#
############################################

cat("Stage 3: J-statistic calculation and plots. \n")

# Initialize dataframe for the Excel document output with the generated statistics for all samples
excel_data<-data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c('Chromosome', 'p.KS', 'p.J', 'Num_CNV', 'SNP_CNV_Association', 'Sample'))))

# Initialize Stage 3 Progress Bar
stage_3_pb <- txtProgressBar(min = 0, max = length(list.files(path = output_directory, pattern= paste0('^(SNP_).*(.csv)$'), full.names=TRUE)), style = 3)
stage_3_pb_counter = 0

# Reads processed input files in folder in a loop
# Generates Rainbow, Rainfall, and J-statistic plots along with summary statistics to test for the existence of SNP clusters and association between SNPs and CNVs 
for (file in list.files(path=output_directory, pattern=".csv", all.files=TRUE, full.names=TRUE)){

  # Read in merged SNP and CNV dataframe
  cancer.data <- read.csv(paste0(file), header=T)
  
  # Initialize dataframe for the Excel document output with the generated statistics for one sample
  excel_table<-data.frame(matrix(ncol=6,nrow=length(as.numeric(na.omit(unique(cancer.data[["SNP.Chromosome"]]))))))
  colnames(excel_table)<- c('Chromosome', 'p.KS', 'p.J', 'Num_CNV', 'SNP_CNV_Association', 'Sample')
  
  # Get sample name
  filename<-gsub('.csv','',file)
  filenamesplit <- strsplit(filename, "/")[[1]]
  fileNoSuffix <- filenamesplit[length(filenamesplit)]
  
  # Create empty output directory
  dir.create(paste0(output_directory, '/', fileNoSuffix), showWarnings = F)

  # Generate statistics and plots for each chromosome with SNPs (above cutoff) and CNVs
  for (i in as.numeric(na.omit(unique(cancer.data[["SNP.Chromosome"]])))){

    #########################
    #Stage 3.1 Rainfall Plot#
    #########################
    
    # Subset merged dataframe for one chromosome
    cancer.data.19 <- subset(cancer.data, SNP.Chromosome == i)
    
    # Get the set of all SNPs (marked with either 0 or 1)
    probe.set.19 <- cancer.data.19$Position
    
    # Get the subset of SNPs of interest (marked with only 1)
    SNP.example <- sort(as.numeric(na.omit(probe.set.19[which(cancer.data.19[,fileNoSuffix] == 1)])))

    # If chromosome has no CNVs, then just make empty folder and note on output summary statistic table 
    if (!(i %in% (unique(cancer.data$CNV.Chromosome)))){
      dir.create(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_Chr_', i,'-NoCNVs'), showWarnings = FALSE) 
      excel_table[i,]<-c(i, "N", "N", 0, "No_CNVs",  fileNoSuffix)
      next
    }
    
    # Number of CNVs for each chromosome
    numberCNVs<- sum((cancer.data$CNV.Chromosome[!is.na(cancer.data$CNV.Chromosome)])==i)
    
    # If less than threshold cutoff of SNPs, skip the statistical analyses and plotting for the chromosome
    if (length(SNP.example)<= heterozygousCallCutoff){
      dir.create(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_Chr_', i,'-LessThanHetCutoff'), showWarnings = FALSE) 
      excel_table[i,]<-c(i, "<", "<", numberCNVs, "< Min_SNP_Cutoff", fileNoSuffix)
      next
    }

    dir.create(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_Chr_', i), showWarnings = FALSE) 

    # Save Rainfall Plot as PDF for each chromosome
    pdf(file=paste0(output_directory,'/',fileNoSuffix,'/',fileNoSuffix,'_Chr_', i,'/',fileNoSuffix,'_Chr_', i,'_Rainfall_Plot.pdf'))
    plot(log10(diff(SNP.example)) ~ SNP.example[-1], pch = 16, ylab="", xlab="")
    
    # Rainfall Plot Label on X-axis
    mtext(1, text = bquote("Chromosome Position (bp)"), cex = 1.2, line = 3)
    
    # Rainfall Plot Label on Y-axis
    mtext(2, text = "Distance between SNPs", cex = 1.2, line = 2.5)
    
    # Rainfall Plot Title
    title(main = paste0(fileNoSuffix,"_Chr_", i))
    
    dev.off()
    
    ###########################################################
    #Stage 3.2 Rainbow Plot and Test for SNP Cluster Existence#
    ###########################################################
    
    # Set up sequence of grid points for SNP cluster test
    grid.points <- (1: (cluster_max_distance/cluster_interval_distance)) * cluster_interval_distance

    # Directory for the output plots and statistics file to be saved
    save.directory.SNP <- paste0(output_directory,'/',fileNoSuffix,'/',fileNoSuffix,'_Chr_', i, '/')
   
    # SNP Cluster RData file name 
    file.name.SNP <- paste0(fileNoSuffix, "_Chr_",i,"_SNP_Cluster")

    # Set seed for reproducibility
    set.seed(seed)
    
    # Test for existence of SNP cluster using Kolmogorov-Smirnov test
    test.result.SNP <- SNP.cluster.test(detected.SNP = SNP.example, grid.d = grid.points, probe.set = probe.set.19, n.null = numRuns, na.rm = T, alpha = alpha)

    # Save the object test.result.SNP as the SNP Cluster RData file
    save(test.result.SNP, file = paste(save.directory.SNP, file.name.SNP, ".RData", sep = ""))

    # Load the SNP Cluster RData file
    load(file = paste(save.directory.SNP, file.name.SNP, ".RData", sep = ""))

    # Hypothesis test p-value of SNP cluster existence using KS test
    test.result.SNP$test.result

    # Subset merged dataframe for chromosome-specific CNVs
    sample.CNVs <- cancer.data
    sample.CNVs.subset <- subset(sample.CNVs, CNV.Chromosome == i)
    sample.CNVs.start <- sample.CNVs.subset$Start
    sample.CNVs.end <-sample.CNVs.subset$End
    sample.CNVs.good <- matrix(c(sample.CNVs.start, sample.CNVs.end), ncol=2)
    sample.CNVs.2 <- Intervals(sample.CNVs.good, closed = c( TRUE, TRUE ), type = "Z")

    # Subset merged dataframe for chromosome-specific SNPs of interest (only 1) for hypothesis test
    sample.SNPs.subset <- cancer.data[((cancer.data[,'SNP.Chromosome'] == i) & (cancer.data[,fileNoSuffix] == 1)),]
    sample.SNPs <- as.numeric(na.omit(sample.SNPs.subset$Position))

    # Subset merged dataframe for all chromosome-specific SNPs (0 and 1) to form null distribution 
    sample.probe.subset <- cancer.data[(cancer.data[,'SNP.Chromosome'] == i),]
    sample.probe <- as.numeric(na.omit(sample.probe.subset$Position))
    
    # Remove nans from CNVs if present 
    sample.CNVs.2<-sample.CNVs.2[!is.na(sample.CNVs.2)]

    # Rainbow Plot
    SNP.CNV.dist <- distance_to_nearest(sample.SNPs, sample.CNVs.2)
    
    pdf(file=paste0(output_directory, "/", fileNoSuffix,'/',fileNoSuffix,'_Chr_', i, '/',fileNoSuffix, '_Chr_', i, '_Rainbow_Plot.pdf'))
    
    plot(sample.CNVs.2, axes = FALSE, y = rep(0.1, dim(sample.CNVs.2)[1]),
         xlim = range(as.numeric(na.omit(sample.SNPs))), 
         ylim = c(1, max(as.numeric(na.omit((ceiling(log10(SNP.CNV.dist))))))), col = 2, cex = 1.25)
    
    points(log10(SNP.CNV.dist) ~ sample.SNPs)
    abline(v = apply(sample.CNVs.2, 1, mean), lty = 2, col =3, lwd = 1.5)
    
    # Rainbow Plot X-axis
    axis(1, at = 2*0:10*max_distance, label = as.character(2*0:10), font = 1)  
    
    # Rainbow Plot Y-axis
    axis(2, at = 0:max(as.numeric(na.omit((ceiling(log10(SNP.CNV.dist)))))), labels = c(1, sapply(1:max(as.numeric(na.omit((ceiling(log10(SNP.CNV.dist)))))), function(i) as.expression(bquote(10^ .(i))))))  # y-axis
    
    # Rainbow Plot Label on X-axis
    mtext(1, text = paste("Chromosome Position (", max_distance, "bp)", sep=""), cex = 1.2, line = 3)
    
    # Rainbow Plot Label on Y-axis
    mtext(2, text = "Distance from SNPs to closest CNV (bp)", cex = 1.2, line = 2.5)
    
    # Rainbow Plot Title
    title(main = paste0(fileNoSuffix,"_Chr_", i))

    dev.off()
    
    ############################################################
    #Stage 3.3 J Function Plot and Test for SNP-CNV Association#
    ############################################################
    
    # Set up sequence of grid points for SNP-CNV association test
    r_use_1 <- 1 : (max_distance/interval_distance) * interval_distance

    # Directory for the output plots and statistics file to be saved
    save.directory.SNP.CNV <- paste0(output_directory,'/',fileNoSuffix,'/',fileNoSuffix,'_Chr_', i, '/')
    
    # SNP-CNV Association RData file name 
    file.name.SNP.CNV <- paste0(fileNoSuffix,"_Chr_",i, "_SNP-CNV_Association")

    # Set seed for reproducibility
    set.seed(seed)
    test.result.SNP.CNV <- SNP.CNV.test(sample.SNPs = sample.SNPs, sample.probe = sample.probe, sample.CNVs = sample.CNVs.2, r_p = r_use_1, n.rep = numRuns, alpha = alpha)

    # Save the object test.result.SNP.CNV as the SNP-CNV Association RData file
    save(test.result.SNP.CNV, file = paste(save.directory.SNP.CNV, file.name.SNP.CNV, ".RData", sep = ""))
    
    # Load the SNP-CNV Association RData file
    load(, file = paste(save.directory.SNP.CNV, file.name.SNP.CNV, ".RData", sep = ""))

    # Negative association, larger value, tend to be farther away
    # Positive association, smaller value, tend to be close together

    # J-function Plot
    pdf(file=paste0(output_directory, '/',fileNoSuffix,'/',fileNoSuffix,'_Chr_', i, '/', fileNoSuffix,'_Chr_',i,'_J-statistic_Plot.pdf'))
    plot.upper <- max(c(test.result.SNP.CNV$sample.info[["J.fun.sample"]], test.result.SNP.CNV$Null.output[["gcb.upper"]]), na.rm =  T) + 0.01
    plot.lower <- min(c(test.result.SNP.CNV$sample.info[["J.fun.sample"]], test.result.SNP.CNV$Null.output[["gcb.lower"]]), na.rm = T) - 0.01
    plot(test.result.SNP.CNV$sample.info[["J.fun.sample"]] ~ r_use_1, type = 'l', ylim = c(plot.lower, plot.upper), lwd = 2, xlab = 'r', ylab = 'J function', main = "J function")
    lines(test.result.SNP.CNV$Null.output[["gcb.upper"]] ~ r_use_1, lty = 2, lwd = 2, col = 4)
    lines(test.result.SNP.CNV$Null.output[["gcb.lower"]] ~ r_use_1, lty = 2, lwd = 2, col = 4)
    legend("topleft", c("Observed J Function", "Global Confidence Bands"), lty = c(1,2), lwd = c(2,2), col = c(1, 4), cex=0.5)

    dev.off()

    # J-statistic association symbols
    # +- : J function crossed both high and low global confidence bands
    # +  : J function crossed low global confidence band
    # -  : J function crossed high global confidence band
    # /  : J function did not cross either high and low global confidence bands
    graphLinePoints <- test.result.SNP.CNV$sample.info[["J.fun.sample"]]
    graphUpperPoints <- test.result.SNP.CNV$Null.output[["gcb.upper"]]
    graphLowerPoints <- test.result.SNP.CNV$Null.output[["gcb.lower"]]
    
    if ((any(graphLinePoints[-1:-100] >graphUpperPoints[-1:-100])) &
        (any(graphLinePoints[-1:-100] <graphLowerPoints[-1:-100]))){
      association<-'+-'

    } else if (any(graphLinePoints[-1:-100] <graphLowerPoints[-1:-100])){
      association<-'+'
    } else if (any(graphLinePoints[-1:-100] >graphUpperPoints[-1:-100])){
      association<-'-'
    } else {
      association<-'/'  # \ is special character, need 2 to produce 1 in end 
    }
    
    #if unadjusted P-value of J-statistic is greater than alpha, no association between SNPs and CNVs
    if ((as.numeric(test.result.SNP.CNV$test.result)[1])>alpha){
      association<-'/'
    }

    # Output statistical results of KS (SNP cluster existence) and J-statistic (SNP-CNV association) tests to file, one output text file per chromosome per sample 
    chars <- capture.output(cat(fileNoSuffix,'Chromosome',i,sep='\n'),print(test.result.SNP$test.result), print(test.result.SNP.CNV$test.result),
                            cat('numCNVs',numberCNVs, sep='\n'), cat('Association',association,sep='\n'))    
    writeLines(chars, con = file(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_Chr_', i, '/', fileNoSuffix,'_Chr_',i,'_Stats.txt')))
    
    # Add statistical results to summary statistics dataframe
    excel_table[i,]<-c(i, (as.numeric(test.result.SNP$test.result)[1]), (as.numeric(test.result.SNP.CNV$test.result)[1]), numberCNVs,association, fileNoSuffix)
    closeAllConnections()
    
  }

  # Merge summary statistics dataframe for each sample into one combined dataframe for all samples  
  colnames(excel_table)<- c('Chromosome', 'p.KS', 'p.J', 'Num_CNV', 'SNP_CNV_Association', 'Sample')
  rownames(excel_table)<-NULL
  excel_data <- rbind(excel_data,excel_table)
  
  # Update Stage 3 progress bar 
  stage_3_pb_counter = stage_3_pb_counter + 1
  setTxtProgressBar(stage_3_pb, stage_3_pb_counter)
}

close(stage_3_pb)

# Update column names for summary statistics dataframe of all samples
colnames(excel_data)<- rep(c('Chromosome', 'p.KS', 'p.J', 'Num_CNV', 'SNP_CNV_Association', 'Sample'))

# Write summary statistics dataframe of all samples to CSV file in user-specified output directory
write.table(excel_data, file=paste(output_directory, 'summary_output_statistics.csv', sep="/"),row.names=FALSE,col.names=T,sep=",",na='')

# Close all connections
closeAllConnections()

