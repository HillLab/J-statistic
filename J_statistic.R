# Load libraries 
# Install libraries if not already installed
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
  make_option(c("-n", "--nrun"), type="integer", default=10, 
              help="Number of permutation runs/", metavar="NUMBER_OF_RUNS"),
  make_option(c("-h", "--het"), type="integer", default=10, 
              help="Minimum number of heterozygous SNVs.", metavar="MIN_NUMBER_OF_SNVS"),
  make_option(c("-x", "--seed"), type="integer", default=12345, 
              help="Random seed.", metavar="SEED"))

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE);
opt = parse_args(opt_parser)

# CHECK: Mandatory SNP argument
if(is.null(opt$snp)) { 
  stop("Missing: Absolute file path of SNP calls (CSV file).")
}

# CHECK: Mandatory CNV argument
if(is.null(opt$cnv)) { 
  stop("Missing: Absolute file path of CNV calls (CSV file).")
}

# CHECK: Mandatory Output Directory argument
if(is.null(opt$output)) { 
  stop("Missing: Absolute file path of output directory.")
}

cat("Stage 1: Setup user-input parameters.\n")

calls = read.csv(opt$snp,check.names=FALSE)
cnvsData = read.csv(opt$cnv,check.names=FALSE)
output_directory = opt$output
numRuns = opt$nrun
heterozygousCallCutoff = opt$het
seed = opt$seed

########################################################
#Stage 2: Pre-process input data for J-statistic script#
########################################################

cat("Stage 2: Pre-process input data for J-statistic script.\n")

#Load Affy position (build 38) reference probe dataset
probeData <- read.csv("./data/GGG.csv", check.names=FALSE) 
colnames(probeData)[1]<-"probeset_id"

#Create output folder
dir.create(output_directory, showWarnings = FALSE) 

#Initialize Stage 2 Progress Bar
stage_2_pb <- txtProgressBar(min = 0, max = length(2:ncol(calls)), style = 3)
stage_2_pb_counter = 0

#Prepare input data by merging CNV and SNV calls into a single dataframe
for (i in 2:ncol(calls)){
  
  sampleData<- calls[,c(1,i)]
  
  completeData<- merge(sampleData,probeData[,c(1,2,3)])# all =TRUE gives you the NA values
  
  sampleName <- tools::file_path_sans_ext(colnames(completeData)[2])
  
  nameData <- cnvsData[grep(sampleName,cnvsData$'Sample'),]
  
  if (nrow(nameData)==0){
    
    write.csv(NA, file.path(output_directory, paste0(sampleName,"-NULL", ".csv")), row.names=FALSE, na='' )
    next
  }
  nameData <- nameData[order(nameData$'Chromosome',nameData$'Start'),]
  nameData<-nameData[,-c(ncol(nameData),ncol(nameData)-1)] # gets rid of last 2 columns
  nameData<- nameData[grep("X|Y|M",nameData$'Chromosome', invert=TRUE),]# gets rid of X/Y/M

  colnames(completeData)[2]<-"Sample"
  colnames(completeData)[3]<-"Chromosome"
  colnames(completeData)[4]<-"Position"

  completeData<- completeData[grep("X|Y|M",completeData$'Chromosome', invert=TRUE),]
  completeData <- completeData[order(as.numeric(as.character(completeData$'Chromosome')),completeData$'Position'),]
  
  completeData$ID <- seq.int(nrow(completeData))
  nameData$ID <- seq.int(nrow(nameData))
  outputtable<-merge(completeData,nameData, by ="ID", all.x=T, check.names=F)
  outputtable$'ID'<-NULL
  names(outputtable) <-gsub(".[x-y]","",names(outputtable))
  
  colnames(outputtable)[5]<-"CNV.Chromosome"
  
  write.csv(outputtable, file.path(output_directory, paste0(sampleName,".csv")), row.names=FALSE, na='' )
  
  stage_2_pb_counter = stage_2_pb_counter + 1
  setTxtProgressBar(stage_2_pb, stage_2_pb_counter)
}

close(stage_2_pb)

############################################
#Stage 3: J-statistic calculation and plots#
############################################

cat("Stage 3: J-statistic calculation and plots. \n")

#Load reference SNP CNV association dataset
load(file = "./data/SNP cluster detection - SNP CNV association test - FUNCTION v4.RData")

#Initialize dataframe for the Excel document output with the generated statistics 
excel_data<-data.frame()

#Initialize Stage 3 Progress Bar
stage_3_pb <- txtProgressBar(min = 0, max = length(list.files(path = output_directory, pattern= paste0('^(SNP_).*(.csv)$'), full.names=TRUE)), style = 3)
stage_3_pb_counter = 0

#Reads files in folder in a loop( files within the homedirectory with the SNP_ as prefix and .csv as suffix)
for (file in list.files(path = output_directory, pattern= paste0('^(SNP_).*(.csv)$'), full.names=TRUE)){

  excel_table<-data.frame(matrix(ncol=5,nrow=19))
  colnames(excel_table)<- c("Chr","P.J.",'Hetero%','Num_CNV','Association')
  
  filename<-gsub('.csv','',file)
  filenamesplit <- strsplit(filename, "/")[[1]]
  fileNoSuffix <- filenamesplit[length(filenamesplit)]
  
  dir.create(paste0(output_directory, '/', fileNoSuffix), showWarnings = F)

  cancer.data <- read.csv(paste0(file), header=T)
  
  for (i in seq(1:19)){

    if (is.na(i)){
      print('i = NA, ending program')
      
      break
    }
    
    ############# SNP cluster detection ################
    
    # Take the subset of a chromosome
    cancer.data.19 <- subset(cancer.data, Chromosome == i)
    
    # Obtain the probe set on chromosome
    probe.set.19 <- cancer.data.19$Position
    
    # Get pseudo SNP sites
    SNP.example <- probe.set.19[which(cancer.data.19[, "Sample"] == "1")]
    
    percentHetero<-(length(SNP.example))/(length(probe.set.19))
    
    #If chromosome has no CNVs, then just make folder 
    if (!(i %in% (unique(cancer.data$CNV.Chromosome)))){
      #print (paste0('Chr',i,' has no CNVs'))
      dir.create(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_chr ', i,'-NoCNVs'), showWarnings = FALSE) #creates a new folders
      excel_table[i,]<-c(i, "NoCNVs",percentHetero, 0,"N")
      
      next
    }
    
    # Number of CNVs for each chr 
    numberCNVs<- sum((cancer.data$CNV.Chromosome[!is.na(cancer.data$CNV.Chromosome)])==i)
    
    
    # If less than desired cutoff of heterozygous calls it skips the chromosome
    if (length(SNP.example)<= heterozygousCallCutoff){
      print(paste0('less than 10 Heterozygous calls in chr',i))
      dir.create(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_chr ', i,'-LessThanCutoff'), showWarnings = FALSE) #creates a new folders
      excel_table[i,]<-c(i, "<CUtoff",percentHetero, numberCNVs,"<")
      next
    }

    dir.create(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_chr ', i), showWarnings = FALSE) #creates a new folders

    plot(log10(diff(SNP.example)) ~ SNP.example[-1], pch = 16)
    dev.copy(pdf, paste0(output_directory,'/',fileNoSuffix,'/',fileNoSuffix,'_chr ', i,'/',fileNoSuffix,' Chr', i,' Rainfall Plot.pdf'))
    dev.off()
    
    
    ###############################################################################
    
    # set grid points R_tilde(d) is evaluated at
    # set d
    grid.points <- 1:20 * 5000

    ######################################################################
    # Hypothesis testing for SNP cluster existence
    
    # Directory the file to be saved atD
    save.directory.SNP <- paste0(output_directory,'/',fileNoSuffix,'/',fileNoSuffix,'_chr ', i, '/')
   
    # File name, RData 
    file.name.SNP <- paste0(fileNoSuffix, "Chr ",i," Heterozygosity Clustering")

    # Set up the number of runs for the null distribution simulation
    n.run.SNP <- numRuns
    
    # Set seed for reproduction purposes
    set.seed(seed)
    
    start.time <- proc.time()[3]
    
    set.seed(seed)
    test.result.SNP <- SNP.cluster.test(detected.SNP = SNP.example, grid.d = grid.points, probe.set = probe.set.19, n.null = n.run.SNP, na.rm = T, alpha = 0.05)
    
    end.time <- proc.time()[3]
    total.time <- end.time - start.time
    total.time
    
    # Save the object test.result.SNP as the result of the run
    save(test.result.SNP, file = paste(save.directory.SNP, file.name.SNP, ".RData", sep = ""))

    # Load the saved object
    load(file = paste(save.directory.SNP, file.name.SNP, ".RData", sep = ""))

    # List out the components in test.result.SNP
    attributes(test.result.SNP)
    
    # Hypothesis test result; main result
    test.result.SNP$test.result
    
    # Calculation related to the test sample
    test.result.SNP$sample.info
   
     # Info on null distribution  
    test.result.SNP$Null.output

    #######################################################################################
    
    ############ SNP CNV association test ####################

    sample.CNVs <- cancer.data
    sample.CNVs.subset <- subset(sample.CNVs, CNV.Chromosome == i)
    sample.CNVs.start <- sample.CNVs.subset$Start
    sample.CNVs.end <-sample.CNVs.subset$End
    
    sample.CNVs.good <- matrix(c(sample.CNVs.start, sample.CNVs.end), ncol=2)
    
    sample.CNVs.2 <- Intervals(sample.CNVs.good, closed = c( TRUE, TRUE ), type = "Z")
    
    sample.SNPs.subset <- subset(cancer.data, (Chromosome == i) & (Sample == 1))
    sample.SNPs <- sample.SNPs.subset$Position
    
    sample.probe.subset <- subset(cancer.data, (Chromosome == i))
    sample.probe <- sample.probe.subset$Position
    
    sample.CNVs.2<-sample.CNVs.2[!is.na(sample.CNVs.2)]
    
    class(sample.CNVs.2) 
    length(sample.SNPs) # 4906
    length(sample.probe) # 40893
    
    #New Rainbow Plot#
    SNP.CNV.dist <- distance_to_nearest(sample.SNPs, sample.CNVs.2)
    
    
    
    plot(sample.CNVs.2, axes = FALSE, y = rep(0.1, dim(sample.CNVs.2)[1]),
         xlim = range(sample.SNPs), 
         ylim = c(0, max(ceiling(log10(SNP.CNV.dist)))), col = 2, cex = 1.25)
    
    points(log10(SNP.CNV.dist) ~ sample.SNPs)
    abline(v = apply(sample.CNVs.2, 1, mean), lty = 2, col =3, lwd = 1.5)
    
    # Parameters below specifying axes may need to be altered based on different chromosome length
    
    # x-axis
    # at = 2*0:10*10^7, indicate where to make ticks on x axis
    # label = as.character(2*0:10), indicate what to label on the ticks specified by "at"
    axis(1, at = 2*0:10*10^7, label = as.character(2*0:10), font = 1)  
    # y-axis
    axis(2, at = 0:max(ceiling(log10(SNP.CNV.dist))), labels = c(1, sapply(1:max(ceiling(log10(SNP.CNV.dist))), function(i) as.expression(bquote(10^ .(i))))))  # y-axis
    # label on x axis
    mtext(1, text = bquote("Chromosome Position (" * 10^7 ~ "bp)"), cex = 1.2, line = 3)
    # label on y axis
    mtext(2, text = "Distance from SNPs to closest CNV (bp)", cex = 1.2, line = 2.5)
    title(main = paste0(fileNoSuffix," Chr", i))
    #End of New Rainbow Plot#
    dev.copy(pdf, paste0(output_directory, "/", fileNoSuffix,'/',fileNoSuffix,'_chr ', i, '/',fileNoSuffix, ' Chr ', i, ' Rainbow.pdf'))
    dev.off()
    
    #Beginning of J statistic Analysis#
    range(distance_to_nearest(sample.SNPs, sample.CNVs.2))
    
    # set grid r 
    #r_use_1 <- sort(unique(c(seq(0,10^6,100), seq(10^6, 10^7, 5 * 10^2), seq(10^7, 2.5*10^7, 10^3))))
    r_use_1 <- 1 : 2000 * 5000
    
    length(r_use_1)
    range(r_use_1)

    # Directory the file to be saved at
    save.directory.SNP.CNV <- paste0(output_directory,'/',fileNoSuffix,'/',fileNoSuffix,'_chr ', i, '/')
    # File name, RData 
    file.name.SNP.CNV <- paste0(fileNoSuffix," chr ",i, " Association Test")

    # Set up the number of runs for the null distribution simulation
    n.run.SNP.CNV <- numRuns
    start.time <- proc.time()[3]
    
    # Set seed for reproduction purpose
    set.seed(seed)
    test.result.SNP.CNV <- SNP.CNV.test(sample.SNPs = sample.SNPs, sample.probe = sample.probe, sample.CNVs = sample.CNVs.2, r_p = r_use_1, n.rep = n.run.SNP.CNV, alpha = 0.05)
    
    end.time <- proc.time()[3]
    total.time <- end.time - start.time
    total.time
    
    # Save the object test.result.SNP as the result of the run
    save(test.result.SNP.CNV, file = paste(save.directory.SNP.CNV, file.name.SNP.CNV, ".RData", sep = ""))
    
    # Load the saved object
    load(, file = paste(save.directory.SNP.CNV, file.name.SNP.CNV, ".RData", sep = ""))

    # List out the components in test.result.SNP
    attributes(test.result.SNP.CNV)

    # Hypothesis test result; main result
    test.result.SNP.CNV$test.result
    
    # Calculation related to the test sample
    test.result.SNP.CNV$sample.info
    
    # Info on null distribution  
    test.result.SNP.CNV$Null.output
    
    ###################################################################################
    
    # Negative association, larger value, tend to be farther away
    # Positive association, smaller value, tend to be close together
    # For small chromosomes: xlim = range(distance_to_nearest(sample.SNPs, sample.CNVs.2))
    
    plot.upper <- max(c(test.result.SNP.CNV$sample.info[["J.fun.sample"]], test.result.SNP.CNV$Null.output[["gcb.upper"]]), na.rm =  T) + 0.01
    plot.lower <- min(c(test.result.SNP.CNV$sample.info[["J.fun.sample"]], test.result.SNP.CNV$Null.output[["gcb.lower"]]), na.rm = T) - 0.01
    plot(test.result.SNP.CNV$sample.info[["J.fun.sample"]] ~ r_use_1, type = 'l', ylim = c(plot.lower, plot.upper), lwd = 2, xlab = 'r', ylab = 'J function', main = "J function")
    lines(test.result.SNP.CNV$Null.output[["gcb.upper"]] ~ r_use_1, lty = 2, lwd = 2, col = 4)
    lines(test.result.SNP.CNV$Null.output[["gcb.lower"]] ~ r_use_1, lty = 2, lwd = 2, col = 4)
    legend("topleft", c("Observed J Function", "Global Confidence Bands"), lty = c(1,2), lwd = c(2,2), col = c(1, 4), cex=0.5)
    dev.copy(pdf, paste0(output_directory, '/',fileNoSuffix,'/',fileNoSuffix,'_chr ', i, '/', fileNoSuffix,' Chr ',i,' J Statistic.pdf'))
    dev.off()
    
    #Print to terminal if 
    graphLinePoints <- test.result.SNP.CNV$sample.info[["J.fun.sample"]]
    graphUpperPoints <- test.result.SNP.CNV$Null.output[["gcb.upper"]]
    graphLowerPoints <- test.result.SNP.CNV$Null.output[["gcb.lower"]]
    
    if ((any(graphLinePoints[-1:-100] >graphUpperPoints[-1:-100])) &
        (any(graphLinePoints[-1:-100] <graphLowerPoints[-1:-100]))){
      association<-'B'
      #print (paste(file, ': J function crossed both high and low global confidence bands.', sep=""))
      
    } else if (any(graphLinePoints[-1:-100] <graphLowerPoints[-1:-100])){
      association<-'+'
      #print (paste(file, ': J function crossed low global confidence band.', sep=""))
    } else if (any(graphLinePoints[-1:-100] >graphUpperPoints[-1:-100])){
      association<-'-'
      #print(paste(file, ': J function crossed high global confidence band.', sep=""))
    } else {
      association<-'/'  # \ is special character, need 2 to produce 1 in end 
      #print(paste(file, ': J function did not cross either high and low global confidence bands.', sep=""))
    }
    
    #if P.J. unadjusted is >0.05, make the association = '\'
    if ((as.numeric(test.result.SNP.CNV$test.result)[1])>0.05){
      association<-'/'
    }

    chars <- capture.output(cat(fileNoSuffix,'Chromosome',i,sep='\n'),print(test.result.SNP$test.result), print(test.result.SNP.CNV$test.result),
                            cat('%Hetero',percentHetero,sep='\n'), cat('numCNVs',numberCNVs, sep='\n'),
                            cat('Association',association,sep='\n'))    
    
    writeLines(chars, con = file(paste0(output_directory, '/', fileNoSuffix,'/',fileNoSuffix,'_chr ', i, '/', fileNoSuffix,' Chr ',i,' stats.txt')))
    
    excel_table[i,]<-c(i, (as.numeric(test.result.SNP.CNV$test.result)[1]),percentHetero, numberCNVs,association)
    closeAllConnections()
    
  }
  
  excel_table[,6]<-''
  colnames(excel_table)<- c("Chr","P.J.",'Hetero%','Num_CNV','Association','-')
  excel_table[(i+1),]<-c(fileNoSuffix,'','','','','')
  
  rownames(excel_table)<-NULL
  
  excel_data<- cbind.fill(excel_data,excel_table,fill=NA)
  
  stage_3_pb_counter = stage_3_pb_counter + 1
  setTxtProgressBar(stage_3_pb, stage_3_pb_counter)
}

close(stage_3_pb)

excel_data[1]<-NULL
colnames(excel_data)<- rep(c("Chr","P.J.",'Hetero%','Num_CNV','Association',' '),length(list.files(pattern= paste0('^(SNP_).*(.csv)$'))))

write.table(excel_data, file=paste(output_directory, 'summary_output_statistics.csv', sep="/"),row.names=FALSE,col.names=T,sep=",",na='')
closeAllConnections()

