#truseq resequencing SNP calling pipeline
setwd("~/Dropbox/syma_wgs")
library(magrittr);library(ggplot2);library(foreach);library(doMC)
registerDoMC(cores=4)

#######################################################################################
# setup: you'll need R, java (v1.8), samtools, picard tools, gatk, and AdapterRemoval #
# Budget around 5x the initial download size for hard drive space

######################################################################################################
#run AdapterRemoval to trim adapters, merge overlapping paired reads, and drop low-quality base calls#
R1 <- list.files("demultiplexed_reads",full.names=T) %>% grep("R1_001",.,value=T)
R2 <- list.files("demultiplexed_reads",full.names=T) %>% grep("R2_001",.,value=T)
commands <- c()
for(i in 1:40){
  commands[i] <- paste0("AdapterRemoval --file1 ",R1[i],
                        " --file2 ",R2[i],
                        " --basename trimmed/", R1[i] %>% basename() %>% strsplit("_") %>% unlist() %>% .[1],
                        " --trimns --trimqualities --collapse --threads 30")
}

for(i in commands){
  system(i)
}

#foreach(i=commands) %dopar% system(i) #untested parallel version