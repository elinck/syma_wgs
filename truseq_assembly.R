############################################################
# truseq resequencing SNP calling pipeline                 #
# modified by E. Linck, original wrappers by C.J. Battey   #
############################################################

setwd("/home/elinck/syma_wgs")
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("foreach", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("doMC", repo = "http://ftp.osuosl.org/pub/cran/")
library(magrittr);library(foreach);library(doMC)
registerDoMC(cores=4)

#######################################################################################
# setup: you'll need R, java (v1.8), samtools, AdapterRemoval; bbmap; and angsd       #
# Budget around 5x the initial download size for hard drive space                     #
#######################################################################################

#run AdapterRemoval to trim adapters, merge overlapping paired reads, and drop low-quality base calls#
R1 <- list.files("/data/jdumbacher/Syma/syma_wgs/raw",full.names=T) %>% grep("R1_001",.,value=T)
R2 <- list.files("/data/jdumbacher/Syma/syma_wgs/raw",full.names=T) %>% grep("R2_001",.,value=T)
commands <- c()
for(i in 1:20){
  commands[i] <- paste0("AdapterRemoval --file1 ", R1[i],
                        " --file2 ",R2[i],
                        " --basename trimmed/", R1[i] %>% basename() %>% strsplit("_L006_R._001") %>% unlist() %>% .[1],
                        " --trimns --trimqualities --collapse --threads 32")
}

for(i in commands){
  system(i)
}
#foreach(i=commands) %dopar% system(i) #untested parallel version

################################################################
# assemble reference sample                                    #
# ELA_S largest file (shrug)                                   #
################################################################

# spades command here

################################################################
# remove duplicate reads (bbmap::dedupe.sh)                    #
################################################################

R1 <- list.files("/data/jdumbacher/Syma/syma_wgs/trimmed",full.names=T) %>% grep("pair1",.,value=T)
R2 <- list.files("/data/jdumbacher/Syma/syma_wgs/trimmed",full.names=T) %>% grep("pair2",.,value=T)
sampleID <- basename(R1[i]) %>% strsplit("_L006_R._001") %>% unlist() %>% .[1]
commands <- c()
for(i in 1:20){
  commands[i] <- paste0("dedupe.sh", 
                        "in1=", R1[i], 
                        "in2=", R2[i],
                        "out1=", sampleID[i], "R1.dedup.fq",
                        "out2=", sampleID[i], "R2.dedup.fq",
                        "ac=f",
                        "t=30")
}
for(i in 1:20){
  system(commands[i])
}

################################################################
# align trimmed reads to reference sample                      #
################################################################

R1 <- list.files("/data/jdumbacher/Syma/syma_wgs/trimmed",full.names=T) %>% grep("R1.dedup",.,value=T)
R2 <- list.files("/data/jdumbacher/Syma/syma_wgs/trimmed",full.names=T) %>% grep("R2.dedup",.,value=T)
commands <- c()
for(i in 1:20){
  commands[i] <- paste0("bbmap.sh",
                        "in1=", R1[i], 
                        "in2=", R2[i],
                        "minratio=0.1",
                        "t=48",
                        "bamscript=bs.sh; sh bs.sh")
}
for(i in 1:20){
  system(commands[i])
}

