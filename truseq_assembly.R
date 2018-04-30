############################################################
# truseq resequencing SNP calling pipeline                 #
# modified by E. Linck, original wrappers by C.J. Battey   #
############################################################

setwd("/data/jdumbacher/Syma/syma_wgs/")
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("foreach", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("doMC", repo = "http://ftp.osuosl.org/pub/cran/")
library(magrittr);library(foreach);library(doMC)
registerDoMC(cores=4)

#######################################################################################
# setup: you'll need R, java (v1.8), samtools, AdapterRemoval; bbmap; NextGenMap;     #
# and angsd                                                                           #
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
# remove duplicate reads (bbmap::dedupe.sh; memory intensive   #
################################################################

R1 <- list.files("/data/jdumbacher/Syma/syma_wgs/trimmed",full.names=T) %>% grep("pair1",.,value=T)
R2 <- list.files("/data/jdumbacher/Syma/syma_wgs/trimmed",full.names=T) %>% grep("pair2",.,value=T)
commands <- c()
for(i in 1:20){
  sampleID[i] <- basename(R1[i]) %>% strsplit("_S...pair1.truncated") %>% unlist() %>% .[1]
  commands[i] <- paste0("dedupe.sh", 
                        " in1=", R1[i], 
                        " in2=", R2[i],
                        " out=", sampleID[i], ".dedup.fq",
                        " ac=f",
                        " t=30;",
                        " reformat.sh in=", sampleID[i], ".dedup.fq", 
                        " out1=", sampleID[i], "_R1.dedup.fq",
                        " out2=", sampleID[i], "_R2.dedup.fq")
}
for(i in 1:20){
  system(commands[i])
}

system("mkdir dedup; mv *dedup* dedup")

################################################################
# align trimmed reads to reference sample                      #
################################################################
R1 <- list.files("/data/jdumbacher/Syma/syma_wgs/dedup",full.names=T) %>% grep("R1.dedup",.,value=T)
R2 <- list.files("/data/jdumbacher/Syma/syma_wgs/dedup",full.names=T) %>% grep("R2.dedup",.,value=T)
commands <- c()
sampleID <- c()
for(i in 1:20){
  sampleID[i] <- basename(R1[i]) %>% strsplit("_R..dedup") %>% unlist() %>% .[1]
  commands[i] <- paste0("ngm -r /data/jdumbacher/Syma/syma_wgs/halcyon_senegalensis/B10K-DU-024-03.genomic.fa", #change this to your reference
                        " -1 ", R1[i], 
                        " -2 ", R2[i],
                        " -o ", sampleID[i],".sam",
                        " -t 48")
}
for(i in 1:20){
  system(commands[i])
}

################################################################
# sort files                                                   #
################################################################
sams <- list.files("/data/jdumbacher/Syma/syma_wgs/alignments",full.names=T)
commands <- c()
sampleID <- c()
for(i in 1:length(sams)){
  sampleID[i] <- basename(sams[i]) %>% strsplit(".sam") %>% unlist() %>% .[1] 
  commands[i] <-paste0("picard SortSam",
                      " I=", sams[i], 
                      " O=", sampleID[i],
                      ".bam",
                      " SORT_ORDER=coordinate")
}

for(i in 1:length(sams)){
  system(commands[i])
}

################################################################
# index files                                                 #
################################################################
bams <- list.files("/data/jdumbacher/Syma/syma_wgs/alignments",full.names=T) %>% grep("bam",.,value=T)
commands <- c()
sampleID <- c()
for(i in 1:length(bams)){
  sampleID[i] <- basename(bams[i]) %>% strsplit(".bam") %>% unlist() %>% .[1] 
  commands[i] <-paste0("samtools index -b ",
                       bams[i], 
                       " ", sampleID[i],
                       ".bam.bai")
}

for(i in 1:length(bams)){
  system(commands[i])
}

