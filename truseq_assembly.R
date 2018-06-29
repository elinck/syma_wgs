############################################################
# truseq resequencing SNP calling pipeline                 #
# modified by E. Linck, wrappers inspired by C.J. Battey   #
############################################################

setwd("/media/burke/bigMac/ethan/")
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("foreach", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("doMC", repo = "http://ftp.osuosl.org/pub/cran/")
library(magrittr);library(foreach)

# setup: you'll need R, java (v1.8), samtools, bbmap, bowtie, picard

# mask repeats from reference using bbmask entropy approach
system("bbmask.sh in=B10K-DU-024-03.genomic.fa out=masked.ref.fa entropy=0.7")

#run bbduk to trim adapters
R1 <- list.files("/media/burke/bigMac/ethan/syma_raw",full.names=T) %>% grep("R1",.,value=T)
R2 <- list.files("/media/burke/bigMac/ethan/syma_raw",full.names=T) %>% grep("R2",.,value=T)
commands <- c()
for(i in 1:length(R1)){
  commands[i] <- paste0("bbduk.sh in1=", R1[i],
                        " in2=",R2[i],
                        " out1=", R1[i] %>% basename() %>% strsplit(".fastq") %>% unlist() %>% .[1],".trimmed.fq",
                        " out2=", R2[i] %>% basename() %>% strsplit(".fastq") %>% unlist() %>% .[1],".trimmed.fq",
                        " ref=adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo")
}
for(i in commands){
  system(i)
}
system("mkdir syma_trimmed; mv *trimmed* syma_trimmed")

#run bbduk to quality trim
R1 <- list.files("/media/burke/bigMac/ethan/syma_trimmed",full.names=T) %>% grep("R1",.,value=T)
R2 <- list.files("/media/burke/bigMac/ethan/syma_trimmed",full.names=T) %>% grep("R2",.,value=T)
commands <- c()
for(i in 1:length(R1)){
  commands[i] <- paste0("bbduk.sh in1=", R1[i],
                        " in2=",R2[i],
                        " out1=", R1[i] %>% basename() %>% strsplit(".trimmed.fq") %>% unlist() %>% .[1],".clean.fq",
                        " out2=", R2[i] %>% basename() %>% strsplit(".trimmed.fq") %>% unlist() %>% .[1],".clean.fq",
                        " qtrim=r trimq=10")
}
for(i in commands){
  system(i)
}
system("mkdir syma_cleaned; mv *clean* syma_cleaned")

# align trimmed reads to reference sample
R1 <- list.files("/media/burke/bigMac/ethan/syma_cleaned",full.names=T) %>% grep("R1",.,value=T)
R2 <- list.files("/media/burke/bigMac/ethan/syma_cleaned",full.names=T) %>% grep("R2",.,value=T)
commands <- c()
sampleID <- c()
for(i in 1:length(R1)){
  commands[i] <- paste0("bbmap.sh", 
                        " in1=", R1[i], 
                        " in2=", R2[i],
                        " out=", R1[i] %>% basename() %>% strsplit(".clean.fq") %>% unlist() %>% .[1], ".sam",
                        " ref=masked.ref.fa", 
                        " vslow minratio=0.1 k=8 maxindel=200 bamscript=bs.sh; sh bs.sh")
}
for(i in commands){
  system(i)
}

system("mkdir syma_alignment; mv *am* syma_alignment")
