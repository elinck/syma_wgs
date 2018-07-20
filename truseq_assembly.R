############################################################
# truseq resequencing SNP calling pipeline                 #
# modified by E. Linck, wrappers inspired by C.J. Battey   #
############################################################

setwd("/media/burke/bigMac/ethan/")
install.packages("magrittr", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("foreach", repo = "http://ftp.osuosl.org/pub/cran/"); install.packages("rebus", repo = "http://ftp.osuosl.org/pub/cran/")
library(magrittr);library(foreach)


install.packages("inline", repo = "http://ftp.osuosl.org/pub/cran/")
install.packages("gam", repo = "http://ftp.osuosl.org/pub/cran/")
install.packages("Rcpp", repo = "http://ftp.osuosl.org/pub/cran/")
install.packages("ggplot2", repo = "http://ftp.osuosl.org/pub/cran/")
install.packages("RcppGSL", "http://ftp.osuosl.org/pub/cran/")


### setup: you'll need R, java (v1.8), samtools, bbmap, bowtie, picard

### mask repeats from reference using bbmask entropy approach
system("bbmask.sh in=B10K-DU-024-03.genomic.fa out=masked.ref.fa entropy=0.7")

### run bbduk to trim adapters
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

### run bbduk to quality trim
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

### align trimmed reads to reference sample
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
                        " vslow minratio=0.1 k=8 maxindel=200 bamscript=bs.sh >> alignment_log.txt 2>&1; sh bs.sh ")
}
for(i in commands){
  system(i)
}

### fuck I guess I have to use GATK after all
setwd("/media/burke/bigMac/ethan/alignment/sorted/") 
# can't figure out a reasonable regex for this because I suck at coding
bam_wgs <- c("EL1_mega_R1_sorted.bam", "EL10_toro_R1_sorted.bam", "EL11_toro_R1_sorted.bam", "EL13_toro_R1_sorted.bam", "EL18_mega_R1_sorted.bam",
             "EL19_mega_R1_sorted.bam", "EL20_mega_R1_sorted.bam", "EL21_toro_R1_sorted.bam", "EL23_mega_R1_sorted.bam", "EL24_mega_R1_sorted.bam",
             "EL27_mega_R1_sorted.bam", "EL29_ochr_R1_sorted.bam", "EL32_toro_R1_sorted.bam", "EL39_toro_R1_sorted.bam", "EL4_mega_R1_sorted.bam",
             "EL40_toro_R1_sorted.bam", "EL6_mega_R1_sorted.bam", "EL8_toro_R1_sorted.bam", "EL9_toro_R1_sorted.bam", "EL5_ochr_R1_sorted.bam")
bam_hyrad <- c("EL45_toro_R1_sorted.bam", "EL46_toro_R1_sorted.bam", "EL47_toro_R1_sorted.bam", "EL48_ochr_R1_sorted.bam", "EL49_toro_R1_sorted.bam",
               "EL41_toro_R1_sorted.bam",  "EL50_toro_R1_sorted.bam", "EL51_toro_R1_sorted.bam", "EL52_toro_R1_sorted.bam", "EL53_toro_R1_sorted.bam",
               "EL54_toro_R1_sorted.bam", "EL55_toro_R1_sorted.bam", "EL56_toro_R1_sorted.bam", "EL57_toro_R1_sorted.bam", "EL58_toro_R1_sorted.bam",
               "EL59_toro_R1_sorted.bam", "EL42_toro_R1_sorted.bam",  "EL60_toro_R1_sorted.bam", "EL43_toro_R1_sorted.bam", "EL44_toro_R1_sorted.bam")
commands <- c()
for(i in 1:length(bam_wgs)){
  sampleID <- basename(bam_wgs[i]) %>% strsplit("_R1_sorted.bam") %>% unlist() %>% .[1]
  commands[i] <- paste0("java -Xmx10g -jar /media/burke/bigMac/Dropbox/tools/picard.jar AddOrReplaceReadGroups",
                        " I=",bam_wgs[i],
                        " O=tmp.bam",
                        " R=/media/burke/bigMac/ethan/masked.ref.fa",
                        " RGID=1",
                        " RGLB=syma_wgs",
                        " RGPL=illumina",
                        " RGPU=pool1",
                        " RGSM=",sampleID,";",
                        
                        "mv tmp.bam ",bam_wgs[i],";",
                        
                        "samtools index ",bam_wgs[i]
  )
}
for(i in 1:length(bam_wgs)){
  system(commands[i])
}
for(i in 1:length(bam_hyrad)){
  sampleID <- basename(bam_hyrad[i]) %>% strsplit("_R1_sorted.bam") %>% unlist() %>% .[1]
  commands[i] <- paste0("java -Xmx10g -jar /media/burke/bigMac/Dropbox/tools/picard.jar AddOrReplaceReadGroups",
                        " I=",bam_hyrad[i],
                        " O=tmp.bam",
                        " R=/media/burke/bigMac/ethan/masked.ref.fa",
                        " RGID=2",
                        " RGLB=syma_hyrad",
                        " RGPL=illumina",
                        " RGPU=pool2",
                        " RGSM=",sampleID,";",
                        
                        "mv tmp.bam ",bam_hyrad[i],";",
                        
                        "samtools index ",bam_hyrad[i]
  )
}

for(i in 1:length(bam_hyrad)){
  system(commands[i])
}

### recalibrate q scores for aDNA w/ mapDamage
bam <- list.files("/media/burke/bigMac/ethan/alignment/sorted/",full.names=F) %>% grep(".bam$",.,value=T)
for(i in 1:length(bam)){
  commands[i] <- paste0("mapDamage -i ",bam[i],
                      " -r /media/burke/bigMac/ethan/masked.ref.fa --rescale")
  }
for(i in 1:length(bam)){
  system(commands[i])
}

### static calls w/ gatk
system("java -jar /media/burke/bigMac/Dropbox/tools/picard.jar CreateSequenceDictionary R=/media/burke/bigMac/ethan/masked.ref.fa O=masked.ref.dict;
java -Xmx40g -jar /media/burke/bigMac/Dropbox/tools/GenomeAnalysisTK.jar -T UnifiedGenotyper\
 -R /media/burke/bigMac/ethan/masked.ref.fa\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL10_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL11_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL13_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL18_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL19_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL1_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL20_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL21_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL23_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL24_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL27_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL29_ochr_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL32_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL39_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL40_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL41_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL42_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL43_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL44_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL45_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL46_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL47_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL48_ochr_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL49_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL4_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL50_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL51_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL52_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL53_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL54_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL55_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL56_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL57_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL58_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL59_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL5_ochr_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL60_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL6_mega_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL8_toro_R1_sorted.bam\
 -I /media/burke/bigMac/ethan/alignment/sorted/EL9_toro_R1_sorted.bam\
 -o syma_gatk.vcf")
system("mkdir syma_alignment; mv *am* syma_alignment")

### vcftools commands
system("vcftools --vcf syma_gatk.vcf --max-missing 0.85 --minDP 3 --minQ 30 --recode --recode-INFO-all --out syma.gatk.85p.d3.vcf")

system("vcftools --vcf syma.gatk.75p.d3.vcf.recode.vcf --missing-indv")

system("vcftools --vcf syma_gatk.vcf --max-missing 0.85 --minDP 5 --minQ 30 --maf 0.05 --recode --recode-INFO-all --out syma.gatk.85p.d5.maf05")

system("vcftools --vcf syma_gatk.vcf --max-missing 1.0 --minDP 5 --minQ 30 --maf 0.05 --recode --recode-INFO-all --out syma.gatk.100p.d5.maf05")
