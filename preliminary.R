setwd("~/Dropbox/syma_wgs")

library(vcfR);library(adegenet);library(ggplot2);
library(ape);library(strataG);library(data.table);
library(pcadapt);library("qvalue");
library("OutFLANK");library("ggplot2")
library(vcfR);library(PopGenome)
library(plyr)

#50% 
syma50 <- read.vcfR("syma.50.recode.vcf")

# 75% complete matrix
syma <- read.vcfR("syma.filtered2.recode.vcf")
syma.genind <- vcfR2genind(syma)
seq.scaled <- scaleGen(syma.genind,NA.method="mean",scale=F)
snpID <- as.data.frame(syma.genind$loc.n.all)

# 85% complete matrix
syma.med <- read.vcfR("syma.med.filtered.recode.vcf")
syma.med.genind <- vcfR2genind(syma)
seq.med.scaled <- scaleGen(syma.genind,NA.method="mean",scale=F)
snpID <- as.data.frame(syma.genind$loc.n.all)

# 100% complete matrix
syma.full <- read.vcfR("syma.hard.filtered.recode.vcf")
syma.full.genind <- vcfR2genind(syma.full)
seq.full.scaled <- scaleGen(syma.full.genind,NA.method="mean",scale=F)

samples <- c("EL10_toro","EL11_toro","EL13_toro","EL18_mega",
             "EL19_mega","EL1_mega","EL20_mega","EL21_toro",
             "EL23_mega","EL24_mega","EL27_mega","EL29_ochr",
             "EL32_toro","EL39_toro","EL40_toro","EL4_mega",
             "EL5_ochr","EL6_mega","EL8_toro","EL9_toro")
pc$sampleID <- samples
rownames(syma.genind@tab) <- samples
rownames(syma.full.genind@tab) <- samples

clust.k1 <- find.clusters(syma.genind,n.pca=95,n.clust = 1,choose.n.clust = F)
clust.k2 <- find.clusters(syma.genind,n.pca=95,n.clust = 2,choose.n.clust = F)
clust.k3 <- find.clusters(syma.genind,n.pca=95,n.clust = 3,choose.n.clust = F)
clust.k4 <- find.clusters(syma.genind,n.pca=95,n.clust = 4,choose.n.clust = F)
clust.k5 <- find.clusters(syma.genind,n.pca=95,n.clust = 5,choose.n.clust = F)
clust.k6 <- find.clusters(syma.genind,n.pca=95,n.clust = 6,choose.n.clust = F)
clust <- cbind(sampleID=rownames(syma.genind@tab),clust.k1=unname(clust.k1$grp),clust.k2=unname(clust.k2$grp),clust.k3=unname(clust.k3$grp),
               clust.k4=unname(clust.k4$grp),clust.k5=unname(clust.k5$grp),clust.k6=unname(clust.k6$grp)) %>% data.frame()

clust.full.k1 <- find.clusters(syma.full.genind,n.pca=95,n.clust = 1,choose.n.clust = F)
clust.full.k2 <- find.clusters(syma.full.genind,n.pca=95,n.clust = 2,choose.n.clust = F)
clust.full.k3 <- find.clusters(syma.full.genind,n.pca=95,n.clust = 3,choose.n.clust = F)
clust.full.k4 <- find.clusters(syma.full.genind,n.pca=95,n.clust = 4,choose.n.clust = F)
clust.full.k5 <- find.clusters(syma.full.genind,n.pca=95,n.clust = 5,choose.n.clust = F)
clust.full.k6 <- find.clusters(syma.full.genind,n.pca=95,n.clust = 6,choose.n.clust = F)
clust.full <- cbind(sampleID=rownames(syma.genind@tab),clust.k1=unname(clust.k1$grp),clust.k2=unname(clust.k2$grp),clust.k3=unname(clust.k3$grp),
               clust.k4=unname(clust.k4$grp),clust.k5=unname(clust.k5$grp),clust.k6=unname(clust.k6$grp)) %>% data.frame()

pca <- prcomp(seq.scaled,center=F,scale=F)
screeplot(pca)
pc <- data.frame(pca$x[,1:3])
ggplot(data=pc,aes(x=PC1,y=PC2,col=clust$clust.k3))+geom_text(aes(label=sampleID))

pca.full <- prcomp(seq.full.scaled,center=F,scale=F)
screeplot(pca)
pc.full <- data.frame(pca.full$x[,1:3])
pc.full$sampleID <- samples
ggplot(data=pc.full,aes(x=PC1,y=PC2,col=clust.full$clust.k3))+geom_text(aes(label=sampleID))

# prelim trees
d <- dist(syma.genind)
plot(nj(d), type="fan")
d2 <- dist(syma.full.genind)
plot(nj(d2))

plot(upgma(d))
plot(upgma(d2))

# convert vcf to pcadapt
pcadapt::vcf2pcadapt("syma.raw.scan.recode.vcf", "pcadapt.snps")
snps <- read.table("pcadapt.snps", head = TRUE)
snps <- read.pcadapt(snps)
x <- pcadapt(snps)
pop <- c(1,1,1,2,2,2,2,1,2,2,2,1,1,1,2,2,1,2,1,1)
plot(x, option = "scores", pop = pop)
plot(x , option = "manhattan", pop = pop, )

#save pval obj
pval <- x$pvalues

#better manhattan plot
ggdf <- as.data.frame(cbind(as.vector(-log10(pval)), c(1:length(pval))))
colnames(ggdf) <- c("pval", "position")
png(width=6,height=1.5,units="in",res=600,file="manhattan.png")
ggplot(ggdf,aes(x=position,y=pval))+theme_bw()+
  theme(legend.position="none",
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.3),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7)) +
  #scale_color_manual(values=rep(c("grey60","grey80"),length(levels(factor(fst$chr)))/2+1))+
  scale_y_continuous(limits=c(0,75),minor_breaks = NULL)+
  geom_point(data = ggdf, size=0.5) +
  geom_point(data=subset(ggdf,pval>=quantile(ggdf$pval,0.975, na.rm = TRUE)),size=0.5,col="red") +
  xlab("position") +
  ylab("-log10 pvalue")
dev.off()

### demographic inference
devtools::install_github('statgenlmu/jaatha');install.packages("coala")
library(jaatha);library(coala)

#write snp file to save position
pcadapt::vcf2pcadapt("syma.filtered2.recode.vcf.gz", "syma75.snps")
pos <- read.table("positions.txt")
pos <- pos$V1
tnex <- read.table("nexus/syma_snp.txt", na.strings = ".")
tnex <- tnex[,1:18581] #sketchy -- why does it drop a column???
df <- rbind(pos, tnex)
df <- df[ , colSums(is.na(df)) == 0]
df <- as.matrix(df)
newpos <- as.vector(df[1:1,])
sub <- df[,sample(ncol(df), 100)]
subpos <- as.vector(sub[1:1,])
df <- df[2:41,]
sub <- sub[2:41,]

# split df by windows of 10
n <- 10
col <- ncol(df)
list.df <- list()
list.df <- split(df, rep(1:ceiling(col/n), each=n, length.out = col))
list.pos <- split(newpos, rep(1:ceiling(col/n), each=n, length.out = col))

# write 300 segsites loci
segsites <- list()
for(i in 1:300){
  a <- matrix(list.df[[i]], nrow=40, ncol=10)
  b <- as.vector(list.pos[[i]])
  segsites[[i]] <- create_segsites(a, b, check = TRUE)
}

#models 
im <- coal_model(c(10, 10), loci_number = 300, loci_length = 1, ploidy = 2) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_migration(par_range("m", 0, 3), symmetric = TRUE) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  sumstat_jsfs()

hybrid <- coal_model(c(10, 10), loci_number = 300, loci_length = 1, ploidy = 2) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_migration(par_range("m", 0, 3), symmetric = TRUE, time = 2) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  sumstat_jsfs()

i <- coal_model(c(10, 10), loci_number = 300, loci_length = 1, ploidy = 2) +
  feat_mutation(par_range("theta", 1, 10)) +
  feat_pop_merge(par_range("t_split", 0.1, 2), 2, 1) + 
  sumstat_jsfs()

## sumstats 
sumstats_im <- calc_sumstats_from_data(im, segsites)
sumstats_hybrid <- calc_sumstats_from_data(hybrid, segsites)
sumstats_i <- calc_sumstats_from_data(i, segsites)

#hybrid model test
hybrid_sim <- create_jaatha_model(hybrid)
hybrid_empirical <- create_jaatha_data(sumstats_hybrid, hybrid_sim)
estimates_hybrid <- jaatha(hybrid_sim, hybrid_empirical, 
                    sim = 100, repetitions = 2, verbose = FALSE)
#im model test
im_sim <- create_jaatha_model(im)
im_empirical <- create_jaatha_data(sumstats_im, jaatha_im)
estimates_im <- jaatha(im_sim, im_empirical, 
                    sim = 100, repetitions = 2, verbose = FALSE)

#isolation model test
i_sim <- create_jaatha_model(i)
i_empirical <- create_jaatha_data(sumstats_i, i_sim)
estimates_i <- jaatha(i_sim, i_empirical, 
                       sim = 100, repetitions = 2, verbose = FALSE)

estimates_hybrid$loglikelihood #[1] -3000.998
estimates_im$loglikelihood #[1] -2136.74
estimates_i$loglikelihood #[1] -3008.813

