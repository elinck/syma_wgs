### create bamlists
ls -1 *mega_R1_sorted.bam > megarhyncha.txt
ls -1 *toro_R1_sorted.bam > torotoro.txt
ls -1 *ochr_R1_sorted.bam > ochracea.txt
ls -1 *sorted.bam > all_bams.txt # manually deleted hyrad files w/ nano 
ls -1 *sorted.bam > combined_bams.txt # all 40 samples

### estimate SFS per species with angsd
angsd -b /home/ubuntu/sorted/megarhyncha.txt -anc /home/ubuntu/masked.ref.fa \
-dosaf 1 -gl 1 -P 24 -doMajorMinor 1 -doMAF 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-2 -out megarhyncha \
 
angsd -b /media/burke/bigMac/ethan/alignment/torotoro.txt -anc /media/burke/bigMac/ethan/masked.ref.fa \
-dosaf 1 -gl 1 -P 16 -doMajorMinor 1 -doMAF 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-2 -out torotoro \

### get 2-d SFS prior
realSFS megarhyncha.saf.idx torotoro.saf.idx -P 24 > meg_tor.ml
 
### estimate Fst by SNP
realSFS fst index megarhyncha.saf.idx torotoro.saf.idx -sfs meg_tor.ml \
-P 16 -fstout meg_tor
 
### global Fst
realSFS fst stats meg_tor.fst.idx -P 16 > global_fst.txt
 
### sliding-window Fst
realSFS fst stats2 meg_tor.fst.idx -P 16 -win 50000 \
-step 10000 > ../meg_tor_fstwindow.txt

angsd -bam all_bams.txt -doSaf 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-2 \
-anc /home/ubuntu/masked.ref.fa -GL 1 -P 24 -out syma_spp 

realSFS syma_spp.saf.idx -P 24 > syma.sfs

angsd -bam all_bams.txt -doSaf 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-2 \
-anc /media/burke/bigMac/ethan/masked.ref.fa -GL 1 -P 24 -doThetas 1 -pest syma.sfs \
-out syma

### windowed thetas
angsd -bam all_bams.txt -doSaf 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-2 \
-anc /home/ubuntu/masked.ref.fa -GL 1 -P 24 -doThetas 1 -pest syma.sfs \
-out new_syma

thetaStat do_stat new_syma.thetas.idx -win 50000 -step 50000 \
-outnames syma.theta_windows.gz

### call genotypes on full dataset, 75% complete matrix 
angsd -bam combined_bams.txt -anc /home/ubuntu/masked.ref.fa -minQ 20 \
-fold 1 -out syma_angsd_snps_75p  -GL 1  -doGeno 2  -doPost 1  -postCutoff 0.95 \
-doCounts 1 -doGlf 2 -geno_minDepth 5 -SNP_pval 0.05 -minMaf 0.025 \
-doMaf 2 -doMajorMinor 1 -doSaf 1 -doPlink 2 -minInd 30 -P 24

