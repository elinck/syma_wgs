### create bamlists

ls -1 *mega_R1.sam > megarhyncha.txt
ls -1 *toro_R1.sam > torotoro.txt
ls -1 *ochr_R1.sam > ochracea.txt
ls -1 ochr_R1.bam > ochracea.txt


### estimate SFS per species with angsd

angsd -b /media/burke/bigMac/ethan/alignment/megarhyncha.txt -anc /media/burke/bigMac/ethan/masked.ref.fa \
-dosaf 1 -gl 1 -P 16 -doMajorMinor 1 -doMAF 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-2 -out megarhyncha \
 
angsd -b /media/burke/bigMac/ethan/alignment/torotoro.txt -anc /media/burke/bigMac/ethan/masked.ref.fa \
-dosaf 1 -gl 1 -P 16 -doMajorMinor 1 -doMAF 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-2 -out torotoro \

### get 2-d SFS prior
realSFS megarhyncha.saf.idx torotoro.saf.idx -P 16 > meg_tor.ml
 
### estimate Fst by SNP
realSFS fst index megarhyncha.saf.idx torotoro.saf.idx -sfs meg_tor.ml \
-P 16 -fstout meg_tor
 
### global Fst
~/angsd/misc/realSFS fst stats meg_tor.fst.idx -P 16
 
### sliding-window Fst
realSFS fst stats2 meg_tor.fst.idx -P 16 -win 50000 \
-step 10000 > ../meg_tor_fstwindow.txt

### sliding-window Tajima's D

angsd -bam all_bams.txt -doSaf 1 -anc /media/burke/bigMac/ethan/masked.ref.fa -GL 1 -P 24 -out global 

ealSFS out.saf.idx -P 24 > global.sfs

angsd -bam all_bams.txt -out global -doThetas 1 -doSaf 1 -pest out.sfs -anc /media/burke/bigMac/ethan/masked.ref.fa -GL 1

thetaStat do_stat global.thetas.idx

thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz

