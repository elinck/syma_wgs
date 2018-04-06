##estimate SFS per species with angsd

angsd -b /data/jdumbacher/Syma/syma_wgs/megarhyncha.txt -anc /data/jdumbacher/Syma/syma_wgs/ELA_S/ELA_S.fasta \
-dosaf 1 -gl 1 -P 30 -doMajorMinor 1 -doMAF 1 -SNP_pval 1e-2 -out megarhyncha \
 
angsd -b /data/jdumbacher/Syma/syma_wgs/torotoro.txt -anc /data/jdumbacher/Syma/syma_wgs/ELA_S/ELA_S.fasta \
-dosaf 1 -gl 1 -P 30 -doMajorMinor 1 -doMAF 1 -SNP_pval 1e-2 -out torotoro \
 
##get 2-d SFS prior (needs >150gb RAM)
~/angsd/misc/realSFS megarhyncha.saf.idx torotoro.saf.idx -P 32 > meg_tor.ml
 
##estimate Fst by SNP
~/angsd/misc/realSFS fst index megarhyncha.saf.idx torotoro.saf.idx -sfs meg_tor.ml \
-P 32 -fstout meg_tor
 
##global Fst
~/angsd/misc/realSFS fst stats meg_tor.fst.idx -P 32
 
##sliding-window Fst
~/angsd/misc/realSFS fst stats2 meg_tor.fst.idx -P 32 -win 50000 \
-step 10000 > ../meg_tor_fstwindow.txt