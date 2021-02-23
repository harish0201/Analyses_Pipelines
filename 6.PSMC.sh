#### PSMC (assumes alignments with BWA and sorted bam (population demography))

#### 1. Run samtools mpileup, bcftools and vcfutils: needs a really old version - attached with script in a zip file
samtools mpileup -C50 -uf genome.fasta -r sorted.bam | bcftools view -c /dev/stdin | vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > genome.fq

#(loop as required)

#### 2. Remove sex chromosomes from the fq list. I used minimap2/mash etc to identify the scaffolds and delete them (check GIB writeup)

#### 3. Diploid fq to psmcfa
fq2psmcfa -q20 pilon.fq > pilon.fa
cat *.fa > diploid.psmcfa

#### 4. Run 100-bootstrapped PSMC 
parallel -j $NSLOTS "psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" diploid.psmcfa -o diplod.round-{}.psmc" :::: <(seq 100)

#### 5. Run psmc-plot to visualize the results (-u parameter is mutation rate, get from client or literature survey. -g is generation time, get this as well)
psmc_plot.pl -u 4e-09 -g 5 combined.plot combined.psmc  
