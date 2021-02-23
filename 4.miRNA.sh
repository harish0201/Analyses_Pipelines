####miRNA analyses

#### 1. Index genome (download preferably from ensembl. ncRNAs from RNACentral, remove miRNA related data)
sed 's/ .*//g' genome.fa #(remove white-space from header)
bowtie-build genome.fa genome
bowtie-build ncRNA.fa ncRNA

#### 2. QC data (change adapter sequences as needed)
#cutadapt
for i in *.fq; do echo "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -m 17 -M 30 "$i" > /data/analysis/dr.nahidAli.iicb_new/smallRNA/clean_data/1.clean/$(basename "$i" .fq).clean.fq 2>  /data/analysis/dr.nahidAli.iicb_new/smallRNA/clean_data/1.clean/$(basename "$i" .fq).log"; done | parallel -j48 &

#fastp
for i in *.gz; do echo "fastp -l 17 --length_limit 30 -p -w 10 --adapter_fasta ../clean/adapter.fa -R $(basename "$i" .fq.gz) -j ../clean/$(basename "$i" .fq.gz).fastp.json -h ../clean/$(basename "$i" .fq.gz).html -i "$i" -o ../clean/"$i""; done | parallel -j10 &

#### 3. align to ncRNA
mkdir ncrna_map
for i in ../clean/*.fq.gz; echo "zcat "$i" | bowtie -p 10 --un ncrna_map/$(basename "$i" .fq.gz).fq ncrna /dev/stdin -S ncrna_map/$(basename "$i" .fq.gz).sam 2> ncrna_map/$(basename "$i" .fq.gz).log"; done | parallel -j10
	
#### 4. miRNA analyses 
#predict denovo miRNA
mkdir unmap_ncrna
cat ../ncrna_map/*fq | fastx_collapser -i /dev/stdin -o collapsed.fa
sed 's/>/>seq_/g' collapsed.fa | sed 's/-/_x/g' > new_lapse.fa
bowtie -p 40 -f /data/ref_genomes/finger_millet/genome new_lapse.fa | convert_bowtie_output.pl /dev/stdin > mapped.arf 1> genome.log #mirdeep needs arf file
nohup /data/analysis/shared_files/Harish/apps/miniconda3/bin/miRDeep2.pl new_lapse.fa /data/ref_genomes/finger_millet/genome.fa mapped.arf none none none 2>mirdeep.log (# for completely unknown species)
nohup /data/analysis/shared_files/Harish/apps/miniconda3/bin/miRDeep2.pl new_lapse.fa /data/ref_genomes/finger_millet/genome.fa mapped.arf mature_ref_miRNAs.fa mature_other_miRNAs.fa hairpin_ref_miRNAs 2>mirdeep.log (# for closer species if present)

#mature mirna analyses
#align to genome using above unmapped_ncrna with bowtie
#get sam files
/apps/subread-1.6.2-source/bin/featureCounts -T 48 -t miRNA_mature -g ID -o counts_profile -a ../mature.gff $(ls *sam)


#### 5. Differential expression is same as for RNASeq. Follow that.

#### 6. Target prediction
#pick significant miRNAs and upload them to psRNAtarget
#pick longest CDS using this: assumes ensembl formatted transcript catalogue.
zcat "$1" | sed 's/ .*gene:/|/g' | sed 's/ .*//g' > tmp.fa
samtools faidx tmp.fa
sed 's/|/\t/g' tmp.fa.fai | sort -k2,2 -rk3 --sort g | awk '!arr[$2]++' | awk '{print $1"|"$2}' | seqtk subseq tmp.fa - > longest_"$1".fa
rm tmp.fa tmp.fa.fai

#### 7. Ontology and pathway enrichment
#pick genes from the psrnatarget output and use them for GO and KEGG
