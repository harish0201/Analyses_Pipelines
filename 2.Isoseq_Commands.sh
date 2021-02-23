####Isoseq commands (followed from this website: https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md) (set polishing on if you user gui)

#### 1. CCS
conda activate pacbio 
ccs barcoded.bam barcoded.ccs.bam --min-rq 0.99 -j 300 --report-file barcoded.report.txt #(loop if needed)


#### 2. Lima
lima barcoded.ccs.bam primer.fasta lima.bam --isoseq --peek-guess -j 100 #(loop if needed)

#### 3. Refine
isoseq3 refine --require-polya lima.bam primer.fasta fl.bam -j 100 #(loop if needed)

#### 4. Cluster
ls fl*.bam > flnc.fofn
isoseq3 cluster flnc.fofn clustered.bam --verbose --use-qvs -j 400 

#### 5. Polish 
dataset create --type SubreadSet merged.subreadset.xml barcoded*.bam
isoseq3 polish clustered.bam merged.subreadset.xml barcoded*.bam -j 400

#### 6. Coding and Non-Conding analyses
#Read the following wiki: https://github.com/GenomeRIK/tama/wiki
#Functional annotation with interproscan, eggNOG and pannzer2
#Non-Coding annotation with PLEK, CNCI, CPAT



###For GUI
#### 1. Create XML that can be loaded into GUI
for i in *bam; do echo "dataset create --force --type SubreadSet $(basename "$i" .bam).subreadset.xml "$i""; done | parallel -j4

#### 2. Import Data using the Data Management Portal

#### 3. Select the data for analyses SMRT Analyses Portal and use Iso-Seq from the dropdown

#### 4. Advanced Parameters: Select Polish CCS on and change compute settings

#### 5. Download all the HQ fasta/q files as this is what is going to be used later


####For TAMA

#### 1. Align with minimap2 to genome (hq.fasta)
minimap2 -ax splice --secondary=no -C5 -uf -t 48 Reference.fasta hq.fasta > hq.sam
samtools view -@ 24 -Sbo hq.bam hq.sam
samtools sort -@ 24 -o hq.sorted.bam hq.bam
samtools view -h -o hq.sorted.sam hq.sorted.bam

#### 2. Tama collapse for capped and uncapped transcripts for degradation signature (https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-Degradation-Signature and report if signature greater than 0.25)
python tama_collapse.py -s hq.sorted.bam -f Reference.fasta -p capped -x capped 
python tama_collapse.py -s hq.sorted.bam -f Reference.fasta -p nocap -x no_cap
python tama_degradation_signature.py -c capped_trans_read.bed -nc nocap_trans_read.bed -o degrade.txt

#### 3. Tama merge models (https://github.com/GenomeRIK/tama/wiki/Tama-Merge create confif file as per here. you can use only capped transcriptome bed from above)
python tama_merge.py -f filelist.txt -p merged_annos

#### 4. Tama ORF-NMD prediction (bed files from above, change threads and number of jobs wherever applicable) (for blast check Uniref and OrthoDBv10 for the clade)
bedtools getfasta -name -split -s -fi Reference.fasta -bed merged_annos.bed -fo merged_annos.fasta
python tama_orf_seeker.py -f merged_annos.fasta -o merged_annos.pep.fasta
seqkit split2 -p20 merged_annos.pep.fasta && cd merged_annos.pep.fasta.split
for i in *fasta; do echo "blastp -evalue 1e-10 -num_threads 16 -db uniref50.fasta -query "$i" > "$i"blast.out"; done | parallel -j5
cat *blast.out > output.bls
python tama_orf_blastp_parser.py -b output.bls -o output.tama.parsed
python tama_cds_regions_bed_add.py -p merged_annos.fasta -a merged_annos.bed -f Reference.fasta -o final.cds.bed
bedtools getfasta -name -split -s -fi Reference.fasta -bed final.cds.bed -fo cds.fasta

####FOR SQUANTI3 (installed on 10.10.10.14 under SQUANTI3.env)
conda activate SQUANTI3.env
export PYTHONPATH=$PYTHONPATH:/apps_940/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/apps_940/cDNA_Cupcake/
export PATH=$PATH:/apps_940/SQUANTI3/utilities/
##gtf conversion help:
/apps_940/stringtie-1*/gffread genome.gff3 -T -o genome.gtf #(if this doesn't work for whatever reason: in monodon case it was augustus gff which is not well formatted)
#cat genome.gff3 | /apps_940/EVidencemodeller/EVM_utils/misc/convert_augustus_to_gff3 /dev/stdin > genome.mod.gff#
python /apps_940/SQUANTI3/sqanti3_qc.py -t 48 -o prefix cds.fasta genome.gtf Reference.fasta
