#### Genome Assembly Pipeline

#1. Convert Pacbio bam to fasta:
bam2fasta -o movie movie.subreads.bam

#2. Genome Survey steps: (for both pacbio and illumina, for 10X based data, complete till scaff_reads in step 5)
ls *.fastq *.fasta > input.fofn #(check for file name)
ntcard -t 48 -k 17,19,21,29,31 -c 10000 -p kmer @input.fofn
sed -i 'g/^F/d' *hist
sed -i 's/\t/ /g' *hist
#upload to genomescope to check the metric: http://qb.cshl.edu/genomescope/
# you can use jellyfish, but it takes a lot of time and space

#3.1 Configure falcon/Unzip for assembly:
#get config from here: https://github.com/PacificBiosciences/pb-assembly and modify for SGE, for both falcon and falcon-unzip
#if the runs for falcon are successful, you should get the following folders: 0-raw*, 1-preads_ovl, 2-asm-falcon
#if unzip is successful, you will get the following folders: 3-unzip and 4-polish
#falcon unzip polishes the assembly to phase it

#3.2 Configure canu for assembly: X=genome size in Gb check canu documentation for configuration and usage with SGE
canu -p asm -d dir  genomeSize=Xg -pacbio-raw *.fast*
#check step 9 for polishing. 
#3.3 WTDBG2 assembly
wtdbg2 -t 96 -x sq -o prefix $(ls *.fasta.gz)
wtpoa-cns -t 96 -i prefix.ctg.lay.gz -fo prefix.raw.fa
#check step 9 for polishing
#please check for HiFi usage for all the tools

#4. Duplicate haplotigs removal:
for i in pacbio.*.fasta.gz; do echo "minimap2 -x map-pb -t 40 assembly.fasta "$i" | gzip -c - > "$i".read.paf.gz"; done | parallel -j4 
split_fa assembly.fasta > assembly.fasta.split
minimap2 -xasm5 -DP -t 48 assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
pbcstat *.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed assembly.fasta
mv purged.fa primary.fasta
mv haps.fa alts.fasta

#5. Scaffolding with 10X data
#goto directory with 10X data &&:
ls *gz | sort -uVk1 > data.fofn
#open the file in terminal and mod the file to append the q1 and q2. Ensure the reads are paired (example below)
``
q1=PT10X19B_PT10X19A-AK1843_S1_L004_R1_001.fastq.gz
q2=PT10X19B_PT10X19A-AK1843_S1_L004_R2_001.fastq.gz
q1=PT10X19B_PT10X19A-AK1844_S2_L004_R1_001.fastq.gz
q2=PT10X19B_PT10X19A-AK1844_S2_L004_R2_001.fastq.gz
q1=PT10X19B_PT10X19A-AK1845_S3_L004_R1_001.fastq.gz
q2=PT10X19B_PT10X19A-AK1845_S3_L004_R2_001.fastq.gz
q1=PT10X19B_PT10X19A-AK1861_S4_L004_R1_001.fastq.gz
q2=PT10X19B_PT10X19A-AK1861_S4_L004_R2_001.fastq.gz
``
scaff_reads -nodes 48 data.fofn genome-BC_1.fastq.gz genome-BC_2.fastq.gz #(debarcodes the data, using 48 threads)
scaff10X -nodes 48 -align bwa -score 20 -matrix 2000 -read-s1 10 -read-s2 10 -longread 1 -gap 100 -edge 50000 -link-s1 8 -link-s2 8 -block 50000 primary.fasta genome-BC_1.fastq.gz genome-BC_2.fastq.gz scaff10x.fasta

#6. Assemble organelles
#assumes closer species mitome are present in mitos.fasta or chloroplast.fasta
for i in pacbio.*.fasta.gz; do echo "minimap2 -ax map-pb -t 40 mitos.fasta "$i" | gzip -c - > "$i".sam"; done | parallel -j4
#extract reads from sam
samtools fasta -F 4 -@ 48 -s file.fasta file.sam
cat file*fasta > mito.reads.fasta
#depending on number of reads sample it to either 1% or 5%
seqtk sample -s 11 mito.reads.fasta 0.05 > 5perc_11seed.fasta
#get -g parameter from closer species 
flye --pacbio-raw 5perc_11seed.fasta -g 20kb -o mito -t 96
#similarly for chloroplast

#7. Scaffold with Optical Maps
#assemble maps check config.xml depending on chemistry and machine
python pipelineCL.py -T 50 -i 5 -b genome.bnx -l genome_assembly -t /path/to/refaligner -a /path/to/config.xml -V 0 -A -m
#generate hybrid scaffolds
perl hybridScaffold.pl -n primary.fasta -b OM.cmap -c hybridScaffold_DLE1_config_3.3.xml -r /path/to/refaligner -B 2 -N 2 -f -o hybrid

#8. Scaffold with HiC
#Alignment to the genome: #combine _HYBRID_SCAFFOLD.fasta and _HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta to create primary.OM.fasta. You'll need the enzyme cut site
samtools faidx primary.OM.fasta
bwa index primary.OM.fasta
nohup 'bwa mem -5SP -t 96 primary.OM.fasta R1.fq R2.fq | samblaster | samtools view -@ 96 -S -h -b -F 2316 -o aligned.bam /dev/stdin' &
#SALSA2: will need convert.sh for further manual curation, never worked properly for me. check phasegenomics github repository for proper converter. https://github.com/phasegenomics/juicebox_scripts from make agp from fasta to the bottom of the page
bamToBed -i alignment.bam && LC_ALL=C sort --parallel=8 -k 4 /dev/stdin > alignment.bed
python run_pipeline.py -a primary.OM.fasta -l primary.OM.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds -p yes -m yes

#ScaffHiC:
scaffhic -nodes 48 -score 200 -depth 50 -length primary.OM.fasta.fai -fq1 R1.fq -R2 R2.fq primary.OM.fasta genome-hic.fasta
scaffhic -nodes 48 -score 200 -depth 50 -length primary.OM.fasta.fai -map genome-hic.map -plot genome-hic-length.png -file 0 -fq1 R1.fq -fq2 R2.fq genome-hic.fasta genome-hic2.fasta #(no improvements in genome-hic2.fasta is usually seen)
#can been seen in PretextView: for further manual curation, never worked properly for me. check phasegenomics github repository for proper converter. https://github.com/phasegenomics/juicebox_scripts from make agp from fasta to the bottom of the page

#9.1 Arrow polishing:
pbmm2 align ref.fasta movie.subreads.bam ref.movie.bam --sort -j 4 -J 2 (change -j && -J as per the system config and number of parallel alignments)
gcpp -j 6 -r ref.fasta -o ref.polish.fastq,ref.polish.vcf ref.movie.bam
##(tips: if multiple bams are present, use bamtools to slice them contig/scaffold wise: bamtools split -reference ref.movie.bam; please check how the files have to be merged and then gcpp can be run for each contig parallelizing it)

#9.2 Pilon polishing:
#check resequencing pipeline for mapping Illumina reads
pilon --genome ref.fasta --fix all --vcfqe --output ref --outdir pilon_seq --changes --threads 4
