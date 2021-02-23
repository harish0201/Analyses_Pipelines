#######Denovo RNASeq analyses

### 1 QC
fastp -h ../clean/sample.fastp.html -j ../clean/sample.fastp.json -w 10 -l 30 -a auto -f 10 -F 10 -p -z 9 -i sample_1.fq.gz -I sample_2.fq.gz -o ../clean/sample_1.clean.fastq.gz -O ../clean/sample_2.clean.fastq.gz

#### 2 Pooling and normalization

cat *_1.clean.fastq.gz > R1.fq.gz
cat *_2.clean.fastq.gz > R2.fq.gz

#### 3 normalization
/apps/bbmap/bbnorm.sh target=40 fixspikes=t k=25 threads=40 rdk=t in=R1.fq.gz in2=R2.fq.gz out=norm.R1.fq.gz out2=norm.R2.fq.gz

#### 4 assemblies (for advanced pick all the assembler, for standard pick trinity)
conda activate assembler
rnaspades.py -1 norm.R1.fq.gz -2 norm.R2.fq.gz -o spades -t 40
rnabloom -ntcard -t 40 -l norm.R1.fq.gz -r norm.R2.fq.gz -k 25-65:20
docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity --seqType fq --left `pwd`/norm.R1.fq.gz --right `pwd`/norm.R2.fq.gz --max_memory 512G --CPU 96 --output `pwd`/trinity_out_dir --no_normalize_reads --full_cleanup

### 5 combining assemblies and removing redundancy
cat spades/transcripts.fasta trinity.fasta rnabloom.transcripts.fa > all_pooled.fa
/data/analysis/shared_files/Harish/apps/evigene17dec14/scripts/prot/tr2aacds2.pl -MINCDS=90 -NCPU=40 -MAXMEM=100000 -tidyup -mrnaseq all_pooled.fa &
cat okayset/*.fasta > okay.fasta
cat okayset/*.pep > okay.pep

#### 6 index for quants for rapclust (https://github.com/COMBINE-lab/RapClust ; locate and add mcl to path)
salmon index -t okay.fasta -p 96 -i okay
mkdir quants
salmon quant -i okay -l A -1 sample_1.clean.fastq.gz -2 sample_1.clean.fastq.gz -o quants/sample_1 --seqBias --gcBias --dumpEq
#prepare config file as per the website and login goto (harish@10.10.10.14)
activate python2
Rapclust --config config.yaml 
awk '{print $2"\t"$1}' mag.flat.clust > gene.map

#### 7 quantify all the samples using the command above in 6
salmon quant -i okay -l A -1 sample_1.clean.fastq.gz -2 sample_1.clean.fastq.gz -o quants/sample_1 --seqBias --gcBias --dumpEq --gene gene.map
salmon quantmerge --column=numReads --gene $(find -type d -maxdepth 1) --output counts.matrix
#will have to round of the values to the closest integer

#### 8 run differential expression after preparing the replicate matrix and contrast matrix
#### filter for log2foldchange > 1 and FDR/p-adj<0.05

#### 9 GO term enrichment
##fetch data from bioMart for the species of interest and if not available download protein sequences and upload on eggNOG mapper.
##for eggNOG: 
awk -F'\t' '{print $1"\t"$7}' <(grep "GO:" eggNOG_out.tsv) > GO.txt
Rscript topGO.r GO.txt gene.list

#### 10 KEGG enrichment
#fetch data from bioMart for the species of interest: entrezID and uniprotID
#grep genelist against it
Rscript Kegg_unknown.R gene.list org_id #(org_id: https://www.genome.jp/kegg/catalog/org_list4.html)


#######Reference RNASeq analyses
##AFTER QC

#### 1. Prepare splice-site files (assumes gtf input from ENSEMBL)
extract_splice_sites.py annotation.gtf > splice.hisat2

#### 2. Index genome with hisat2 (-p => threads)
hisat2-build -p ${nproc} genome.fa genome

#### 3. hisat2 alignment ($1==directory where QC PASS files are present; $2==splice.hisat2 from above step. if annotation is not available, remove -x "$2"; $3==genome index; $4==output directory)
for i in "$1"/*_R1.clean.fq.gz; do echo "/apps/hisat2-2.1.0/hisat2 -p 10 --dta --new-summary --summary-file "$4"/$(basename "$i" _R1.clean.fq.gz).log --mm -x "$2" --known-splicesite-infile "$3" -1 "$i" -2 "$1"/$(basename "$i" _R1.clean.fq.gz)_R2.clean.fq.gz | sambamba view -S /dev/stdin -f bam -t 10 -o /dev/stdout | sambamba sort --memory-limit=3GB -t 10 -n --tmpdir="$4"/tmp -o "$4"/$(basename "$i" _R1.clean.fq.gz).sort.bam /dev/stdin" ;done | parallel -j6 &

#### 4. feature counting (might need to change -g gene_id as per requirement; -T threads)
/apps/subread-1.5.0-p3-source/bin/featureCounts -p -B -T ${nproc} -g gene_id -a annotation.gtf -o counts.profile *bam
grep -v "#" counts.profile | cut -f1,7- > counts.matrix
sed -i 's/.sorted.bam//g' counts.matrix

#### 5. DEGS and the following pipeline from point 8,9,10 from above.





