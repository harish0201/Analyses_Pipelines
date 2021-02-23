#### Annotation Pipeline (Assumes isoseq data and rnaseq data)

#1. Align RNASeq data to genome (final curated genome=> genome.fasta)
cat *R1* > R1.fastq.gz
cat *R2* > R2.fastq.gz
##align with hisat2 using earlier commands. you will not have splice sites
stringtie aligned.bam -p 96 -o stringtie.gtf
gffread stringtie.gtf -g genome.fasta -w cds.fasta

#2. Denovo assembly of RNASeq reads
bbnorm.sh target=40 fixspikes=t k=25 threads=40 rdk=t in=R1.fastq.gz in2=R2.fastq.gz out=norm.R1.fq.gz out2=norm.R2.fq.gz
Trinity --seqType fq --left norm.R1.fq.gz --right norm.R2.fq.gz --max_memory 160G --CPU 96 --output trinity_out --no_normalize_reads

#3. PASA analyses to create gene-models
#Combine all types of cDNA evidences #Follow isoseq protocol to get gene-models from TAMA (I'll call it tama.cds.fasta)
cat tama.cds.fasta cds.fasta Trinity.fasta > transcripts.fasta
#create template file as per pasa documentation, it will either be called template.txt or alignAssembly.conf. change only the Database path. Put absolute path. eg template below for GIB. we had no isoseq here, but only trinity transcripts
``
DATABASE=annot.sqlite3
``
Launch_PASA_pipeline.pl -c alignAssembly.conf -C -R --ALIGNER gmap -g genome.fa -t transcripts.fasta
#required files will be : genome.sqlite3.pasa_assemblies*

#4. Mask the genome:
#abinitio repeats (-pa => number of threads)
BuildDatabase -name genome genome.fasta
RepeatModeler -database genome -pa 96 -engine ncbi
#library based -species will take genus, family etc also
queryRepeatDatabase.pl -species family > family.fasta
cat family.fasta *-families.fa > repeat_classes.fasta
#Masking
RepeatMasker -pa 96 -lib repeat_classes.fasta -s -xsmall -gff -excln genome.fasta

#5. Augustus and Genemark training
#with busco train augustus (much faster and easier results)
python /apps/busco/scripts/run_BUSCO.py -i genome.fasta.masked -o scaffolds_long -l /apps/busco/lineage -m genome -c 96 --long -z #(using BUSCO, copy retraining folder to Augustus config directory:augustus/config/species)
augustus --gff3=on --species=scaffolds_long genome.fasta.masked > scaffolds.augustus.gff3 

#with proteins from orthodb train genemark (takes maximum 64 threads)
##first Cluster using mmseqs2 at 90% coverage and 90% identity:
mmseqs createdb orthos.aa proteins_mmseqs
mmseqs cluster --cov-mode 1 -c 0.9 --min-seq-id 0.9 proteins_mmseqs cluster $PWD
mmseqs result2flat proteins_mmseqs proteins_mmseqs cluster cluster_out --use-fasta-header
grep "^>" cluster_out | sed 's/>//g' | seqtk subseq orthos.aa - | seqkit seq -m 10 > uniq_90_proteins.fasta
mkdir proteins && mv orthos.aa uniq_90_proteins.fasta proteins
##predict now
prothint.py --threads 64 scaffolds.fasta uniq_90_proteins.fasta

#predict augustus with BRAKER (much slower and might be slightly accurate than busco, prefer busco way as its faster and generally reliable. use this if busco doesn't feel good) (prothint_augustus.gff from Prothint, if using PASA replace CDNA_match with CDSPart)
nohup BRAKER_v2.1.0/braker.pl --genome=genome.fasta.masked --hints=prothint_augustus.gff --softmasking --gff3 --AUGUSTUS_CONFIG_PATH=Augustus/config/ --AUGUSTUS_BIN_PATH=Augustus/bin --AUGUSTUS_SCRIPTS_PATH=Augutus/scripts --cores=64 --AUGUSTUS_ab_initio --species=genome_braker --GENEMARK_PATH=GeneMarkES/gmes_linux_64/ --hints=genome.sqlite3.pasa_assemblies.gff3 --geneMarkGtf=genemark/genemark.gtf --bam=aligned.bam

#6. Protein alignments with Genomethreader (using internal scripts from BRAKER for parallel processing)
BRAKER_v2.1.0/startAlign.pl --genome=genome.fa --prot=uniq_90_proteins.fasta --prg=gth --CPU=64

#7. Accumulate evidences for EVM and run it
mkdir evm
#syslink/symlink/copy all the gff3 files from above to this folder including PASA and repeatmasker gff
#convert to EVM compatible format
EVidenceModeler/EvmUtils/misc/genomeThreader_to_evm_gff3.pl gth.aln > gth.gff3 #(or whatever is genomethreader's output) 
EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl augustus.gtf > augustus.gff3 #(use augustus_GFF... if your output is gff from augustus)
gt gtf_to_gff3 -o Genemark.gff3 genemark.gtf
cat Genemark.gff3 augustus.gff3 > Abinits.gff3
# create weights
EVidenceModeler/EvmUtils/create_weights_file.pl -A Abinits.gff3 -P gth.gff3 -T genome.sqlite3.pasa_assemblies.gff3 > weights.txt
#give largest weights to PASA, moderate/lower weights to Genemark & gth and lowest to augustus. will have to play around for choosing the best predicted set
#run EVM
EVidenceModeler//EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions Abinits.gff3 --transcript_alignments genome.sqlite3.pasa_assemblies.gff3 --protein_alignments gth.gff3 --segmentSize 1000000 --overlapSize 10000 --partition_listing partitions_list.out --repeats repeats.gff
EVidenceModeler//EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights `pwd`/weights.txt --gene_predictions Abinits.gff3 --protein_alignments gth.gff3 --repeats repeats.gff --transcript_alignments genome.sqlite3.pasa_assemblies.gff3 --output_file_name evm.out  --partition
s partitions_list.out >  commands.list
cat commands.list | parallel -j200
EVidenceModeler//EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
EVidenceModeler//EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome genome.fasta
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.gff3
EVidenceModeler/EvmUtils/gff3_file_fix_CDS_phases.pl EVM.gff3 genome.fasta.masked > fixed.EVM.gff3
python2.7 genomeannotation-GAG-997e384/gag.py -f genome.fasta -g fixed.EVM.gff3

#8. Functional annotation
#submit proteins from gag_output to eggNOG mapper and pannzer online
#interpro submission script is with suresh sir

#9. Tertiary analyses shared with suresh sir

