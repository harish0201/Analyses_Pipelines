#### Assumes QC and indexing is done###

#### 0. index a genome
bwa index genome.fasta
samtools faidx genome.fasta
samtools dict genome.fasta > genome.dict

#### 1. Align data to the reference: (assumes files are labelled as *_R1.clean.fq.gz and *_R2.clean.fq.gz: $1==path of clean data; $2==path to genome index; $3==path of output 
mkdir -p "$3"/sorted "$3"/disc "$3"/split && for i in "$1"/*_R1.clean.fq.gz; do echo "bwa mem -t 24 -R '@RG\tID:$(basename "$i" _R1.clean.fq.gz)\tSM:$(basename "$i" _R1.clean.fq.gz)\tLB:$(basename "$i" _R1.clean.fq.gz)' "$2" "$i" "$1"/$(basename "$i" _R1.clean.fq.gz)_R2.clean.fq.gz | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d "$3"/disc/$(basename "$i" _R1.clean.fq.gz).disc.sam -s "$3"/split/$(basename "$i" _R1.clean.fq.gz).split.sam -o /dev/stdout |sambamba view -S /dev/stdin -f bam -t 24 -o /dev/stdout | sambamba sort -m 3GB --tmpdir="$1"/tmp -o "$3"/sorted/$(basename "$i" _R1.clean.fq.gz).sorted.bam -l 9 -t 24 /dev/stdin"; done | parallel -j4
#generate csi index if the bam file is > 100Gb or if deepvariant etc complains
#samtools index -@ 100 sample.bam

#### 2. Variant calling via multiple tools (pick your tool of interest). For join calling prefer GATK or Deepvariant. For single samples you can use freebayes
	#a) deepvariant sample command: change input and output path; this will run singly for mutlple samples. cant get looped. gvcfs have to be joined
BIN_VERSION="0.10.0"
docker run -v "/Analysis/dr.joykumar/clean2/JKRB1/remap/sorted":"/input" -v "/Analysis/dr.joykumar/clean2/JKRB1/remap/sorted":"/output" google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/input/
Triticum_aestivum.IWGSC.dna.toplevel.fa.gz --reads=/input/JKRB1.sorted.bam --output_vcf=/output/JKRB1.vcf.gz  --output_gvcf=/output/JKRB1.g.vcf.gz --tmpdir=/input/tmp --num_shards=96 2>log
		#combine all the gvcf for joint calling
docker run -v "/Analysis/dr.joykumar/clean2/JKRB1/remap/sorted":"/data" quay.io/mlin/glnexus:v1.2.6  /usr/local/bin/glnexus_cli --config DeepVariantWGS *.g.vcf.gz |  | bcftools view - bgzip -c >  /Analysis/dr.joykumar/clean2/JKRB1/remap/sorted/deepvariant.cohort.vcf.gz


	#b) freebayes (need to filter with various approaches):
freebayes-parallel <(fasta_generate_regions.py ../genome.fasta.fai 1000000) 36 -f ../genome.fasta --ploidy 2 -m 30 -q 30 --no-mnps --no-complex --min-coverage 10 --haplotype-length 5 --strict-vcf $(ls *.bam) > Raw_variants.vcf
bcftools query -f '%DP\n' split_variants.vcf | datamash mean 1 > meanDP
MAX_Depth_filter = DP_av+3*sqrt(DP_av) or DP_av+4*sqrt(DP_av)
Qual thresh= 2*MAX_Depth_filter
bcftools filter -g 35 -G 35 -i -o filtered_variants.vcf -O z --soft-filter FBQualDepth -e '(AF[0] <= 0.5 && (FORMAT/DP < 4 || (FORMAT/DP < 13 && %QUAL < 10))) || (AF[0] > 0.5 && (FORMAT/DP < 4 && %QUAL < 50)) || (%QUAL < qual_thresh && FORMAT/DP > depth_thresh && AF[0] <= 0.5)' -m + Raw_variants.vcf

	#c) gatk (comes with filtering)
REFERENCE="$1"
name=$(echo $REFERENCE | sed 's/.fa.*\|.fna.*//g')
for i in *.bam; do echo "gatk4 --java-options "-Xmx40g" HaplotypeCaller -R $REFERENCE -I "$i" -O $(basename "$i" .bam).g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation"; done | parallel -j2
ls *g.vcf.gz > gvcf.list
awk '{print $1"\t0\t"$2-1}' $REFERENCE.fai > reference.bed
picard BedToIntervalList I=reference.bed SD=$name.dict O=interval.list
gatk4 GenomicsDBImport -V gvcf.list --genomicsdb-workspace-path my_database -L interval.list
for i in SNP INDEL; do echo "gatk4 SelectVariants -R $REFERENCE -V genotyped.vcf.gz --select-type-to-include "$i" -O "$i".vcf.gz"; done | parallel -j2
gatk4 VariantFiltration -R $REFERENCE -V SNP.vcf.gz --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O filtered_snps.vcf.gz
gatk4 VariantFiltration -R $REFERENCE -V INDEL.vcf.gz --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O filtered_indels.vcf.gz

gatk4 SortVcf -I filtered_snps.vcf.gz -O sort.SNP.filt.vcf.gz
gatk4 RemoveNearbyIndels -V filtered_indels.vcf.gz -O distance.INDEL.filt.vcf.gz
gatk4 SortVcf -I distance.INDEL.filt.vcf.gz -O sort.INDEL.filt.vcf.gz

ls sort*.vcf.gz | grep -v "tbi" > filtered.list
gatk4 MergeVcfs -R $REFERENCE -I filtered.list -O filtered.vcf.gz -D $name.dict
bcftools view -f PASS -O z -o passed.vcf.gz filtered.vcf.gz


##TIPS: When combining from multiple sources, prefer GATK as the primary, it will make filtering easier

#### 3. For SV calling:
mkdir outs
docker run -v $PWD:/home/dnanexus/in -v $PWD/outs:/home/dnanexus/out dnanexus/parliament2:latest --bam input.bam --bai input.bai -r ref_genome.fa.gz --breakdancer --cnvnator --breakseq --lumpy --delly_deletion --delly_insertion --delly_inversion  --delly_duplication --genotype --fai ref_genome.fai --filter_short_contigs --svviz_only_validated_candidates --svviz


#### 4. Annotation of variants:
#with snpeff
for i in *.vcf;do echo "snpEff ann -csvStats $(basename "$i" .vcf).csv.stats genome "$i" > $(basename "$i" .vcf).Annotated.vcf"; done | parallel -j2 &
#with VEP
bcftools view pass.snp.vcf.gz | vep --fork 4 --everything --species gallus_gallus --format vcf --cache --dir_cache /data/analysis/shared_files/Harish/apps/miniconda3/envs/variants/vep/  --fasta /data/ref_genomes/chicken/chicken.fasta -i STDIN --vcf --compress_output bgzip
