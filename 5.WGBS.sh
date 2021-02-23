#####WGBS Analyses and assumes QC Passed data: explore wg-blimp tool as it does all the below. it is a new tool.
###### CGMAPtools (https://cgmaptools.github.io/cgmaptools_documentation/what-is-cgmaptools.html) and defiant (https://github.com/hhg7/defiant)
##### change parameters as needed. WGBS is really funky.

#### 0. Index the genome with bismark
bismark_genome_preparation --parallel 4 --path_to_aligner /apps_940/hisat2*/ --hisat2 --genomic_composition <path_to_genome_folder>/genome.fasta

#### 1. Align with Bismark (will be present on /apps_940/WGBS) (will use 20*4=80 threads)
bismark --hisat2 --known-splicesite-infile splice.hisat2 -p 20 --nucleotide_coverage -B sample1 --bam --gzip --parallel 4 --path_to_hisat2 /apps_940/hisat2*/ <path_to_genome_folder> -1 sample1_r1.fq.gz -2 sample2_r2.fq.gz

#### 2. Deduplicate bismark alignments
deduplicate_bismark --bam -o sample1.dedupe.bam sample1.bam #(replace with whatever bismark bam is)

#### 3. Extract methylation:
bismark_methylation_extractor --genome_folder <path_to_genome_folder> --CX_context --cytosine_report --ample_memory --parallel 10 --gzip --no_header --merge_non_CpG --comprehensive sample1.dedupe.bam

#### 4. From here will use CGMAPtools (defiant is a good alternative)
cgmaptools convert bam2cgmap -b sample1.dedupe.bam -g genome.fasta --rmOverlap -o sample.cgmap

#### 5. Merge replicates (ignore if no replicates and skip to 6)
cgmaptools merge2 cgmap -1 sample1_1.cgmap -2 sample1_2.cgmap
#if more than 2 samples then prepare a list of files
ls sample1_*.cgmap > sample1.list
cgmaptools mergelist tosingle -i sample1.list -o sample1.merged.cgmap

#### 6. Intersect the comparisons (use multiple context if needed: replace CG with CH, CHG, CHH, CA, CC, CT, CW etc)
cgmaptools intersect -1 sample1.merged.cgmap -2 sample2.merged.cgmap -o sample1_2.intersected.cgmap -C CG

#### 7. Call DMRs:
cgmaptools dmr -i sample1_2.intersected.cgmap -o sample1_2.dmr.gz
Rscript DMR_padj_cgmap.R sample1_2.dmr.gz

#### 8. Call DMSs:
cgmaptools dms -i sample1_2.intersected.cgmap -o sample1_2.dms.gz
Rscript DMS_padj_cgmap.R sample1_2.dms.gz

#### 9. Heatmap plots and others:
#1. heatmap
cgmaptools mmbin -l sample1.CGmap,sample2.CGmap,sample3.CGmap > mmbin.tab
cgmaptools heatmap -i mmbin.tab -c -o cluster.pdf -f pdf

#2. Genebody: (needs ref-flat file, please check the internet)
cgmaptools bed2fragreg -i genome.annot.bed -F 50,50,50 -T 50,50,50 -n 1 -o fragreg.bed
cgmaptools mfg -i sample1.cgmap -r fragreg.bed -c 2 -x CG > S1.mfg
cgmaptools mfg -i sample2.cgmap -r fragreg.bed -c 2 -x CG > S2.mfg
(head -1 S1.mfg | gawk '{$1="Sample"; print $0;}'; for F in *.mfg; do gawk -vSampleName=`echo $F | sed s/.mfg//g` '/total_ave_mC/{$1=SampleName; print $0;}'  done) > mfg_merge.xls
cgmaptools fragreg -i mfg_merge.xls -o merge.fragreg.pdf -f pdf

#3. OVerall coverage:
cgmaptools oac bin -i sample.cgmap -B 1000 -f png -p sample -t sample > sample.oac_bin.data

#4. Methylation effective coverage:
cgmaptools mec stat -i sample.cgmap -p sample -f png > sample.mec_stat.data


#### 10. Gene identification with DMRs and Ontology/Pathway analyses
#filter for regions/bases with Padj < 0.05 from step 7.
#bedtools merge can be used for merging consecutive regions
#bedtools intersect to intersect with the refflat file from #### 9/#2.
#follow GO and KEGG pathway from RNASeq steps
