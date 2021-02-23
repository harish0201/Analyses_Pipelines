####Segmental duplications. Find a better pipeline, I think mashmap can be used according to their paper. This is with asgart, works faster and has integrated plotting approach.

#1. Asgart lookup
asgart --threads 96 -RCv --compute-score genome.fasta

#2. Asgard convert to gff
asgart-slice -f gff3 -F --min-length 1000 -o genome.gff3 genome.json
grep -v "COLLAPSED" genome.gff3 | grep -v "^#" | awk -F'\t' '{if($5-$4>=1000 && $6>=90)print $0}' > 1000_90.uncollapsed.gff3
grep "COLLAPSED" genome.gff3 | grep -v "^#" | awk -F'\t' '{if($5-$4>=1000 && $6>=90)print $0}' > 1000_90.collapsed.gff3

#3. For plotting: needs to explore. check asgart-plot
asgart-plot 1000_score.json chord --min-identity 90 -C by-type --min-length 1000
