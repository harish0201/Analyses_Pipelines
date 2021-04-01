##### Gene loss/gain evolution analyses #### (CAFE)

#### 1. Get protein sequences of the species of interest including yours (prefer ensembl for other species)

#### 2. Run orthofinder to get two things: 
	## 1. Orthogroups.GeneCount.tsv
	## 2. SpeciesTree_rooted.txt

orthofinder -t 192 -a 192 -M msa -S diamond -A mafft -T fasttree -f $DIR ##($DIR=folder with protein sequences)

#### 3. Convert SpeciesTree_rooted.txt to ultrametric. Find time of outgroup or root node from TimeTree.org: age-> of root node
python2.7 OrthoFinder/tools/make_ultrametric.py -r age SpeciesTree_rooted.txt

#### 4. Tips:
	## 1. Check the cafe tutorial pdf
	## 2. Estimate using multi-lambda-mu model always
	## 3. Round off the branch lengths in ultrametric tree so that it matches the total age inputted from timetree to create ultrametric tree
	## 4. iTOL website is your best help to visualize and check

#### 5. Analyses:
awk -F'\t' '{print "(null)\t"$0}' Orthogroups.GeneCount.tsv > tmp.tsv
``
Orthogroup      Bombyx_mori     Chilo_suppressalis      Drosophila_melanogaster Helicoverpa_armigera    Leucinodes_orbonalis    Manduca_sexta   Papilio_machaon Plutella_xylostella     Spodoptera_litura       Tribolium_castaneum     Trichoplusia_ni        Total
OG0000000       2       7       0       9       450     0       0       0       0       1       0       469
OG0000001       0       0       0       0       459     0       1       0       1       0       0       461
OG0000002       0       3       0       2       397     0       0       0       0       0       0       402
OG0000003       0       0       0       1       349     0       0       0       0       0       1       351
OG0000004       50      0
``
#Change the header (null) to Desc and save
``
Desc    Orthogroup      Bombyx  Chilo   Drosophila      Helicoverpa     Leucinodes      Manduca Papilio Plutella        Spodoptera      Tribolium       Trichoplusia
(null)  OG0000064       18      5       9       2       1       16      3       12      11      6       16
(null)  OG0000065       0       5       0       15      74      0       0       5       0       0       0
(null)  OG0000066       22      1       13      1       1       3       7       1       9       20      20
(null)  OG0000067       0       3       0       0       93      0       2       0       0       0       0
(null)  OG0000068       5       4       0       4       33      0       1       24      21      0       5
(null)  OG0000069       12      9       7       8       8       11      10      8       8
``

#run size filter finally
#filter the Orthogroups.GeneCount.tsv file to remove OG that have more than 100 proteins in a particular species:
python2.7 python_scripts/cafetutorial_clade_and_size_filter.py -i mod.tsv -s -o cafe.input.tsv

# Prepare lamda tree structure: will look like this:
#(6,(5,(4,(3,((2,2)2,((1,1)1,((,),1)1)2)3)4)5)6)
# run cafe: inputs the tsv output from above; and using the rounded off tree.nwk

#prepare the report file: (input is the cafe output)
##### inputs for -d: (last line in the report.txt generated below: will look something liek: ((Aedes<0>,Drosophila<2>)<1>,((Galleria<4>,(Amyelois<6>,Plodia<8>)<7>)<5>,(((Chilo<10>,Myelobia<12>)<11>,(Falsella<14>,Oregonicus<16>)<15>)<13>,(Spilomela<18>,(Cirrhochrista<20>,Maruca<22>)<21>)<19>)<17>)<9>)<3>
##### inputs for -t: rounded off ultrametric tree
##### enclose -d and -t in single quotes ''
python2.7 python_scripts/cafetutorial_report_analysis.py -r 0 -o report -i try2.report.cafe > report.txt
python2.7 python_scripts/cafetutorial_draw_tree.py -i report_node.txt -o exansions.png -y Expansions -t -d 
