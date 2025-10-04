# KPC-producing *K. pneumoniae* ST11 genomic analysis in China #   
This repository outlines the analysis pipeline for the paper *Z. Ma. et al. Emergence and regional evolution of Klebsiella pneumoniae carbapenemase (KPC)-producing Klebsiella pneumoniae ST11 in China (BMC Genomics)*

## Genome assembly and QC assessment ##
High-quality reads were assembled using SPAdes (https://github.com/ablab/spades)  
```
python spades.py --pe1-1 file1 --pe1-2 file2 -o assmebly 
```
Species confirmation by GTDB-Tk (https://github.com/Ecogenomics/GTDBTk)  
```
gtdbtk classify_wf --genome_dir genomes --out_dir gtdbtk/classify --cpus 10 --skip_ani_screen
```
QC assessment by checkM (https://github.com/Ecogenomics/CheckM)  
```
checkm lineage_wf -x fasta input_bins output_folder
```
Sequence statistics analysis
```
seqkit stats -a input.fasta > stats.tsv
```
## Genome annotation ##
The MLST and capsule type were annotated using Kleborate v2.3.1 (https://github.com/klebgenomics/Kleborate)
```
kleborate --all -a input.fasta -o kaptive.results.tsv
```
The ARG and VF were predicted using Abricate (https://github.com/tseemann/abricate). The databases used for annotation are ResFinder (https://github.com/cadms/resfinder) and VFDB (https://www.mgc.ac.cn/VFs/), respectively.
```
abricate --db resfinder/VF input.fasta > results.tsv
```
The plasmid type of our genomes was determined using KleTy (https://github.com/zheminzhou/KleTy).
```
KleTy.py -q input.fasta -o predix -n 8 -g
```
## Phylogenetic analysis ##
Core genome alignment was generated using snippy (https://github.com/tseemann/snippy), and recombination sites were removed with ClonalFrameML (https://github.com/xavierdidelot/ClonalFrameML). A maximum-likelihood phylogenetic tree was then constructed using IQ-TREE (http://www.iqtree.org/) based on clean core genome SNP alignments.
```
snippy --outdir mut1 --ref ref.gbk --ctgs mut1.fasta
iqtree -s core.aln --boot-trees --wbtl -m GTR+I+G -B 1000 -nt 18
ClonalFrameML tree_file recom_filter
iqtree -s clonalframe.ML_sequence.fasta --boot-trees --wbtl -m GTR+I+G -B 1000 -nt 16
```
## Genomic population ##
The core gene allele calling was performed using chewBBACA (https://github.com/B-UMMI/chewBBACA)
```
## create schema
chewBBACA.py CreateSchema -i InputAssembliesFolder -o KPC-ST11_schema --ptf KPC-ST11

## perform allele calling
chewBBACA.py AlleleCall -i genomes/ -g KPC-ST11_schema/schema_seed -o KPC-ST11_wgMLST --cpu 6

## Paralog detection to determine if some of the loci can be considered paralogs on the wgMLST allele calling.
chewBBACA.py RemoveGenes -i KPC-ST11_wgMLST/results_alleles.tsv -g KPC-ST11_wgMLST/paralogous_counts.tsv -o KPC-ST11_wgMLST/results_alleles_NoParalogs.tsv

## cgMLST schema determination
chewBBACA.py ExtractCgMLST -i KPC-ST11_wgMLST/results_alleles_NoParalogs.tsv -o KPC-ST11_wgMLST/cgMLST
```
All cgMLST loci were then grouped into multi-level clusters using pHierCC (https://github.com/zheminzhou/pHierCC)
```
pHierCC -p KPC-ST11_cgMLST -o KPC-ST11.cgMLST.HierCC
HCCeval -p KPC-ST11.cgMLST -c KPC-ST11.cgMLST.HierCC.HierCC -o KPC-ST11.cgMLST.HierCC.eval
```
The minimum spanning trees were constructed using GrapeTree (https://github.com/achtman-lab/GrapeTree)
```
python grapetree.py -p cgMLST.profile -m MSTreeV2
```
## Bayesian inference
Bayesian dating of the phylogenetic tree after recombination analysis was performed by BactDating (https://github.com/xavierdidelot/BactDating)
```
library(BactDating)
library(ape)
t=loadCFML("clonalframe_output")
dates=c(date1,date2...)
tree=bactdate(t, dates, useRec=T)
```
Phylodynamic inference of effective population size was conducted using Skygrowth (https://github.com/mrc-ide/skygrowth)
library(skygrowth)
library(ape)
skygrowth.map(tree_file, res=24*13, tau0=.1)
```
Geographical transmission patterns were estimated using TreeTime (https://treetime.readthedocs.io/en/latest/tutorials/mugration.html)
```
treetime mugration --tree tree_file --states metadata.tsv --attribute country
```

