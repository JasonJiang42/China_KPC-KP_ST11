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
