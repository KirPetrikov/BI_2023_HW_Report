## Homework 7 Report
> *This is the repo for the homeworks of the BI-2023 Practicum Course*

Supplementary materials for project:

## Metagenomic analysis of the ancient oral microbiome from the oral cavity of fossil human remains

***by Dmitriy  Matach and Kirill Petrikov***

### Content

[**16S**](#16S)

[**MAG**](#MAG)

[**Terminal commands**](#Terminal-commands-used-for-MAG-analisys)

&nbsp;  
#### 16S
Additional files for 16S rRNA amplicone metagenome sequence:

- `teeth.Rmd` - `R Markdown`-script for `DADA2`-analisys
- `reads_quality_control.png` - `DADA2` quality control results for each sample
- `reads_quality_control_aggregate.png` - aggregated `DADA2` quality control results
- `asv_table.csv` - final ASV-tabel from `DADA2`-analisys
- `tax_table.csv` - ASV taxonomy assignment for `MicrobiomeAnalyst`
- `metadata.csv` - samples metadata for `MicrobiomeAnalyst`
- `norm_libsizes_0.png` - `MicrobiomeAnalyst` library size plot for input samples

Additional R-script for visualization of ASVs corresponding to "red complex" bacteria

```R
library(dplyr)
library(reshape2)

barplotdata <- ps@otu_table[, c(150, 135, 243, 129, 177, 327, 397)]
bpdata <- melt(barplotdata)

bpdata <- bpdata %>%
  mutate(affliction = case_when(
    Var1 == 'SRR986773.fastq' ~ "Peridontitis and the red complex",
    Var1 == 'SRR986774.fastq' ~ "Peridontitis and the red complex",
    Var1 == 'SRR986778.fastq' ~ "Peridontitis",
    Var1 == 'SRR986779.fastq' ~ "Peridontitis",
    Var1 == 'SRR986782.fastq' ~ "Peridontitis",
    .default = 'Unafflicted')
    )

p <- ggplot(bpdata, aes(x = Var1, y = Var2)) + 
  geom_point(aes(size = value, col = affliction)) + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  xlab('read') + 
  ylab('ASV')

p
```

&nbsp;  
#### MAG
Additional files for metagenome assembled genomes analisys:
- `sankey_plot.html` - Sankey diagram by `Pavian` for MAGs taxonomy assignment results
- `ref_unique_genes.gff` - Features selected after intersection as unique for *Tannerella forsythia* modern reference genome
- `KofamKOALA_result.txt` - results of `KofamKOALA` search on selected proteins unique for reference genome

&nbsp;  
#### Terminal commands used for MAG analisys
Create index of referense genome, MAGs contigs alignment
```bash
bwa index NC_016610.fasta && \
bwa mem -t 8 NC_016610.fasta G12_assembly.fna | \
samtools view -b --threads 8 | \
samtools sort --threads 8 > alignment.sorted.bam
```

Get alignment statistics
```bash
samtools flagstat alignment.sorted.bam
```

Create `bed`-file from `bam`-file
```bash
bedtools bamtobed -i alignment.sorted.bam > alignment.bed
```

Intersect referense genome and MAGs contigs, keeping only unique for referense genome
```bash
bedtools intersect -v -a NC_016610.gff3 -b alignment.bed > intersect.gff
```

Select only `CDS`, filter out «hypothetical» and «pseudo» proteins, keep only transposases and all but transposases
```bash
awk '{FS="\t";OFS="\t"} $3 ~ "CDS"' intersect.gff | grep -v 'product=hypothetical' | grep -v 'pseudo=true' | grep 'transposase' > transp.gff

awk '{FS="\t";OFS="\t"} $3 ~ "CDS"' intersect.gff | grep -v 'product=hypothetical' | grep -v 'pseudo=true' | grep -v 'transposase' > cds_inters.gff
```

Select proteins IDs for further `NCBI Entrez` batch download and `KofamKOALA`-analisys
```bash
cut -f 9  cds_inters.gff | cut -f 1 -d ";" | cut -f 2 -d "=" | grep -o -P "WP_\d+.\d" > prot_idxs.txt
```
