## Homework 5 Report
> *This is the repo for the homeworks of the BI-2023 Practicum Course*

Supplementary materials for project:

## Title
***by Dmitriy  Matach and Kirill Petrikov***

Here you can find:
- scripts to reproduse analysis
- detailed intermediate and final results

### Instruction to reproduse analysis

- Clone repo
```bash
git@github.com:KirPetrikov/BI_2023_HW_Report.git && \
cd BI_2023_HW_Report/
```

- Create new environment `Project_yeast`
```bash
mamba activate && \
mamba env create -f environment.yml
```

- Run `de_analysis.sh`
```
bash de_analysis.sh
```

### Repo content

**Root directory files**:
- `environment.yml` - file for setting up the conda/mamba virtual environment

- `de_analysis.sh` - script to reproduse analysis from files downloads to DESeq2 results

    *If you want to create reads quality control reports just uncomment line 14.*

    Files required to run the script:

  - `links.txt` - list of download links

  - `deseq2.r`, `draw-heatmap.r` - `R`-scripts to run `DESeq2` and draw a heatmap

&nbsp;  
**`data` directory files**
- `hisat.log` - log of `HISAT2` alignment
- `featureCounts_results.txt`, `featureCounts.log` - results and log of 'featureCounts' processing
- `count_table.txt` - count table for `DESeq2`
- `DESeq2_results.tsv` - list of genes with assessment of their differential expression by `DESeq2`

&nbsp;  
**`QC` directory files**
- quality control results by `FastQC` and `MultiQC`
