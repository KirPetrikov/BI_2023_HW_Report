# BI_2023_HW_2_Report
> *This is the repo for the 2th homework of the BI-2023 Practicum Course*

Supplementary files for project:

### Title

**Files**:
- `Snakefile`  - rules file  for Snakemake workflow to reproduce pipeline
- `environment.yml` - file for setting up the conda/mamba virtual environment
- `vcf_parser.py` - Python script to extract values from vcf-file columns `REF`, `POS`, `ALT,` and `FORMAT:FREQ` to tsv-file.

**The script can be used standalone for any vcf.**

It is designed to be run from the command line. You must specify path to input file as first arg, optional: path to output file as second arg.

Default output filename: [INPUT_FILENAME]_selected.tsv.

*Example*
```bash
>>> python vcf_parser.py exmpl.vcf
>>> cat exmpl_selected.tsv
REF	COORD	ALT	FREQ,%
T	1458	C	0.84
```

### Instruction

**Don't forget to download files:**

**snakefile's folder must contain `{reference}.fna` - reference sequence, `{sample}.fastq` - sample reads.**

- Clone repo
```bash
git clone git@github.com:KirPetrikov/BI_2023_HW_Report.git
```

- Create new environment `RareSNP` or specified any name you want
```bash
mamba env create -f environment.yml
#
mamba env create -f environment.yml -n [ENV_NAME]
```

- To run Snakemake workflow specify output tsv-file as `{reference}.{sample}.tsv`. 
```bash
snakemake --cores=all -p reference_HA.sample_1.tsv
```

For `samtools mpileup` parameter `--max-depth` is set up to 50'000.

For `VarScan mpileup2snp` parameter `--min-var-freq` is set up to 0.001.


