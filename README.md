## Homework 4 Report
> *This is the repo for the homeworks of the BI-2023 Practicum Course*

Supplementary files for project:

### Search for proteins with nuclear localization potentially responsible for DNA radiation protection in the tardigrade *Ramazzottius varieornatus* strain YOKOZUNA-1
***Dmitriy  Matach and Kirill Petrikov***

**Files**:
- `environment.yml` - file for setting up the conda/mamba virtual environment
- `Fasta_extract.sh` - script for fasta extraction

#### Instruction

- Clone repo
```bash
git@github.com:KirPetrikov/BI_2023_HW_Report.git
```

- Create new environment `Project`
```bash
mamba env create -f environment.yml
```

**Using `Fasta_extract.sh`**

This script will allow you to extract only selected records from the fasta-file.

You need  a txt-file containing the headers of the required fasta-entries, which is on a new line, without >. Run the script. A new file `new_[INPUT_FILE]` will be created with the selected records.

```bash
bash Fasta_extract.sh [INPUT_FILE]
```


