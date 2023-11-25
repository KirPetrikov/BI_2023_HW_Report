## Homework 3 Report WIP
> *This is the repo for the homeworks of the BI-2023 Practicum Course*

Supplementary files for project:

#### A case of the emergence of a pathogenic E. coli strain caused by horizontal transfer of enterotoxin and antibiotic resistance genes

***Ariuna Aiusheeva and Kirill Petrikov***

**Files**:
- `environment.yml` - file for setting up the conda/mamba virtual environment with most of packages exept Quast
- `environment_q.yml` - file for setting up the conda/mamba virtual environment with Quast

#### Instruction

- Clone repo
```bash
git@github.com:KirPetrikov/BI_2023_HW_Report.git
```

- Create new environment `Project`
```bash
mamba env create -f environment.yml
```

You can try to install Quast in this environment.

In case of error you can create environment `QUAST` from second yaml-file:

```bash
mamba env create -f environment_q.yml
```
