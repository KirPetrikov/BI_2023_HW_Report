### Data quality control

```bash
fastqc raw_data/amp_res_1.fastq.gz raw_data/amp_res_2.fastq.gz
```

### Raw reads trimming

```bash
trimmomatic PE -phred33 \
raw_data/amp_res_1.fastq.gz \
raw_data/amp_res_2.fastq.gz trim_reads_1/Reads1_F_p.fq.gz \
trim_reads_1/Reads1_F_u.fq.gz trim_reads_1/Reads1_R_p.fq.gz \
trim_reads_1/Reads1_R_u.fq.gz SLIDINGWINDOW:10:20 LEADING:20 TRAILING:20 MINLEN:20
```

**Output**

Input Read Pairs: 455876 Both Surviving: 445689 (97.77%) Forward Only Surviving: 9758 (2.14%) Reverse Only Surviving: 284 (0.06%) Dropped: 145 (0.03%)

```bash
code
```
