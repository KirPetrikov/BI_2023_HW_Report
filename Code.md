
1. [Raw reads quality control](#raw-reads-quality-control)
2. [Raw reads trimming](#raw-reads-trimming)
3. [Trimmed reads quality control](#trimmed-reads-quality-control)
4. [Indexing, alignment](#indexing-alignment)
5. [Alignment preprocessing, variant calling](#alignment-preprocessing-variant-calling)
6. [Variant annotation](#variant-annotation)

### Raw reads quality control

```bash
fastqc raw_data/amp_res_1.fastq.gz raw_data/amp_res_2.fastq.gz
zcat raw_data/amp_res_1.fastq.gz | wc -l
zcat raw_data/amp_res_2.fastq.gz | wc -l
```

Output:
```bash
1823504
1823504
```

### Raw reads trimming

With threshold 20

```bash
trimmomatic PE -phred33 \
raw_data/amp_res_1.fastq.gz \
raw_data/amp_res_2.fastq.gz trim_reads_1/Reads1_F_p.fq.gz \
trim_reads_1/Reads1_F_u.fq.gz trim_reads_1/Reads1_R_p.fq.gz \
trim_reads_1/Reads1_R_u.fq.gz SLIDINGWINDOW:10:20 LEADING:20 TRAILING:20 MINLEN:20
```

Output:
```
Input Read Pairs: 455876
Both Surviving: 445689 (97.77%)
Forward Only Surviving: 9758 (2.14%)
Reverse Only Surviving: 284 (0.06%)
Dropped: 145 (0.03%)
```

With threshold 30

```bash
trimmomatic PE -phred33 \
raw_data/amp_res_1.fastq.gz raw_data/amp_res_2.fastq.gz \
trim_reads_2/Reads2_F_p.fq.gz trim_reads_2/Reads2_F_u.fq.gz \
trim_reads_2/Reads2_R_p.fq.gz trim_reads_2/Reads2_R_u.fq.gz \
SLIDINGWINDOW:10:30 LEADING:30 TRAILING:30 MINLEN:20
```

Output:

```
Input Read Pairs: 455876
Both Surviving: 370302 (81.23%)
Forward Only Surviving: 35688 (7.83%)
Reverse Only Surviving: 26523 (5.82%)
Dropped: 23363 (5.12%)
```

### Trimmed reads quality control

```bash
fastqc trim_reads_1/Reads1_F_p.fq.gz trim_reads_1/Reads1_R_p.fq.gz
fastqc trim_reads_2/Reads2_F_p.fq.gz trim_reads_2/Reads2_R_p.fq.gz
```

### Indexing, alignment

```bash
bwa index genome_ref.fna.gz

bwa mem -t 8 genome_ref.fna.gz \
trim_reads_1/Reads1_F_p.fq.gz \
trim_reads_1/Reads1_R_p.fq.gz > alignment.sam

samtools view -@ 8 -b alignment.sam > alignment.bam
```

Alignment statistic

```bash
samtools flagstat alignment.bam
```
Output:

```bash
891635 + 0 in total (QC-passed reads + QC-failed reads)
891378 + 0 primary
0 + 0 secondary
257 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
890570 + 0 mapped (99.88% : N/A)
890313 + 0 primary mapped (99.88% : N/A)
891378 + 0 paired in sequencing
445689 + 0 read1
445689 + 0 read2
887614 + 0 properly paired (99.58% : N/A)
889386 + 0 with itself and mate mapped
927 + 0 singletons (0.10% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Alignment preprocessing, variant calling

```bash

samtools sort -@ 8 alignment.bam -o alignment.sort.bam
samtools index alignment.sort.bam

samtools mpileup -f genome_ref.fna alignment.sort.bam > alignment.mpileup

# threshold 50%
java -jar ~/opt/VarScan/VarScan.v2.3.9.jar mpileup2snp alignment.mpileup \
--min-var-freq 0.50 --output-vcf 1 --variants > VarScan_results_50.vcf

# threshold 70%
java -jar ~/opt/VarScan/VarScan.v2.3.9.jar mpileup2snp alignment.mpileup \
--min-var-freq 0.70 --output-vcf 1 --variants > VarScan_results_70.vcf
```

### Variant annotation

```bash
# database building
java -Xmx8g -jar ~/opt/snpEff/snpEff.jar \
build -c ~/Main/IB/Projects/1/snpEff_db/snpEff.config -v k12

# annotation
java -Xmx8g -jar ~/opt/snpEff/snpEff.jar \
ann -c ~/Main/IB/Projects/1/snpEff_db/snpEff.config -v \
k12 ~/Main/IB/Projects/1/VarScan_results_70.vcf >VarScan_results_ann.vcf
```

```bash
code
```
