#!/bin/bash
# Script to reproduse analysis of yeast's genes differential expression

# ----- Download reads, reference genome with annotation
for VAR in $(cat links.txt)
do
wget ${VAR}
done
echo -e "\n ---- Downloads completed ----- \n"


# ----- If you want create reads quality control reports just uncomment next line:

# for VAR in $(ls | grep '^SRR'); do fastqc ${VAR} -o QC/; done && multiqc QC/ -o MultiQC/ && echo -e "\n ---- Quality control complete ----- \n"

# ----- Reporst will created in QC/ and MultiQC/ directiries


# ----- Unpack archives with genome and annotation
for VAR in $(ls | grep '^GCF')
do
gzip -d ${VAR}
done
echo -e "\n ---- Archives unpacked ----- \n"


# ----- Confert gff-annotation to gtf for featureCounts
gffread GCF_000146045.2_R64_genomic.gff -T -o GCF_000146045.2_R64_genomic.gtf
# ----- Remove two lines which provoked featureCounts error
cat GCF_000146045.2_R64_genomic.gtf | grep 'gene_id' > GCF_000146045.2_R64_genomic_ed.gtf


# ----- Alignment with HISAT2
# ----- Create index
hisat2-build *fna hisat_idx

# -----Align, log-file will created
for VAR in $(ls | grep '^SRR')
do
TITLE=${VAR::-9}
echo -e "${TITLE}" >>./Alig.log
hisat2 -p 16 -x hisat_idx -U ${VAR} 2>>hisat.log | samtools sort --threads 16 > ${TITLE}.sort.bam
echo -e "\n" >> hisat.log
done
echo -e "\n ---- Alignments completed ----- \n"


# ----- Create count table with featureCounts
# ----- Use featureCounts co count number of aligned reads per gene
featureCounts -g gene_id -a GCF_000146045.2_R64_genomic_ed.gtf -o featureCounts_results.txt SRR*.sort.bam 2>>featureCounts.log
# ----- Rename samples and select columns for further analysis
cat featureCounts_results.txt | grep -v "^#" | \
sed '1s/.*/Geneid\tChr\tStart\tEnd\tStrand\tLength\t0_min_rep_1\t0_min_rep_2\t30_min_rep_1\t30_min_rep_2/' | \
cut -f 1,7-10 > count_table.tsv
echo -e "\n ---- Count table created ----- \n"


# ----- Calculate differentialy expressed genes via DESeq2
cat count_table.tsv | R -f deseq2.r
# ---- Draw heatmap


cat norm-matrix-deseq2.txt | R -f draw-heatmap.r
echo -e "\n ---- Analysis completed ----- \n"

