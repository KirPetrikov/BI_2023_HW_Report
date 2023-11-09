
rule bwa_index:
    input:
        "{reference}.fna"
    output:
        "{reference}.fna.amb",
        "{reference}.fna.ann",
        "{reference}.fna.bwt",
        "{reference}.fna.pac",
        "{reference}.fna.sa"
    shell:
        "bwa index {input}"

rule bwa_align:
    input:
        rules.bwa_index.output,
        ref="{reference}.fna",
        reads="{sample}.fastq"
    threads: 8
    log:
        "logs/bwa.{reference}.{sample}.log"
    output:
        "{reference}.{sample}.bam"
    shell:
        "bwa mem -t {threads} {input.ref} {input.reads} 2>{log} | samtools view -b -o {output}"

rule bam_sort:
    input:
        rules.bwa_align.output
    threads: 8
    output:
        "{reference}.{sample}.sort.bam"
    shell:
        "samtools sort --threads {threads} {input} -o {output}"

rule bam_idx:
    input:
        rules.bam_sort.output
    threads: 8
    output:
        "{reference}.{sample}.sort.bai"
    shell:
        "samtools index -@ {threads} {input} -o {output}"

rule mpileup:
    input:
        rules.bam_idx.output,
        ref="{reference}.fna",
        bam="{reference}.{sample}.sort.bam"
    output:
        "{reference}.{sample}.mpileup"
    shell:
        "samtools mpileup -d 50000 -f {input.ref} {input.bam} -o {output}"

rule varscan:
    input:
        rules.mpileup.output
    log:
        "logs/vcf.{reference}.{sample}.log"
    output:
        "{reference}.{sample}.vcf"
    shell:
        "java -jar ~/opt/VarScan/VarScan.v2.3.9.jar mpileup2snp {input} --min-var-freq 0.001 --variants --output-vcf 1 > {output} 2>{log}"

rule parse_vcf:
    input:
        rules.varscan.output
    output:
        "{reference}.{sample}.tsv"
    shell:
        "python vcf_parser.py {input} {output}"

