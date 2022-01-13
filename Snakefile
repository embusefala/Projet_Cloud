configfile: "config/config.yaml"

##### load rules #####
rule all:
	input:
		expand("results/fastqc_init/{sample}_fastqc.zip", sample = config["samples"]),
		expand("results/fastqc_init/{sample}_fastqc.html", sample = config["samples"]),
		expand("results/cutadapt/{samplename}1_cleaned.fastq.gz", sample = config["samples"]),
		expand("results/cutadapt/{samplename}2_cleaned.fastq.gz", sample = config["samples"])

rule unzip:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		temp("tmp/{sample}.fastq")
	shell:
		"""
		mkdir -p tmp
		gunzip -c {input} > {output}
		"""

rule fastqc_init:
	input:
		"tmp/{sample}.fastq"
	output:
		html="results/fastqc_init/{sample}_fastqc.html",
		zip="results/fastqc_init/{sample}_fastqc.zip"
	conda:
		"env.yaml"
	threads: 2
	shell:
		"""
		mkdir -p results/fastqc_init
		fastqc {input} -o "results/fastqc_init" -t {threads}
		"""

rule cutadapt:
	input:
		lambda wildcards: config["name_sample"==*1.fastq.gz][wildcards.sample],
		lambda wildcards: config["name_sample"==*2.fastq.gz][wildcards.sample]
	output:
		fastq1="results/cutadapt/{sample}1_cleaned.fastq.gz",
		fastq2="results/cutadapt/{sample}2_cleaned.fastq.gz"
	conda:
		"env.yaml"
	threads: 6

	shell:
		"""
		mkdir -p results/cutadapt
		cutadapt -j 1 -a R1=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A R2=CTGTCTCTTATACACATCTGACGCTGCCGACGA \
		--output={output} --paired-output={output} --error-rate0.1 --times=1 --overlap=3 --minimum-length=20 --pair-filter=any --quality-cutoff=20 \
		-j {threads} {input}
		"""

rule bowtie2 alignement:
	input:
		genome="data/mydatalocal/fasta/all.fasta"
		r1="results/cutadapt/{sample}1_*"
		r2="results/cutadapt/{sample}2_*"
	output:
		"results/bowtie2/{sample}_cleaned_mapped_sorted_q2.bam"

	conda:
		"env.yaml"
	threads: 8
	shell:
		"""
		mkdir -p results/bowtie2
		bowtie2 --very-sensitive -p {threads} -k 10 -x {input.genome} -1 {input.r1} -2 \
		| samtools view -q 2 -bs - | samtools sort - -o {output}
		"""

rule samtools Indexing mapped sorted bam file:
	input:
		index="data/mydatalocal/bowtie2"
		bam="results/bowtie2/{sample}_cleaned_mapped_sorted_q2.bam"
	output:


	conda:
		"env.yaml"
	shell:
		"""
		samtools index -b {input.bam}
		"""
rule samtools Statistics on mapping:
	input:
		bam="results/bowtie2/{sample}_cleaned_mapped_sorted_q2.bam"
	output:
		"results/bowtie2/{sample}.log"
	conda:
		"env.yaml"
	shell:
		"""
		samtools idxstats {input} > {output}
		"""
rule picard:
	input:
	output:
		"results/picard/{sample}.nondup.bam"
	conda:
		"env.yaml"
	threads: 6
	shell:
		"""
		java --
		"""
rule deeptools:
	input:
	output:
	conda:
		"env.yaml"
	threads:
	shell:

rule MACS2:
	input:
	output:
	conda:
		"env.yaml"
	threads:
	shell:

rule bedtools:
	input:
	output:
	conda:
		"env.yaml"
	threads:
	shell: