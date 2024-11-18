from os.path import join
import yaml, sys, os, csv

snakemake_dir = '/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/'
configfile: join(snakemake_dir, 'config/config.yaml')

#make directories if they haven't already been found
step_dirs = ['bam', 'fastqc_out', 'trimmed_fastqc_out', 'junction_counts', 'kallisto_out', 'bam_stats', 'trimmed_fastq'] 

for step_dir in step_dirs:
    if os.path.exists(join(snakemake_dir, step_dir)) == False:
        os.mkdir(join(snakemake_dir, step_dir))

fastq_dir = config['fastq_dir']
fasta_file = config['fasta_file']

def scheduler_time_hours(hour):
    return hour*60

def scheduler_memory(memory_gb):
    return memory_gb*1000

splicing_workflow = config['splicing_workflow']
rule all:
    input:
        expand(join(snakemake_dir, "fastqc_out/{sample}_R1_fastqc.html"), sample = config['samples']),
        expand(join(snakemake_dir, "fastqc_out/{sample}_R2_fastqc.html"), sample = config['samples']),
        expand(join(snakemake_dir, "trimmed_fastqc_out/{sample}_trimmed_R1_fastqc.html"), sample = config['samples']),
        expand(join(snakemake_dir, "trimmed_fastqc_out/{sample}_trimmed_R2_fastqc.html"), sample = config['samples']),
        expand(join(snakemake_dir, "kallisto_out/{sample}/abundance.tsv"), sample = config['samples']),
        expand(join(snakemake_dir, "bam_stats/{sample}_stats.txt"), sample = config['samples']),
        expand(join(snakemake_dir, "bam_stats/{sample}.metrics.tsv"), sample = config['samples'])
        if splicing_workflow:
            expand(join(snakemake_dir, "junction_counts/{sample}_junctions.out"), sample = config['samples'])

rule trim_adapter:
    input:
        R1 = join(snakemake_dir, "fastq/{sample}_R1.fastq.gz"),
        R2 = join(snakemake_dir, "fastq/{sample}_R2.fastq.gz")
    output:
        R1 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"),
        R2 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz") 
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 40 -j {threads} -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        """

rule fastqc_R1:
    input:
        join(snakemake_dir, "fastq/{sample}_R1.fastq.gz")
    output:
        join(snakemake_dir, "fastqc_out/{sample}_R1_fastqc.zip"),
        join(snakemake_dir, "fastqc_out/{sample}_R1_fastqc.html")
    params:
        out_dir = join(snakemake_dir, 'fastqc_out/')
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        fastqc {input} -t {threads} -o {params.out_dir}
        """    

rule fastqc_R2:
    input:
        join(snakemake_dir, "fastq/{sample}_R2.fastq.gz")
    output:
        join(snakemake_dir, "fastqc_out/{sample}_R2_fastqc.zip"),
        join(snakemake_dir, "fastqc_out/{sample}_R2_fastqc.html")
    params:
        out_dir = join(snakemake_dir, 'fastqc_out/')
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        fastqc {input} -t {threads} -o {params.out_dir}
        """    

rule fastqc_trimmed_R1:
    input:
        temp(join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"))
    output:
        join(snakemake_dir, "trimmed_fastqc_out/{sample}_trimmed_R1_fastqc.zip"),
        join(snakemake_dir, "trimmed_fastqc_out/{sample}_trimmed_R1_fastqc.html")
    params:
        out_dir = join(snakemake_dir, 'trimmed_fastqc_out/')
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        fastqc {input} -t {threads} -o {params.out_dir}
        """    

rule fastqc_trimmed_R2:
    input:
        temp(join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz"))
    output:
        join(snakemake_dir, "trimmed_fastqc_out/{sample}_trimmed_R2_fastqc.zip"),
        join(snakemake_dir, "trimmed_fastqc_out/{sample}_trimmed_R2_fastqc.html")
    params:
        out_dir = join(snakemake_dir, 'trimmed_fastqc_out/')
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        fastqc {input} -t {threads} -o {params.out_dir}
        """    

rule star_map:
    input:
        R1 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"),
        R2 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz")
    output:
        bam_file = temp(join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam")),
        log_out = join(snakemake_dir, "bam/{sample}_Log.out"),
        final_out = join(snakemake_dir, "bam/{sample}_Log.final.out"),
        progress_out = temp(join(snakemake_dir, "bam/{sample}_Log.progress.out")),
        SJ_out = temp(join(snakemake_dir, "bam/{sample}_SJ.out.tab")),
        genome_folder = temp(directory(join(snakemake_dir, "bam/{sample}__STARgenome"))),
        pass_folder = temp(directory(join(snakemake_dir, "bam/{sample}__STARpass1")))
    threads:
        8
    resources:
        mem = scheduler_memory(41),
        time = scheduler_time_hours(8)
    message:
        "mapping with star"
    params:
        genome_index = config['star_index'],
        output_prefix = join(snakemake_dir, "bam/{sample}_"),
    shell:
        """
        STAR --genomeDir {params.genome_index} --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --runMode alignReads --outSAMattributes All --outSAMstrandField intronMotif --limitBAMsortRAM 40000000000 --outSAMtype BAM SortedByCoordinate --runThreadN {threads} --twopassMode Basic --outFileNamePrefix {params.output_prefix} 
        """

rule sam_view_strand_fix:
    input:
        join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam")
    output:
        join(snakemake_dir, "bam/{sample}.bam")
    params:
        awk_script = join(snakemake_dir, 'scripts/tagXSstrandedData.awk')
    resources:
        mem = scheduler_memory(5),
        time = scheduler_time_hours(2)
    shell:
        """
        samtools view -h {input} | awk -v strType=1 -f {params.awk_script} | samtools view -hb > {output}
        """

rule bam_index:
    input:
        join(snakemake_dir, "bam/{sample}.bam")
    output:
        join(snakemake_dir, "bam/{sample}.bam.bai")
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)
    message:
        "indexing bam file"
    shell:
        """    
        samtools index -@ 7 {input} 
        """

rule stats:
    input:
        bam = join(snakemake_dir, "bam/{sample}.bam"),
        bai = join(snakemake_dir, "bam/{sample}.bam.bai")
    output:
        join(snakemake_dir, "bam_stats/{sample}_stats.txt")
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)  
    message:
        "getting bam stats for bam"
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output}
        """

rule kallisto_quant:
    input:
        R1 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"),
        R2 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz")
    output:
        abundance_h5 = join(snakemake_dir, "kallisto_out/{sample}/abundance.h5"),
        abundance_tsv = join(snakemake_dir, "kallisto_out/{sample}/abundance.tsv"),
        run_info = join(snakemake_dir, "kallisto_out/{sample}/run_info.json")
    params:
        kallisto_index = join(snakemake_dir, "references/kallisto.index"),
        out_dir = join(snakemake_dir, "kallisto_out/{sample}")
    threads:
        1
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(6)
    message:
        "quantifying with kallisto"
    shell:
        """
        kallisto quant --fr-stranded --bias -b 1 -i {params.kallisto_index} -t {threads} -o {params.out_dir} {input.R1} {input.R2}
        """

rule junction_count:
    input:
        bam = join(snakemake_dir, "bam/{sample}.bam"),
        bai = join(snakemake_dir, "bam/{sample}.bam.bai")
    output:
        join(snakemake_dir, "junction_counts/{sample}_junctions.out")
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(6)
    message:
        "Counting junctions"
    shell:
        """
        regtools junctions extract -a 8 -m 50 -M 500000 {input.bam} -o {output} -s XS
        """    

rule rnaseqc:
    input:
        collapsed_gtf = join(snakemake_dir, "references/gencode.vM35.basic.annotation_rnaseqc_collapsed.gtf"),
        bam = join(snakemake_dir, "bam/{sample}.bam"),
        bai = join(snakemake_dir, "bam/{sample}.bam.bai")
    output:
        temp(join(snakemake_dir, "bam_stats/{sample}.exon_reads.gct")),
        temp(join(snakemake_dir, "bam_stats/{sample}.gene_fragments.gct")),
        temp(join(snakemake_dir, "bam_stats/{sample}.gene_reads.gct")),
        temp(join(snakemake_dir, "bam_stats/{sample}.gene_tpm.gct")),
        join(snakemake_dir, "bam_stats/{sample}.metrics.tsv")
    params:
        out_dir = join(snakemake_dir, "bam_stats/"),
        sample_id = "{sample}"
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(2)
    message:
        "running rnaseqc samples"
    shell:
        """
        rnaseqc -s {params.sample_id} {input.collapsed_gtf} {input.bam} {params.out_dir} 
        """    

