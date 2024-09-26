# Define common variables for STAR parameters, log paths, and threads
common_params = {
    "star": "--outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout --outSJfilterCountUniqueMin -1 2 2 2 \
            --outSJfilterCountTotalMin -1 2 2 2 --outFilterIntronMotifs RemoveNoncanonical \
            --twopassMode Basic --outSAMattributes NH NM MD --outSAMunmapped Within \
            --outSAMmapqUnique 50 --limitBAMsortRAM 16000000000 --readFilesCommand 'gzip -dc' \
            --quantMode GeneCounts",
    "index": config["ref"]["star"],
    "gtf": config["ref"]["star_gtf"]
}

common_log = {
    "detailed_log": "pre-analysis/{sample}/STAR/Log.out",
    "mapping_summary": "pre-analysis/{sample}/STAR/Log.final.out"
}

threads = 4

# STAR rule for Single-End Reads
rule STAR_SE:
    input:
        "data/{sample}_R1.fastq.gz"

    output:
        aligned_bam="pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"

    log:
        **common_log

    threads: threads

    params:
        **common_params,
        outdir="pre-analysis/{sample}/STAR/"

    conda:
        "star_env"

    shell:
        """
        STAR --genomeDir {params.index} --sjdbGTFfile {params.gtf} --readFilesIn {input} \
         --runThreadN {threads} {params.star} --outFileNamePrefix {params.outdir}
        """

# STAR rule for Paired-End Reads
rule STAR_PE:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"

    output:
        aligned_bam="pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"

    log:
        **common_log

    threads: threads

    params:
        **common_params,
        outdir="pre-analysis/{sample}/STAR/"

    conda:
        "star_env"

    shell:
        """
        STAR --genomeDir {params.index} --sjdbGTFfile {params.gtf} --readFilesIn {input.r1} {input.r2} \
         --runThreadN {threads} {params.star} --outFileNamePrefix {params.outdir}
        """

# Ensure rule order for paired-end over single-end
ruleorder: STAR_PE > STAR_SE