ruleorder: STAR_PE > STAR_SE
rule STAR_SE:
    input:
        "data/{sample}_R1.fastq.gz"

    output:
        aligned_bam="pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"

    log:
        detailed_log="pre-analysis/{sample}/STAR/Log.out",
        mapping_summary="pre-analysis/{sample}/STAR/Log.final.out"

    
    threads: 12

    params:
        star = "--outSAMtype BAM SortedByCoordinate \
			  --outFilterType BySJout --outSJfilterCountUniqueMin -1 2 2 2 --outSJfilterCountTotalMin -1 2 2 2 \
			  --outFilterIntronMotifs RemoveNoncanonical --twopassMode Basic \
			  --outSAMattributes NH NM MD --outSAMunmapped Within --outSAMmapqUnique 50 --limitBAMsortRAM 24000000000 \
              --readFilesCommand zcat --quantMode GeneCounts",
        index = config["ref"]["star"],
        outdir = "pre-analysis/{sample}/STAR/",
        gtf = config["ref"]["star_gtf"]
    conda:
        "ngsmo"

    shell:
        """
        STAR --genomeDir {params.index} --sjdbGTFfile {params.gtf} --readFilesIn {input} \
         --runThreadN {threads} {params.star} --outFileNamePrefix {params.outdir} 
    
        """

rule STAR_PE:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"

    output:
        aligned_bam="pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"

    log:
        detailed_log="pre-analysis/{sample}/STAR/Log.out",
        mapping_summary="pre-analysis/{sample}/STAR/Log.final.out"

    
    threads: 12

    params:
        star = "--outSAMtype BAM SortedByCoordinate \
			  --outFilterType BySJout --outSJfilterCountUniqueMin -1 2 2 2 --outSJfilterCountTotalMin -1 2 2 2 \
			  --outFilterIntronMotifs RemoveNoncanonical --twopassMode Basic \
			  --outSAMattributes NH NM MD --outSAMunmapped Within --outSAMmapqUnique 50 --limitBAMsortRAM 24000000000 \
              --readFilesCommand zcat --quantMode GeneCounts",
        index = config["ref"]["star"],
        outdir = "pre-analysis/{sample}/STAR/",
        gtf = config["ref"]["star_gtf"]
    conda:
        "ngsmo"

    shell:
        """
        STAR --genomeDir {params.index} --sjdbGTFfile {params.gtf} --readFilesIn {input.r1} {input.r2} \
         --runThreadN {threads} {params.star} --outFileNamePrefix {params.outdir} 
    
        """