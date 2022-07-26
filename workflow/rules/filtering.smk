rule markDup:
    input:
        bamFile= "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"     
    output:
        bamMarked= "pre-analysis/{sample}/STAR/aligned_markDup.bam",
        # bamIndex= "pre-analysis/{sample}/STAR/aligned_markDup.bam.bai",
        dupStats="pre-analysis/{sample}/STAR/dupMetrics.tsv" 
    params:
        samtools="-ASO=coordinate --TAGGING_POLICY All --CREATE_INDEX true"
    threads: 4

    conda:
        "ngsmo"
    log:
        "pre-analysis/{sample}/logs/MarkDuplicates.log"

    shell:
        "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"




rule index:
    input:
        "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.markdup.out.bam"        
    output:
        index="pre-analysis/{sample}/STAR/Aligned.sortedByCoord.markdup.out.bam.bai"  
    params:
        index="index --show-progress "
    threads: 4

    conda:
        "sambamba"

    shell:
        """
        sambamba {params.index} -t {threads} {input} {output.index}
        
        """

rule samtools:
    input:
        "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.markdup.out.bam" 
        
    output:
        "pre-analysis/{sample}/STAR/unfiltered_aligned_stats.tsv"
    threads: 4

    conda:
        "ngsmo"

    shell:
        """
        samtools flagstat -O tsv -@ {threads} {input} > {output}
        
        """

rule mark_dups:
    input:
        "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"       
    output:
        "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.markdup.out.bam" 
    params:
        "markdup --show-progress "
    threads: 4

    conda:
        "sambamba"

    shell:
        """
        sambamba {params} -t {threads} {input} {output}
        
        """