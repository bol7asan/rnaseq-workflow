# rule htseq:
#     input:
#         bam = "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.markdup.out.bam",
#         gtf = config["ref"]["star_gtf"]      
#     output:
#         forward = "pre-analysis/{sample}/HTSEQ/gene_count.txt",
#         rev = "pre-analysis/{sample}/HTSEQ/gene_countR.txt"
#     params:
#         rev = "-s reverse -r pos  -m union -i gene_id --additional-attr chr --additional-attr start --additional-attr end --additional-attr gene_name --additional-attr gene_type",
#         forward = "-s yes -r pos  -m union -i gene_id --additional-attr chr --additional-attr start --additional-attr end --additional-attr gene_name --additional-attr gene_type" 

#     threads: 4

#     log:
#         "pre-analysis/{sample}/logs/HtSeq.log"

#     conda:
#         "ngsmo"

#     shell:
#         """
#         htseq-count {params.forward} {input.bam} {input.gtf} > {output.forward}
#         htseq-count {params.rev} {input.bam} {input.gtf} > {output.rev}
        
#         """

rule featureCounts:
    input:
        bam = "pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam",
        gtf = config["ref"]["star_gtf"]      
    output:
        "pre-analysis/{sample}/featureCounts/gene_count.txt"
    params:
        "-s 2 --extraAttributes  chr,gene_name,gene_type -p -M -O"

    threads: 4

    log:
        "pre-analysis/{sample}/logs/featureCounts.log"

    conda:
        "ngsmo"

    shell:
        """
        featureCounts -T {threads} -a {input.gtf} {params} -o {output} {input.bam} &> {log}
        
        """
rule make_raw_table:
    output:
        "pre-analysis/count-tables/raw_counts.tab"
    params:
        countDir = "pre-analysis"

    threads: 4

    conda:
        "pyngs"

    script:
        "../scripts/assemble_featureCounts.py"