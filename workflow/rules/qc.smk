rule fastqc:
    input:
        rawread="data/{sample}_{read}.fastq.gz"
    output:
        zip ="pre-analysis/{sample}/fastqc/{sample}_{read}_fastqc.zip",
        html="pre-analysis/{sample}/fastqc/{sample}_{read}_fastqc.html"
    threads: 1
    
    conda:
        "ngsmo"

    params:
        path="pre-analysis/{sample}/fastqc/"
    shell:
        "fastqc {input.rawread} --threads {threads} -o {params.path}"



rule alignment_qc:   
    output:
        trackDb = "pre-analysis/alignment_stats.html",
    params:
        countTables= "pre-analysis/count-tables"

    script:
        "../scripts/alignment_stats.R"