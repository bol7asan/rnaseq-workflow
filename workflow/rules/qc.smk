rule fastqc:
    input:
        rawread="data/{sample}_{read}.fastq.gz"
    output:
        zip="pre-analysis/{sample}/fastqc/{sample}_{read}_fastqc.zip",
        html="pre-analysis/{sample}/fastqc/{sample}_{read}_fastqc.html"
    threads: 1
    
    container:
        "docker://biocontainers/fastqc:v0.11.9_cv8"  # Apptainer will pull this Docker image
    conda:
        "ngs-env"
    shell:
        "fastqc {input.rawread} --threads {threads} -o pre-analysis/{wildcards.sample}/fastqc/"



rule alignment_qc:   
    output:
        trackDb = "pre-analysis/alignment_stats.html",
    params:
        countTables= "pre-analysis/count-tables"

    script:
        "../scripts/alignment_stats.R"