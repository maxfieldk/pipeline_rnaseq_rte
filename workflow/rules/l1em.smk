rule alignBWAtohg38:
    input:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz"
    output:
        sam = "outs/{sample}/hg38alignments/{sample}_hg38.sam"
    threads: 8
    resources:
        mem_mb  = 60000,
        runtime = 1000
    conda:
        "L1EM"
    shell:
        """
bwa mem -T 0 /oscar/data/jsedivy/mkelsey/tools/L1EM/DEPENDENCIES/hg38.fa {input.r1} {input.r2} > {output.sam}
        """


rule L1EM:
    input:
        bam = "outs/{sample}/hg38alignments/{sample}_hg38.sorted.bam"
    output:
        outfile = "outfiles/l1em{sample}.txt"
    threads: 20
    resources:
        mem_mb  = 100000,
        runtime = 1000
    conda:
        "L1EM"
    shell:
        """
wd=$(pwd)
mkdir -p outs/{wildcards.sample}/L1EM
cd outs/{wildcards.sample}/L1EM
bash -e /oscar/data/jsedivy/mkelsey/tools/L1EM/run_L1EM.sh $wd/{input.bam} /users/mkelsey/data/tools/L1EM /oscar/data/jsedivy/mkelsey/tools/L1EM/DEPENDENCIES/hg38.fa
touch $wd/{output.outfile}
        """