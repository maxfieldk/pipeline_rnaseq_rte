rule common_sortIndexBam:
    input:
        sam = "outs/{sample}/star_output/Aligned.out.sam"
    output:
        sortedbam =  "outs/{sample}/star_output/{sample}.sorted.bam",
        stats =  "outs/{sample}/star_output/{sample}.sorted.bam.stats.txt",
        index = "outs/{sample}/star_output/{sample}.sorted.bam.bai"
    resources:
        cpus_per_task =10,
        mem_mb = 64000
    conda:
        "omics"
    shell:
        """
samtools view -@8 -b {input.sam} | \
samtools sort -@8 -m4g - > {output.sortedbam}
samtools index  -@6 {output.sortedbam}
samtools stats {output.sortedbam} > {output.stats}
        """


rule common_filterForPrimaryAlignments:
    input:
        "{path}.sorted.bam"
    output:
        bam =  "{path}.sorted.primary.bam",
        bamindex = "{path}.sorted.primary.bam.bai"
    threads: 4
    conda:
        "omics"
    shell: 
        """
samtools view -b -F 0x800 -F 0x100 -F 0x400 {input} > {output.bam}
samtools index {output.bam}
        """ 