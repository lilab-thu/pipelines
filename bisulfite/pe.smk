SAMPLE_LIST = config["sample"]
GENOME_VERSION = config["genome"]
OUT = config["output"]
IN = config["input"]

import random
import subprocess

if GENOME_VERSION == "hg38":
    GENOME_INDEX = "/conglilab/shared/genomics/pipeline_atac/hg38/bwameth_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    GENOME_BLACKLIST = "/conglilab/shared/genomics/pipeline_atac/hg38/hg38.blacklist.withJDB.sorted.bed"
    CHROM_FILTER_1 = "v"
    CHROM_FILTER_2 = "V"
    TSS_BED = "/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/main/hg38.v32.pc.lv12.tss.bed.gz"
    TSS_EXTEND_2000 = "/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/main/hg38.v32.pc.lv12.tss.ext2k.bed"
    CHROMSIZES = "/conglilab/shared/genomics/pipeline_atac/hg38/hg38.chrom.sizes"
    MACS2_GENOME_SIZE = "hs"
elif GENOME_VERSION == "mm10":
    GENOME_INDEX = "/conglilab/shared/genomics/pipeline_atac/mm10/bwameth_index/mm10_no_alt_analysis_set_ENCODE.fasta"
    GENOME_BLACKLIST = "/conglilab/shared/genomics/pipeline_atac/mm10/mm10.blacklist.withJDB.sorted.bed"
    CHROM_FILTER_1 = "random"
    CHROM_FILTER_2 = "chrUn"
    TSS_BED = "/conglilab/shared/genomics/mouse/GRCm38_2019/Annotations/GRCm38_gencode.m23/main/mm10.vM23.pc.lv12.tss.bed.gz"
    TSS_EXTEND_2000 = "/conglilab/shared/genomics/mouse/GRCm38_2019/Annotations/GRCm38_gencode.m23/main/mm10.vM23.pc.lv12.tss.ext2k.bed"
    CHROMSIZES = "/conglilab/shared/genomics/pipeline_atac/mm10/mm10.chrom.sizes"
    MACS2_GENOME_SIZE = "mm"
else:
    print("Invalid Genome Version.")


rule all:
    input:
        # expand(OUT + "/qc/summary.csv"),
        expand(OUT + "/bam/{sample}.bam", sample=SAMPLE_LIST),
        expand(OUT + "/bam/{sample}.bam.bai", sample=SAMPLE_LIST),
        expand(OUT + "/bw/binned/genome/{sample}.rpkm.bw", sample=SAMPLE_LIST),
        expand(OUT + "/bw/unbinned/meth/{sample}.bdg", sample=SAMPLE_LIST),
        expand(OUT + "/bw/unbinned/meth/{sample}.bw", sample=SAMPLE_LIST),
        expand(OUT + "/qc/fragsize/{sample}.fragsize.pdf", sample=SAMPLE_LIST),
        expand(OUT + "/qc/fragsize/{sample}.fragsize.txt", sample=SAMPLE_LIST),
        expand(OUT + "/qc/libComplexity/{sample}.pbc_qc.csv", sample=SAMPLE_LIST),
        expand(OUT + "/qc/fastp/{sample}.json", sample=SAMPLE_LIST),
        expand(OUT + "/qc/fastp/{sample}.html", sample=SAMPLE_LIST),


rule fastp:
    input:
        fq1=IN + "/{sample}_R1.fastq.gz",
        fq2=IN + "/{sample}_R2.fastq.gz",
    output:
        fq1=OUT + "/fastq/{sample}_R1.fastq.gz",
        fq2=OUT + "/fastq/{sample}_R2.fastq.gz",
        json_report=OUT + "/qc/fastp/{sample}.json",
        html_report=OUT + "/qc/fastp/{sample}.html",
    threads: 16
    shell:
        """
        fastp \
        --in1 {input.fq1}  --in2 {input.fq2} \
        --out1 {output.fq1} --out2 {output.fq2} \
        -j {output.json_report} -h {output.html_report} \
        -w {threads}
        """


rule align:
    input:
        fq1=rules.fastp.output.fq1,
        fq2=rules.fastp.output.fq2,
    output:
        bam=OUT + "/bam/intermediate/{sample}.raw.bam",
        bai=OUT + "/bam/intermediate/{sample}.raw.bam.bai",
    threads: 32
    shell:
        """
        bwameth.py \
        --threads {threads} \
        --reference {GENOME_INDEX} \
        {input.fq1} {input.fq2} |\
        samtools sort -@ 4 -m 1G -o {output.bam}

        samtools index -@ {threads} {output.bam}
        """


rule filter_1:
    input:
        bam=rules.align.output.bam,
    output:
        bam=temp(OUT + "/bam/intermediate/{sample}.filter_1.bam"),
        bai=temp(OUT + "/bam/intermediate/{sample}.filter_1.bam.bai"),
    threads: 8
    shell:
        """
        ## filter 1/3    
        # (use MAPQ instead of processing uniquely mapped reads;  uniquely mapping rarely mentioned today )
        # flag: filter 1804=1024+512+256+8+4 ; get 2
        # MAPQ > 30
        samtools view -@ {threads} -F 1804 -f 2 -q 30 -u {input.bam} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """


rule fixmate:
    input:
        bam=rules.filter_1.output.bam,
    output:
        bam=temp(OUT + "/bam/intermediate/{sample}.fixmate.bam"),
    threads: 8
    shell:
        """
        # sort by name 
        samtools sort -@ {threads} -n {input.bam} -o {output.bam}.name_sort.bam

        samtools fixmate -m -r -@ {threads} {output.bam}.name_sort.bam {output.bam}.fixmate.bam

        samtools sort -@ {threads} -o {output.bam} {output.bam}.fixmate.bam

        rm {output.bam}.name_sort.bam {output.bam}.fixmate.bam
        """


# use samtools markdup and samtools view to remove duplicates in this pipeline
rule dedup:
    input:
        bam=rules.fixmate.output.bam,
    output:
        mark_bam=temp(OUT + "/bam/intermediate/{sample}.markdup.bam"),
        mark_bai=temp(OUT + "/bam/intermediate/{sample}.markdup.bam.bai"),
        bam=temp(OUT + "/bam/intermediate/{sample}.dedup.bam"),
        bai=temp(OUT + "/bam/intermediate/{sample}.dedup.bam.bai"),
    threads: 16
    shell:
        """
        samtools markdup -@ {threads} -r {input.bam} {output.mark_bam} 
        asmtools index -@ {threads} {output.mark_bam}

        samtools view -u -F 1804 -f 2 -@4 {output.mark_bam} |\
        samtools sort -@ {threads} -o {output.bam} /dev/stdin

        samtools index -@ {threads} {output.bam}
        """


rule mask_blacklist_region:
    input:
        bam=rules.dedup.output.bam,
        bai=rules.dedup.output.bai,
    output:
        bam=OUT + "/bam/{sample}.bam",
        bai=OUT + "/bam/{sample}.bam.bai",
    params:
        GENOME_BLACKLIST=GENOME_BLACKLIST,
    threads: 8
    shell:
        """
        ## add one more step to filter GENOME_BLACKLIST
        bedtools intersect -v -abam {input.bam} -b {params.GENOME_BLACKLIST} |\
        samtools view -F 1804 -f 2 -@ {threads} -S -h -b |\
        samtools sort -@ {threads} /dev/stdin -o {output.bam}

        samtools index -@ {threads} {output.bam}
        """


rule genome_bw:
    input:
        bam=rules.mask_blacklist_region.output.bam,
        bai=rules.mask_blacklist_region.output.bai,
    output:
        bigWig=OUT + "/bw/binned/genome/{sample}.rpkm.bw",
    threads: 8
    conda:
        "envs/conda.yaml"
    shell:
        """
        bamCoverage \
        -p {threads} \
        -bs 100 \
        --normalizeUsing RPKM \
        -b {input.bam} \
        -o {output.bigWig}
        """


# use methylDackel to get CpG methylation bedgraph and bigwig
rule meth_bw:
    input:
        bam=rules.mask_blacklist_region.output.bam,
        bai=rules.mask_blacklist_region.output.bai,
    output:
        bdg=OUT + "/bw/unbinned/meth/{sample}.bdg",
        bw=OUT + "/bw/unbinned/meth/{sample}.bw",
    params:
        GENOME_INDEX=GENOME_INDEX,
        CHROMSIZES=CHROMSIZES,
    threads: 8
    shell:
        """
        MethylDackel \
        extract \
        {params.GENOME_INDEX} \
        {input.bam} \
        > {output.bdg}

        bedGraphToBigWig \
        {output.bdg} \
        {params.CHROMSIZES} \
        {output.bw}
        """


rule flagstat:
    input:
        nodup=rules.dedup.output.bam,
        marked=rules.dedup.output.mark_bam,
    output:
        nodupMetric=OUT + "/qc/flagstat/{sample}.nodup.flagstat",
        markedMetric=OUT + "/qc/flagstat/{sample}.marked.flagstat",
    threads: 8
    shell:
        """
        samtools flagstat -@ {threads} {input.nodup} > {output.nodupMetric}
        samtools flagstat -@ {threads} {input.marked} > {output.markedMetric}
        """


rule libComplexity:
    input:
        bam=rules.dedup.output.mark_bam,
    output:
        qc=OUT + "/qc/libComplexity/{sample}.pbc_qc.csv",
    threads: 1
    shell:
        """
        echo "TotalPair,DictinctPair,OnePair,TwoPair,NRF=Distinct/Total,PBC1=OnePair/Distinct,PBC2=OnePair/TwoPair" > {output.qc} 
        bedtools bamtobed -i {input.bam} |\
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}' |\
        grep -v 'chrM' |sort |uniq -c |\
        awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} \
        {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; \
        printf "%d,%d,%d,%d,%f,%f,%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}}' >>  {output.qc}
        """


rule fragsize:
    input:
        bam=rules.mask_blacklist_region.output.bam,
    output:
        pdf=OUT + "/qc/fragsize/{sample}.fragsize.pdf",
        txt=OUT + "/qc/fragsize/{sample}.fragsize.txt",
    threads: 1
    resources:
        memPerThread="4G",
    conda:
        "envs/conda.yaml"
    shell:
        """
        _JAVA_OPTIONS="-Xmx4G" \
        picard CollectInsertSizeMetrics \
        I={input.bam} \
        O={output.txt} \
        H={output.pdf} \
        VERBOSITY=ERROR QUIET=TRUE \
        W=1000
        """
