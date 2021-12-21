SAMPLE_LIST = config["sample"]
GENOME_VERSION = config["genome"]
OUT = config["output"]
IN = config["input"]

import random

## parameter initialization
if GENOME_VERSION == "hg38":
    GENOME_INDEX = "/conglilab/shared/genomics/pipeline_atac/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    GENOME_BLACKLIST = "/conglilab/shared/genomics/pipeline_atac/hg38/hg38.blacklist.withJDB.sorted.bed"
    CHROM_FILTER_1 = "v"
    CHROM_FILTER_2 = "V"
    TSS_BED = "/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/main/hg38.v32.pc.lv12.tss.bed.gz"
    TSS_EXTEND_2000 = "/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/main/hg38.v32.pc.lv12.tss.ext2k.bed"
    CHROMSIZES = "/conglilab/shared/genomics/pipeline_atac/hg38/hg38.chrom.sizes"
    MACS2_GENOME_SIZE = "hs"
elif GENOME_VERSION == "mm10":
    GENOME_INDEX = "/conglilab/shared/genomics/pipeline_atac/mm10/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.fasta"
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
        expand(OUT + "/qc/summary.csv"),
        expand(OUT + "/bam/{sample}.nodup.clean.bam", sample=SAMPLE_LIST),
        expand(OUT + "/bigWig/binned/{sample}.no_extend.rpkm.bw", sample=SAMPLE_LIST),
        expand(OUT + "/qc/libComplexity/{sample}.pbc_qc.csv", sample=SAMPLE_LIST),
        expand(OUT + "/qc/fastqc/{sample}", sample=SAMPLE_LIST),


checkpoint detectAdapter:
    input:
        fq1=IN + "/{sample}.fastq.gz",
    output:
        adapter1=OUT + "/qc/adapter/{sample}.adapter.log",
    conda:
        "envs/conda.yaml"
    shell:
        """
        python3 \
        /conglilab/shared/pipelines/atacseq_pipelines/atac_dnase_pipelines/utils/detect_adapter.py \
        {input.fq1} \
        > {output.adapter1}
        """


rule cutAdapter:
    input:
        fq1=IN + "/{sample}.fastq.gz",
        adp1=rules.detectAdapter.output.adapter1,
    output:
        fq1=OUT + "/fastq/{sample}.trimmed.fastq.gz",
        report=OUT + "/qc/cutadapt/{sample}.cutadapt_report.txt",
    threads: 6
    conda:
        "envs/conda.yaml"
    shell:
        """
        adaptor_seq1=$(cat {input.adp1} |sed -n 9p |cut -f 3 )

        ## 
        # cutadapt
        ## make threads half of specified value, because of pigz will eat as many cpu resources
        cutadapt \
        -j `echo "scale=0;{threads}/2"|bc` -m 20 -e 0.1 -O 3 \
        -q 20 --quality-base=33 \
        -a ${{adaptor_seq1}} \
        -o {output.fq1} \
        {input.fq1} \
        > {output.report}
        """


rule fastqc:
    input:
        fq1=IN + "/{sample}.fastq.gz",
    output:
        out=directory(OUT + "/qc/fastqc/{sample}"),
    threads: 1  #do not increase this
    conda:
        "envs/conda.yaml"
    shell:
        """
        mkdir -p {output.out}
        fastqc --quiet \
        -t {threads} \
        {input.fq1} \
        -o {output.out} \
        """


def cut_if_adapter(wildcards):
    r1_adap = checkpoints.detectAdapter.get(sample=wildcards.sample).output[0]
    r1_adap = (
        subprocess.check_output(f"cat {r1_adap} | sed -n 9p | cut -f 3", shell=True)
        .decode("utf-8")
        .strip()
    )
    if r1_adap == "":
        return {
            "fq1": IN + "/{sample}.fastq.gz",
        }
    else:
        return {
            "fq1": OUT + "/fastq/{sample}.trimmed.fastq.gz",
        }


rule bowtie2:
    input:
        unpack(cut_if_adapter),
    output:
        bam=OUT + "/bam/intermediate/{sample}.bam",
        bai=OUT + "/bam/intermediate/{sample}.bam.bai",
    params:
        GENOME_INDEX=GENOME_INDEX,
    threads: 16
    log:
        bowtie2=OUT + "/logs/bowtie2/{sample}.log",
    conda:
        "envs/conda.yaml"
    shell:
        """
        bowtie2 \
        --mm \
        -t -q -N1 -L 25 --no-mixed --no-discordant \
        --threads {threads} \
        -x {params.GENOME_INDEX} \
        -U {input.fq1} \
        2>{log.bowtie2} |\
        samtools view -@ {threads} -Su /dev/stdin |\
        samtools sort -@ {threads} -m 2G - > {output.bam}

        samtools index -@ {threads} {output.bam}
        """


rule filter_1:
    input:
        bam=rules.bowtie2.output.bam,
    output:
        bam=temp(OUT + "/bam/intermediate/{sample}.filter_1.bam"),
    threads: 8
    conda:
        "envs/conda.yaml"
    shell:
        """
        ## filter 1/3    
        # (use MAPQ instead of processing uniquely mapped reads;  uniquely mapping rarely mentioned today )
        # flag: filter 1804=1024+512+256+8+4 ; get 2
        # MAPQ > 30
        # sort by name 
        # 2820 stands for read unmapped, not primary alignment, failed platform/vendor quality checks, read is PCR or optical duplicate
        samtools view -F 2820 -q 30 -@ {threads} -u {input.bam} |\
        samtools sort -@ {threads} -m 2G -n /dev/stdin -o {output.bam}
        """


rule filter_2:
    input:
        bam=rules.filter_1.output.bam,
    output:
        bam=temp(OUT + "/bam/intermediate/{sample}.filter_2.bam"),
    threads: 8
    conda:
        "envs/conda.yaml"
    shell:
        """
        ## filter 2/3
        # fix mate info of name sorted bam
        # sort by coordinate again 
        # flag is the sample as above
        samtools sort -@ {threads} -m 2G {input.bam} -o {output.bam}
        """


rule markdup:
    input:
        bam=rules.filter_2.output.bam,
    output:
        marked=temp(OUT + "/bam/intermediate/{sample}.marked.bam"),
        nodup=temp(OUT + "/bam/{sample}.nodup.bam"),
        bai=OUT + "/bam/{sample}.nodup.bam.bai",
        metrics=OUT + "/qc/markdup/{sample}.dup.qc",
    params:
        tempdir=".tmp/markdup/{sample}/" + str(random.randint(1000000, 9999999)),
    threads: 50
    conda:
        "envs/conda.yaml"
    shell:
        """
        ## filter 3/3
        # picard mark duplicates (not remove) 
        # use samtools view -F 1024(in 1804) to filter, better than picard ?

        # had better specify a java temp path for Markduplicates 
        #     or it might cause error when the default system path is full
        # use the lastest version picard

        mkdir -p {params.tempdir}

        _JAVA_OPTIONS="-Xmx10G -Djava.io.tmpdir={params.tempdir}" \
        picard MarkDuplicates \
        INPUT={input.bam} \
        OUTPUT={output.marked} \
        METRICS_FILE={output.metrics} \
        VALIDATION_STRINGENCY=LENIENT \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=false

        samtools view -F 2820 -@ 6 -b -u {output.marked} |\
        samtools sort -@ 6 -m 2G /dev/stdin -o {output.nodup} 

        samtools index -@ 6 {output.nodup}
        rm -rvf {params.tempdir}
        """


# this is where nodup.clean.bam gnerated
rule mask_blacklist_region:
    input:
        bam=rules.markdup.output.nodup,
        bai=rules.markdup.output.bai,
    output:
        bam=OUT + "/bam/{sample}.nodup.clean.bam",
        bai=OUT + "/bam/{sample}.nodup.clean.bam.bai",
    params:
        GENOME_BLACKLIST=GENOME_BLACKLIST,
    threads: 8
    conda:
        "envs/conda.yaml"
    shell:
        """
        ## add one more step to filter GENOME_BLACKLIST
        # 2820 for single end
        bedtools intersect -v -abam {input.bam} -b {params.GENOME_BLACKLIST} |\
        samtools view -F 2820 -@ {threads} -S -h -b |\
        samtools sort -@ {threads} -m 2G  /dev/stdin -o  {output.bam}

        samtools index -@ {threads} {output.bam}
        """


rule flagstat:
    input:
        nodup=rules.markdup.output.nodup,
        marked=rules.markdup.output.marked,
    output:
        nodupMetric=OUT + "/qc/flagstat/{sample}.nodup.flagstat",
        markedMetric=OUT + "/qc/flagstat/{sample}.marked.flagstat",
    threads: 8
    conda:
        "envs/conda.yaml"
    shell:
        """
        samtools flagstat -@ {threads} {input.nodup} > {output.nodupMetric}
        samtools flagstat -@ {threads} {input.marked} > {output.markedMetric}
        """


rule libComplexity:
    input:
        bam=rules.markdup.output.marked,
    output:
        qc=OUT + "/qc/libComplexity/{sample}.pbc_qc.csv",
    threads: 1
    conda:
        "envs/conda.yaml"
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


rule binnedBigWig:
    input:
        bam=rules.mask_blacklist_region.output.bam,
        bai=rules.mask_blacklist_region.output.bai,
    output:
        bigWig=OUT + "/bigWig/binned/{sample}.no_extend.rpkm.bw",
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


rule atacSummary:
    input:
        rawbam=rules.bowtie2.output.bam,
        nodup=rules.markdup.output.nodup,
        bowtieLog=rules.bowtie2.log.bowtie2,
    output:
        summary=OUT + "/qc/summary/{sample}/{sample}.csv",
    threads: 1
    params:
        TSS_EXTEND_2000=TSS_EXTEND_2000,
    conda:
        "envs/conda.yaml"
    shell:
        """
          echo -e "GENOME_INDEX,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag,Frag_mt,Frag_mt%,Frag_noMT,Alignment_rate,dup%,genomic_dup%" \
          > {output.summary}
          echo "write header ok"

          Reads=$(samtools view -c {input.rawbam} )
          Reads_m=$(samtools view -c {input.rawbam} chrM)
          Reads_m_r=$(echo "scale=4;${{Reads_m}}/${{Reads}}" |bc)
          Frag=$(samtools view -c {input.nodup})
          Frag_m=$(samtools view -c {input.nodup} chrM)
          Frag_m_r=$(echo "scale=4;${{Frag_m}}/${{Frag}}" |bc)
          Frag_n=$(echo "${{Frag}}-${{Frag_m}}" |bc)
          Align=$(cat {input.bowtieLog} |grep "alignment rate" |gawk '{{print $1}}' )
          Align_d=$(echo 0.01*${{Align}} |cut -d "%" -f 1 |bc)
          Reads_n=$(printf "%.0f\\n" `echo "${{Reads}}*${{Align_d}}-${{Reads_m}}" |bc`)
          dup_r=$(echo "scale=4;1-${{Frag}}/(${{Reads}}*${{Align_d}})" |bc )
          genomic_dup_r=$(echo "scale=4;(${{Reads_n}}-${{Frag_n}})/${{Reads_n}}" |bc )
        echo "other stat ok"
          echo -e {wildcards.sample}","${{Reads}}","${{Reads_m}}","${{Reads_m_r}}","${{Reads_n}}","${{Frag}}","${{Frag_m}}","${{Frag_m_r}}","${{Frag_n}}","${{Align}}","${{dup_r}}","${{genomic_dup_r}} \
          >> {output.summary}
        echo "final echo ok"
        """


rule summaryConcat:
    input:
        csvList=expand(OUT + "/qc/summary/{sample}/{sample}.csv", sample=SAMPLE_LIST),
    output:
        csv=OUT + "/qc/summary.csv",
    threads: 1
    run:
        headerLine = "GENOME_INDEX,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag,Frag_mt,Frag_mt%,Frag_noMT,Alignment_rate,dup%,genomic_dup%\n"
        with open(output.csv, "w") as out:
            out.write(headerLine)
        with open(output.csv, "a") as out:
            for summary in input.csvList:
                with open(summary, "r") as __in:
                    summaryLine = (__in.readlines())[1].strip()
                    out.write(summaryLine + "\n")
