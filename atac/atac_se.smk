sampleList=config["sample"]
genomeVersion=config["genome"]
OUT=config["output"]
IN=config["input"]
shiftExtendSize=config.get("shiftExtendSize", [125, 250])


import random

## parameter initialization
if genomeVersion=="hg38":
  Index='/conglilab/shared/genomics/pipeline_atac/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta'
  blacklist='/conglilab/shared/genomics/pipeline_atac/hg38/hg38.blacklist.withJDB.sorted.bed'
  filter1='v'
  filter2='V'
  TSS_BED='/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/main/hg38.v32.pc.lv12.tss.bed.gz'
  TSS_extend='/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/main/hg38.v32.pc.lv12.tss.ext2k.bed'
  CHROMSIZES="/conglilab/shared/genomics/pipeline_atac/hg38/hg38.chrom.sizes"
  macs2_genome_size="hs"
elif genomeVersion=="mm10":
  Index='/conglilab/shared/genomics/pipeline_atac/mm10/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.fasta'
  blacklist='/conglilab/shared/genomics/pipeline_atac/mm10/mm10.blacklist.withJDB.sorted.bed'
  filter1='random'
  filter2='chrUn'
  TSS_BED='/conglilab/shared/genomics/mouse/GRCm38_2019/Annotations/GRCm38_gencode.m23/main/mm10.vM23.pc.lv12.tss.bed.gz'
  TSS_extend='/conglilab/shared/genomics/mouse/GRCm38_2019/Annotations/GRCm38_gencode.m23/main/mm10.vM23.pc.lv12.tss.ext2k.bed'
  CHROMSIZES="/conglilab/shared/genomics/pipeline_atac/mm10/mm10.chrom.sizes"
  macs2_genome_size="mm"
else:
  print("Invalid Genome Version.")

rule all:
  input:
    expand(OUT+"/qc/summary.csv"),
    expand(OUT+"/bam/{sample}.nodup.clean.bam",sample=sampleList),
    expand(OUT+"/bigWig/binned/{sample}.no_extend.rpkm.bw",sample=sampleList),

    expand(OUT+"/qc/libComplexity/{sample}.pbc_qc.csv",sample=sampleList),
    expand(OUT+"/qc/tssplot/{sample}",sample=sampleList),
    expand(OUT+"/qc/tssplot/{sample}/{sample}_tss-enrich.txt",sample=sampleList),

    expand(OUT+"/qc/fastqc/{sample}",sample=sampleList),

    expand(OUT+"/qc/tss/{sample}_reads_in_tss.txt",sample=sampleList),
    expand(OUT+"/qc/tss/{sample}_reads_catched_tss.txt",sample=sampleList),

    expand(OUT+"/bed/{sample}.shift.bed",sample=sampleList),
    expand(OUT+"/bed/{sample}.insert.bed",sample=sampleList),
    expand(OUT+"/bed/ext{ext_n}/{sample}.ext{ext_n}.bed", sample=sampleList, ext_n=shiftExtendSize),
    expand(OUT+"/bigWig/unbinned/ext{ext_n}/{sample}.ext{ext_n}.unbinned.cpm.bw",sample=sampleList, ext_n=shiftExtendSize),

    # expand(OUT+"/peak/ext{ext_n}", ext_n=shiftExtendSize),
    expand(OUT+"/peak/ext{ext_n}/{sample}_ext{ext_n}_peaks.narrowPeak", sample=sampleList, ext_n=shiftExtendSize),
    expand(OUT+"/bam/ext{ext_n}/{sample}.ext{ext_n}.bam", sample=sampleList, ext_n=shiftExtendSize),
    expand(OUT+"/bam/ext{ext_n}/{sample}.ext{ext_n}.bam.bai", sample=sampleList, ext_n=shiftExtendSize),

rule detectAdapter:
  input:
    fq1=IN+"/{sample}.fastq.gz",
  output:
    adapter1=OUT+"/qc/adapter/{sample}.adapter.log",
  resources:
    memPerThread= "100m"
  shell:"""
    python3 \
    /conglilab/shared/pipelines/atacseq_pipelines/atac_dnase_pipelines/utils/detect_adapter.py \
    {input.fq1} \
    > {output.adapter1}
"""

rule cutAdapter:
  input:
    fq1=IN+"/{sample}.fastq.gz",
    adp1=rules.detectAdapter.output.adapter1,
  output:
    fq1=OUT+"/fastq/{sample}.trimmed.fastq.gz",
    report=OUT+"/qc/cutadapt/{sample}.cutadapt_report.txt"
  threads: 6
  resources:
    memPerThread= "100m"
  shell:"""
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
    fq1=IN+"/{sample}.fastq.gz",
  output:
    out=directory(OUT+"/qc/fastqc/{sample}")
  threads: 2 #do not increase this
  resources:
    memPerThread= "500m"
  shell:"""
  mkdir -p {output.out}
  fastqc --quiet \
  -t {threads} \
  {input.fq1} \
  -o {output.out} \
  -d {output.out} 
  """

rule bowtie2:
  input:
    fq1=rules.cutAdapter.output.fq1,
  output:
    bam=OUT+"/bam/intermediate/{sample}.bam",
    bai=OUT+"/bam/intermediate/{sample}.bam.bai"
  params:
    Index=Index,
  threads: 16 
  resources:
    memPerThread= "2G"
  log:
    bowtie2=OUT+"/logs/bowtie2/{sample}.log"
  shell:"""
  bowtie2 \
  -X2000 \
  --mm \
  -t -q -N1 -L 25 --no-mixed --no-discordant \
  --threads {threads} \
  -x {params.Index} \
  -U {input.fq1} \
  2>{log.bowtie2} |\
  samtools view -@ {threads} -Su /dev/stdin |\
  samtools sort -@ {threads} -m 2G - > {output.bam}
  
  samtools index -@ {threads} {output.bam}
  """ 

rule filter1:
  input:
    bam=rules.bowtie2.output.bam
  output:
    bam=temp(OUT+"/bam/intermediate/{sample}.filter1.bam")
  threads: 8
  resources:
    memPerThread= "2G"
  shell:"""
  ## filter 1/3    
  # MAPQ > 30
  # sort by name 
  # 2820 stands for read unmapped, not primary alignment, failed platform/vendor quality checks, read is PCR or optical duplicate
  samtools view -F 2820 -q 30 -@ {threads} -u {input.bam} |\
  samtools sort -@ {threads} -m 2G -n /dev/stdin -o {output.bam}
  """

rule filter2:
  input:
    bam=rules.filter1.output.bam
  output:
    bam=temp(OUT+"/bam/intermediate/{sample}.filter2.bam")
  threads:8
  resources:
    memPerThread= "2G"
  shell:"""
  samtools sort -@ {threads} -m 2G {input.bam} -o {output.bam}
  """

rule markdup:
  input:
    bam=rules.filter2.output.bam
  output:
    marked=temp(OUT+"/bam/intermediate/{sample}.marked.bam"),
    nodup=temp(OUT+"/bam/{sample}.nodup.bam"),
    bai=OUT+"/bam/{sample}.nodup.bam.bai",
    metrics=OUT+"/qc/markdup/{sample}.dup.qc"
  params:
    tempdir=".tmp/markdup/{sample}/"+str(random.randint(1000000,9999999))
  threads: 50
  resources:
    memPerThread= "2G"
  shell:"""
  ## filter 3/3
  # picard mark duplicates (not remove) 
  # use samtools view -F 1024(in 1804) to filter, better than picard ?
  
  # had better specify a java temp path for Markduplicates 
  #     or it might cause error when the default system path is full
  # use the lastest version picard
  
  mkdir -p {params.tempdir}
  
  java -Xmx10G -Djava.io.tmpdir={params.tempdir} \
  -jar /conglilab/shared/applications/picard_2.20.4/picard.jar MarkDuplicates \
  INPUT={input.bam} \
  OUTPUT={output.marked} \
  METRICS_FILE={output.metrics} \
  VALIDATION_STRINGENCY=LENIENT \
  ASSUME_SORTED=true \
  REMOVE_DUPLICATES=false
  
  samtools view -F 2820 -@ 5 -b -u {output.marked} |\
  samtools sort -@ 5 -m 2G /dev/stdin -o {output.nodup} 
  
  samtools index -@ 5 {output.nodup}
  rm -rvf {params.tempdir}
  """

rule blacklist:
  input:
    bam=rules.markdup.output.nodup,
    bai=rules.markdup.output.bai,
  output:
    bam=OUT+"/bam/{sample}.nodup.clean.bam",
    bai=OUT+"/bam/{sample}.nodup.clean.bam.bai",
  params:
    blacklist=blacklist
  threads: 8
  resources:
    memPerThread= "2G"
  shell:"""
  ## add one more step to filter blacklist
  bedtools intersect -v -abam {input.bam} -b {params.blacklist} |\
  samtools view -F 2820 -@ {threads} -S -h -b |\
  samtools sort -@ {threads} -m 2G  /dev/stdin -o  {output.bam}
  
  samtools index -@ {threads} {output.bam}
  """

rule flagstat:
  input:
    nodup=rules.markdup.output.nodup,
    marked=rules.markdup.output.marked,
  output:
    nodupMetric=OUT+"/qc/flagstat/{sample}.nodup.flagstat",
    markedMetric=OUT+"/qc/flagstat/{sample}.marked.flagstat",
  threads: 8
  resources:
    memPerThread= "100m"
  shell:"""
  samtools flagstat -@ {threads} {input.nodup} > {output.nodupMetric}
  samtools flagstat -@ {threads} {input.marked} > {output.markedMetric}
  """

rule libComplexity:
  input:
    bam=rules.markdup.output.marked
  output:
    qc=OUT+"/qc/libComplexity/{sample}.pbc_qc.csv"
  threads: 1
  resources:
    memPerThread= "10m"
  shell:"""
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
    bam=rules.blacklist.output.bam,
    bai=rules.blacklist.output.bai,
  output:
    bigWig=OUT+"/bigWig/binned/{sample}.no_extend.rpkm.bw"
  threads: 8
  resources:
    memPerThread= "1G"
  shell:"""
  bamCoverage \
  -p {threads} \
  -bs 100 \
  --normalizeUsing RPKM \
  -b {input.bam} \
  -o {output.bigWig}
  """

rule bedGeneration:
  input:
    bam=rules.blacklist.output.bam,
  output:
    shift=temp(OUT+"/bed/{sample}.shift.bed"),
    insert__=temp(OUT+"/bed/{sample}.insert.bed"),
  threads: 1
  resources:
    memPerThread= "1G"
  params:
    f1=filter1,
    f2=filter2,
    blacklist=blacklist
  shell:"""
  bedtools intersect -v -abam {input.bam} -b {params.blacklist} |\
  bedtools bamtobed -i - |gawk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' |\
  grep -v chrM |grep -v {params.f1} |grep -v {params.f2} |\
  gawk -F "\\t" 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2=$2+4}} else if($6=="-"){{$3=$3-5}} print $0}}' \
  > {output.shift}
  cat {output.shift} |gawk '{{ if($6=="+"){{print $1"\\t"$2-1"\\t"$2"\\t"$4;}} else{{print $1"\\t"$3-1"\\t"$3"\\t"$4}} }}' \
  > {output.insert__}
  """

rule bedExtension:
  input:
    insert_bed=rules.bedGeneration.output.insert__
  output:
    ext_bed=OUT+"/bed/ext{ext_n}/{sample}.ext{ext_n}.bed"
  threads:1
  shell: """ 
  cat {input.insert_bed} |gawk '{{print $1"\\t"$3-{wildcards.ext_n}"\\t"$3+{wildcards.ext_n}"\\t"$4}}' \
  > {output.ext_bed}
  """

rule extened_bed_to_bam:
  input:
    bed=rules.bedExtension.output.ext_bed,
  output:
    bam=OUT+"/bam/ext{ext_n}/{sample}.ext{ext_n}.bam",
    bai=OUT+"/bam/ext{ext_n}/{sample}.ext{ext_n}.bam.bai",
  params:
    genome_size=CHROMSIZES
  threads:1
  shell: """ 
  bedToBam -i {input.bed} -g {params.genome_size}|samtools sort /dev/stdin > {output.bam} 
  samtools index {output.bam}
  """

rule unbinned_shift_extended_bigwig:
  input:
    bed=rules.bedExtension.output.ext_bed
  output:
    bw=OUT+"/bigWig/unbinned/ext{ext_n}/{sample}.ext{ext_n}.unbinned.cpm.bw"
  params:
    genome_size=CHROMSIZES
  shell:"""
  bed={input.bed}
  genome_size={params.genome_size}
  bigWig={output.bw}
  
  scalingFactor=$(echo "scale=10;1000000/`cat $bed|wc -l`"|bc)
  bedtools genomecov \
    -i $bed \
    -g $genome_size\
    -scale $scalingFactor \
    -bga |\
  sort --parallel 1 -k1,1 -k2,2n |\
  grep -v Un|grep -v v|grep -v chrM|grep -v chrEBV > ${{bed}}.bedGraph

  bedGraphToBigWig ${{bed}}.bedGraph $genome_size $bigWig
  rm ${{bed}}.bedGraph
  """

rule peak_calling:
  input:
    bed=rules.bedExtension.output.ext_bed,
  output:
    narrow_peak=OUT+"/peak/ext{ext_n}/{sample}_ext{ext_n}_peaks.narrowPeak",
  threads: 1
  params:
    peak_dir=OUT+"/peak/ext{ext_n}",
    genome_size=macs2_genome_size
  shell:"""
    macs2 callpeak \
    -t {input.bed} \
    -f BED \
    -g {params.genome_size} \
    --nomodel \
    --bdg \
    --call-summits \
    --outdir {params.peak_dir} \
    --name {wildcards.sample}_ext{wildcards.ext_n}
    """

rule tss:
  input:
    bam=rules.blacklist.output.bam
  output:
    plotdir=directory(OUT+"/qc/tssplot/{sample}"),
    enrichTSS=OUT+"/qc/tssplot/{sample}/{sample}_tss-enrich.txt",
    qc1=OUT+"/qc/tss/{sample}_reads_in_tss.txt",
    qc2=OUT+"/qc/tss/{sample}_reads_catched_tss.txt",
  threads: 1
  resources:
    memPerThread= "1G"
  params:
    chromsize=CHROMSIZES,
    tssbed=TSS_BED,
    TSS_extend=TSS_extend
  shell:"""
  ## tss
  # load conda env built for kundaje pipeline
  #     dependency in this python script hasn't been independent yet
  # alias conda=/ds918_208/shared/applications/conda/miniconda3/condabin/conda
  set +u; source /ds918_208/shared/applications/conda/miniconda3/etc/profile.d/conda.sh; set -u
  set +u;conda activate bds_atac;set -u
  
  
  tss_plot='/conglilab/shared/projects/sc_omics/scripts/utils/tss_plot.py'
  
  python ${{tss_plot}} \
  --outdir {output.plotdir} \
  --outprefix {wildcards.sample} \
  --tss {params.tssbed} \
  --finalbam {input.bam} \
  --chromsizes {params.chromsize}
  
  intersectBed -a {params.TSS_extend} -b {input.bam} |wc -l > {output.qc1}
  intersectBed -a {params.TSS_extend} -b {input.bam} -wa |sort -u |wc -l > {output.qc2}
  """

rule atacSummary:
  input:
    rawbam=rules.bowtie2.output.bam,
    nodup=rules.markdup.output.nodup,
    bowtieLog=rules.bowtie2.log.bowtie2,
    enrichTSS=rules.tss.output.enrichTSS,
    inTSS=rules.tss.output.qc1,
    catchTSS=rules.tss.output.qc2,

  output:
    summary=OUT+"/qc/summary/{sample}/{sample}.csv"
  threads: 1
  resources:
    memPerThread= "10m"
  params:
    TSS_extend=TSS_extend
  shell:"""
  # alias conda=/ds918_208/shared/applications/conda/miniconda3/condabin/conda
  set +u; source /ds918_208/shared/applications/conda/miniconda3/etc/profile.d/conda.sh; set -u
  echo "conda source ok"
  set +u; conda activate base; set -u
  echo "conda activate ok"

  echo -e "Index,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag,Frag_mt,Frag_mt%,Frag_noMT,Alignment_rate,dup%,genomic_dup%,Tss_fold_enrichment,TSS_updown_2kb_reads,inTSS_ratio,TSS_updown_2kb_catched,catchTSS_ratio" \
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
  TSS=$(printf $(cat {input.enrichTSS} )"\\n" )
  inTSS=$(cat {input.inTSS} )
  inTSS_r=$(echo "scale=4;${{inTSS}}/${{Frag_n}}" |bc)
  catchTSS=$(cat {input.catchTSS} )
  TSS_n=$(cat {params.TSS_extend} |wc -l)
  catchTSS_r=$(echo "scale=4;${{catchTSS}}/${{TSS_n}}" |bc)
echo "other stat ok"
  echo -e {wildcards.sample}","${{Reads}}","${{Reads_m}}","${{Reads_m_r}}","${{Reads_n}}","${{Frag}}","${{Frag_m}}","${{Frag_m_r}}","${{Frag_n}}","${{Align}}","${{dup_r}}","${{genomic_dup_r}}","${{TSS}}","${{inTSS}}","${{inTSS_r}}","${{catchTSS}}","${{catchTSS_r}} \
  >> {output.summary}
echo "final echo ok"
  """

rule summaryConcat:
  input:
    csvList=expand(OUT+"/qc/summary/{sample}/{sample}.csv",sample=sampleList)
  output:
    csv=OUT+"/qc/summary.csv"
  threads: 1
  resources:
    memPerThread= "10m"
  run:
    headerLine="Index,Reads,Reads_mt,Reads_mt%,Reads_noMT,Frag,Frag_mt,Frag_mt%,Frag_noMT,Alignment_rate,dup%,genomic_dup%,Tss_fold_enrichment,TSS_updown_2kb_reads,inTSS_ratio,TSS_updown_2kb_catched,catchTSS_ratio\n"
    with open(output.csv,"w") as out:
      out.write(headerLine)
    with open(output.csv,"a") as out:
      for summary in input.csvList:
        with open(summary,"r") as __in:
          summaryLine=(__in.readlines())[1].strip()
          out.write(summaryLine+"\n")
