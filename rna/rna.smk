sampleList=config["sample"]
genomeVersion=config["genome"]
OUT=config["output"]
IN=config["input"]


import random


if genomeVersion=="hg38":
  RRNA_BOWTIE_INDEX="/ds918_208/shared/genomics/FOR_PIPELINE/HG38/RRNA/human.ribo.rna"
  STAR_INDEX='/conglilab/shared/genomics/human/GRCh38_2019'
  RSEM_TRANS='/conglilab/shared/genomics/human/GRCh38_2019/Annotations/GRCh38_gencode.v32/rsem_trans_index'

elif genomeVersion=="mm10":
  RRNA_BOWTIE_INDEX="/ds918_208/shared/genomics/FOR_PIPELINE/MM10/RRNA/mouse.ribo.rna"
  STAR_INDEX='/conglilab/shared/genomics/mouse/GRCm38_2019'
  RSEM_TRANS='/conglilab/shared/genomics/mouse/GRCm38_2019/Annotations/GRCm38_gencode.m23/rsem_trans_index'
else:
  print("Invalid Genome Version.")


rule all:
  input:
    dataSummary=expand("{OUT}/qc/summary.csv",OUT=OUT),
    rawCount=expand(OUT+"/matrix/counts.matrix"),
    tpm=expand(OUT+"/matrix/tpm.matrix"),
    fpkm=expand(OUT+"/matrix/fpkm.matrix"),

    RRNA=expand(OUT+"/qc/RRNA/{sample}/rRNA.sam.frag_count",sample=sampleList),
    transcriptomeBam=expand(OUT+"/bam/star/{sample}.Aligned.toTranscriptome.out.bam",sample=sampleList),
    genomeBam=expand(OUT+"/bam/star/{sample}.Aligned.sortedByCoord.out.bam",sample=sampleList),
    genomeBai=expand(OUT+"/bam/star/{sample}.Aligned.sortedByCoord.out.bam.bai",sample=sampleList),
    genomeBigWig=expand(OUT+"/bigWig/{sample}.bw",sample=sampleList),
    
    fastqc=expand("{OUT}/qc/fastqc/{sample}",OUT=OUT,sample=sampleList),



checkpoint detectAdapter:
  input:
    fq1=IN + "/{sample}_R1.fastq.gz",
    fq2=IN + "/{sample}_R2.fastq.gz",
  output:
    adapter1=OUT + "/qc/adapter/{sample}_R1.adapter.log",
    adapter2=OUT + "/qc/adapter/{sample}_R2.adapter.log",
  conda:
    "envs/conda.yaml"
  shell:
    """
    python3 \
    /conglilab/shared/pipelines/atacseq_pipelines/atac_dnase_pipelines/utils/detect_adapter.py \
    {input.fq1} \
    > {output.adapter1}
    python3 \
    /conglilab/shared/pipelines/atacseq_pipelines/atac_dnase_pipelines/utils/detect_adapter.py \
    {input.fq2} \
    > {output.adapter2}
    """

def cut_if_adapter(wildcards):
  r1_adap = checkpoints.detectAdapter.get(sample=wildcards.sample).output[0]
  r2_adap = checkpoints.detectAdapter.get(sample=wildcards.sample).output[1]
  r1_adap = (
    subprocess.check_output(f"cat {r1_adap} | sed -n 9p | cut -f 3", shell=True)
    .decode("utf-8")
    .strip()
  )
  r2_adap = (
    subprocess.check_output(f"cat {r2_adap} | sed -n 9p | cut -f 3", shell=True)
    .decode("utf-8")
    .strip()
  )
  if any([r1_adap == "", r2_adap == ""]):
    return {
      "fq1": IN + "/{sample}_R1.fastq.gz",
      "fq2": IN + "/{sample}_R2.fastq.gz",
    }
  else:
    return {
      "fq1": OUT + "/fastq/{sample}_R1.trimmed.fastq.gz",
      "fq2": OUT + "/fastq/{sample}_R2.trimmed.fastq.gz",
    }

rule cutAdapter:
  input:
    fq1=IN+"/{sample}_R1.fastq.gz",
    fq2=IN+"/{sample}_R2.fastq.gz",
    adp1=rules.detectAdapter.output.adapter1,
    adp2=rules.detectAdapter.output.adapter2,
  output:
    fq1=OUT+"/fastq/{sample}_R1.trimmed.fastq.gz",
    fq2=OUT+"/fastq/{sample}_R2.trimmed.fastq.gz",
    report=OUT+"/qc/cutadapt/{sample}.cutadapt_report.txt"
  threads: 6
  shell:"""
    adaptor_seq1=$(cat {input.adp1} |sed -n 9p |cut -f 3 )
    adaptor_seq2=$(cat {input.adp2} |sed -n 9p |cut -f 3 )
    
    ## 
    # cutadapt
    ## make thread half of specified, becuase of pigz will eat as many cpu resources 
    cutadapt \
    -j `echo "scale=0;{threads}/2"|bc` -m 20 -e 0.1 -O 3 \
    -q 20 --quality-base=33 \
    -a ${{adaptor_seq1}} \
    -A ${{adaptor_seq2}} \
    -o {output.fq1} \
    -p {output.fq2} \
    {input.fq1} \
    {input.fq2} \
    > {output.report}
    """

rule fastqc:
  input:
    fq1=IN+"/{sample}_R1.fastq.gz",
    fq2=IN+"/{sample}_R2.fastq.gz",
  output:
    out=directory(OUT+"/qc/fastqc/{sample}")
  threads: 2 #do not increase this
  shell:"""
  mkdir -p {output.out}
  fastqc \
  -t {threads} \
  {input.fq1} \
  {input.fq2} \
  -o {output.out} \
  -d {output.out} 
  """


## note for rule rrna
# params.qc specifies the output dir of the perl script
# output.countResult is the output FILE of the perl script, but cannot be set directly from the shell
# so i devided it 
rule rrna:
  input:
    unpack(cut_if_adapter)
  output:
    countResult=OUT+"/qc/RRNA/{sample}/rRNA.sam.frag_count"
  threads: 16
  params:
    RRNA_BOWTIE_INDEX=RRNA_BOWTIE_INDEX
  run:
    import subprocess
    bowtie2Log=subprocess.run(
      " ".join([
        '/usr/local/bin/bowtie2',
        '--mm',
        '-t',
        '-q',
        '--no-mixed',
        '--no-discordant',
        '--threads '+str(threads),
        '-X2000',
        '-x '+params.RRNA_BOWTIE_INDEX,
        '-1 '+input.fq1,
        '-2 '+input.fq2,
        '1>/dev/null'
      ]),
      capture_output=True,
      shell=True,
      text=True
    ).stderr
    print(bowtie2Log)
    t1=re.search("(\d+) \((.*)\) aligned concordantly exactly 1 time",bowtie2Log)
    rrnaCount=t1.group(1)
    with open(output.countResult,'w') as f:
      f.write(rrnaCount)



## note for rule rrna
# params.prefix specifies the output dir of STAR
# output.countResult is the output FILE of STAR, if you need to change the cmd of STAR, please change the output file aswell
# so i devided it 
rule star:
  input:
    unpack(cut_if_adapter)
  output:
    coordsortBam=OUT+"/bam/star/{sample}.Aligned.sortedByCoord.out.bam",
    transcriptomeBam=OUT+"/bam/star/{sample}.Aligned.toTranscriptome.out.bam",
    logFinal=OUT+"/bam/star/{sample}.Log.final.out",
  params:
    STAR_INDEX=STAR_INDEX,
    prefix=OUT+"/bam/star/{sample}."
  threads: 16
  resources:
    mem_mb=30*1024
  shell: """
  /conglilab/shared/applications/staraln/STAR-2.5.3a/bin/Linux_x86_64_static/STAR \
  --genomeDir {params.STAR_INDEX} \
  --genomeLoad LoadAndKeep \
  --runThreadN {threads} \
  --readFilesIn \
  {input.fq1} \
  {input.fq2} \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --limitBAMsortRAM 30000000000 \
  --outFileNamePrefix {params.prefix}  \
  --outSAMstrandField intronMotif \
  --quantMode TranscriptomeSAM
  """

rule bamIndex:
  input:
    bam=rules.star.output.coordsortBam,
  output:
    bai=rules.star.output.coordsortBam+".bai",
  threads: 6
  shell: """
  samtools index -@ {threads} {input.bam}
  """


rule bamCoverage:
  input:
    bam=rules.star.output.coordsortBam,
    bai=rules.bamIndex.output.bai
  output:
    bw=OUT+"/bigWig/{sample}.bw"
  threads: 6
  shell: """
  bamCoverage \
  -p {threads} \
  -b {input.bam}  \
  -o {output.bw}
  """

rule rsem:
  input:
    bam=rules.star.output.transcriptomeBam,
  output:
    results=OUT+"/rsem/{sample}.genes.results",
    isoform=OUT+"/rsem/{sample}.isoforms.results"
  params:
    RSEM_TRANS=RSEM_TRANS,
    prefix=OUT+"/rsem/{sample}"
  threads: 8
  shell: """
  /conglilab/shared/applications/rsem/RSEM-1.3.0/rsem-calculate-expression \
  -p {threads} \
  --no-bam-output \
  --paired-end \
  --alignments {input.bam} \
  {params.RSEM_TRANS} \
  {params.prefix}
  """

rule summaryConcat:
  input:
    starLog=expand(rules.star.output.logFinal,sample=sampleList),
    geneResult=expand(rules.rsem.output.results,sample=sampleList),
    rrnaCount=expand(rules.rrna.output.countResult,sample=sampleList)
  output:
    csv=OUT+"/qc/summary.csv"
  run:
    import re
    import pandas as pd
    def writeInfo(i,fileHandle): 
      starLog=input.starLog[i]    
      rrnaCount=input.rrnaCount[i]    
      geneResultFile=input.geneResult[i]


      with open(starLog,"r") as f:
        starLogDict=[x.strip() for x in f.readlines()]
      starLogDict=dict([x for x in [re.split("\ \|\t",x) for x in starLogDict] if len(x)==2])
      Reads=int(starLogDict["Number of input reads"])
      Unique=int(starLogDict["Uniquely mapped reads number"])
      Unique_r=str(starLogDict["Uniquely mapped reads %"])
      Multiple=int(starLogDict['Number of reads mapped to multiple loci'])
      Multiple_r=str(starLogDict['% of reads mapped to multiple loci'])
      
      with open(rrnaCount,'r') as f:
        rRNA=int(f.readline().strip())
        rRNA_r=str(100*rRNA/(Unique+Multiple))+"%"

      
      geneResult=pd.read_table(geneResultFile) 
      genes_1=geneResult.query("TPM>0").shape[0]
      genes_2=geneResult.query("TPM>2").shape[0]

      sampleName=re.search("(?<=\/rsem\/).*(?=.genes.results)",geneResultFile).group(0)

      fileHandle.write(",".join([str(x) for x in[
        sampleName,
        Reads,
        Unique,
        Unique_r,
        Multiple,
        Multiple_r,
        rRNA,
        rRNA_r,
        genes_1,
        genes_2
      ]])+"\n")

      # fileHandle.write(f"{sample:wildcards.sample},{reads:Reads},{unique:Unique},{unique_r:Unique_r},{multiple:Multiple},{multiple_r:Multiple_r},{rRNA:rRNA},{rRNA_r:rRNA_r},{genes_1:genes_1},{genes_2:genes_2}")


    headerLine="Sample,TotalReads,Unique,UniqueRatio,Mutiple,MutipleRatio,rRNA,rRNARatio,genes(TPM>0),genes(TPM>2)\n"
    with open(output.csv,"w") as out:
      out.write(headerLine)
    with open(output.csv,"a") as out:
      for i in range(len(input.starLog)):
        writeInfo(i,out)

rule countMatrix:
  input:
    geneResult=expand(rules.rsem.output.results,sample=sampleList),
  output:
    rawCount=OUT+"/matrix/counts.matrix",
    tpm=OUT+"/matrix/tpm.matrix",
    fpkm=OUT+"/matrix/fpkm.matrix",
  run:
    import pandas as pd
    from functools import reduce
    import re
    
    pathList=input.geneResult
    nameList=[re.search("(?<=\/rsem\/).*(?=.genes.results)",x).group(0) for x in pathList]
    
    dfList=list(map(pd.read_table,pathList))
    
    def getCount(df):
      return(df.expected_count)
    def getTPM(df):
      return(df.TPM)
    def getFPKM(df):
      return(df.FPKM)
    
    countList=list(map(getCount,dfList))
    countList.append(dfList[0].gene_id)
    countDF=pd.concat(countList,axis=1)
    countDF=countDF.set_index("gene_id")
    countDF.columns=nameList
    countDF.to_csv(output.rawCount,sep="\t",index=True,header=True,float_format="%.2f")
    
    tpmList=list(map(getTPM,dfList))
    tpmList.append(dfList[0].gene_id)
    tpmDF=pd.concat(tpmList,axis=1)
    tpmDF=tpmDF.set_index("gene_id")
    tpmDF.columns=nameList
    tpmDF.to_csv(output.tpm,sep="\t",index=True,header=True,float_format="%.2f")

    fpkmList=list(map(getFPKM,dfList))
    fpkmList.append(dfList[0].gene_id)
    fpkmDF=pd.concat(fpkmList,axis=1)
    fpkmDF=fpkmDF.set_index("gene_id")
    fpkmDF.columns=nameList
    fpkmDF.to_csv(output.fpkm,sep="\t",index=True,header=True,float_format="%.2f")
