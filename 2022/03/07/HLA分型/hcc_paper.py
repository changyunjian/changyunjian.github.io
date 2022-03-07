'''
Author: your name
Date: 2022-02-18 18:57:03
LastEditTime: 2022-03-06 16:43:44
LastEditors: Please set LastEditors
Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
FilePath: /kingdom/code/2022_neoantigen/hcc_paper.py
'''



sample1_list = ['tumor']
sample2_list = ['tumor','normal']

rule all:
    input:
        expand("STAR_result/{sample1}Aligned.out.sam",sample1=sample1_list),
        expand("mpileup_result/{sample2}.mpileup",sample2=sample2_list)


rule bwa_mapping:
    input:
        'raw_data/{sample}_R1.fastq',
        'raw_data/{sample}_R2.fastq'
    params:
        '/sibcb1/wuliganglab1/changyunjian/annotation/human/bwa_index/index'
    output:
        'bwa_result/{sample}.bam'
    log:
        'logs/quality_filter/{sample}.log'
    shell:
        "/sibcb/program/bin/bwa mem -t 32 -R '@RG\tID:foo_lane\tSM:sample_name\tLB:length' "
        "{params} {input[0]} {input[1]} > {output[0]} 2>{log} "


rule RNA_mapping:
    input:
        "raw_data/{sample1}_RNA_R1.fastq",
        "raw_data/{sample1}_RNA_R2.fastq"
    params:
        index='/sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/STAR_index',
        dir_out='STAR_result/{sample1}'
    output:
        'STAR_result/{sample1}Aligned.out.sam'
    shell:
        '~/software/STAR-2.7.5a/bin/Linux_x86_64/STAR --runMode alignReads --runThreadN 32 --genomeDir {params.index} '
        '--readFilesIn {input} --outFileNamePrefix {params.dir_out} --outFilterMultimapNmax 10000 '
        '--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.06 '
        '--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --twopassMode Basic '
        '--outSAMprimaryFlag AllBestScore --outReadsUnmapped Fastx'

rule sort_bam:
    input:
        'wes_bam/{sample2}-wes.bam'
    output:
        'wes_bam/{sample2}_sort.bam'
    shell:
        'samtools sort -o {output} {input}'

rule mpileup:
    input:
        'wes_bam/{sample2}_sort.bam'
    output:
        'mpileup_result/{sample2}.mpileup'
    params:
        '/sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa'
    shell:
        '/sibcb/program/bin/samtools mpileup -q 1 -f {params} {input} > {output} '

###########################################################################
# gatk

'''
下载dbsnp ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
下载clinvar ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
下载gnomAD http://hgdownload.cse.ucsc.edu/gbdb/hg19/gnomAD/vcf/
下载ExAC http://hgdownload.cse.ucsc.edu/gbdb/hg19/ExAC/
'''
#index dict
'samtools dict -o  red.dict ref.fasta'

# conda enviroment
"conda activate gatk"

# deduplicate
"java -jar /sibcb1/wuliganglab1/changyunjian/software/picard/picard.jar MarkDuplicates \
    INPUT=normal_sort.bam OUTPUT=normal_mark_dup.bam METRICS_FILE=normal_dup.txt"
#碱基质量重矫正（BQSR）
# normal and tumor
"gatk --java-options '-XX:ParallelGCThreads=20 -Xmx30G -Djava.io.tmpdir=./tmp' BaseRecalibrator \
    -I normal_mark_dup_fix.bam \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    --known-sites /sibcb1/wuliganglab1/changyunjian/annotation/human/gatk_snp_data/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    --known-sites /sibcb1/wuliganglab1/changyunjian/annotation/human/gatk_snp_data/hg19/dbsnp_138.hg19.vcf \
    --known-sites /sibcb1/wuliganglab1/changyunjian/annotation/human/gatk_snp_data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -O normal_recal_data.table1"

"gatk  --java-options '-XX:ParallelGCThreads=20 -Xmx30G -Djava.io.tmpdir=./tmp' ApplyBQSR \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    -I normal_mark_dup_fix.bam \
    --bqsr-recal-file  normal_recal_data.table1 \
    -O normal_mark_dup_fix_bqsr.bam"

#mutect2 call varient

"gatk'-XX:ParallelGCThreads=20 -Xmx30G -Djava.io.tmpdir=./tmp' Mutect2 \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    -I wes_bam/normal_mark_dup_fix_bqsr.bam \
    -max-mnp-distance 0 \
    -O ./GATK_result/normal.vcf.gz"

#
"gatk --java-options '-XX:ParallelGCThreads=20 -Xmx30G -Djava.io.tmpdir=./tmp' GenomicsDBImport  \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    -L /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19_interval.list \
    --genomicsdb-workspace-path GATK_result/pon_db \
    -V GATK_result/normal.vcf.gz "

#
"gatk --java-options '-XX:ParallelGCThreads=20 -Xmx30G -Djava.io.tmpdir=./tmp' CreateSomaticPanelOfNormals \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    -V gendb://GATK_result/pon_db \
    - O GATK_result/pon.vcf.gz "
#

"bcftools annotate -x ^INFO/AF gnomad.exomes.r2.0.2.sites.vcf.gz > af-only-gnomad.hg19.vcf"
#mutect2 call varient
"gatk --java-options '-XX:ParallelGCThreads=30 -Xmx30G -Djava.io.tmpdir=./tmp' Mutect2 \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    -I wes_bam/normal_mark_dup_fix_bqsr.bam \
    -I wes_bam/tumor_mark_dup_fix_bqsr.bam \
    -normal WHYO_AC \
    --germline-resource /sibcb1/wuliganglab1/changyunjian/annotation/human/gatk_snp_data/Mutect2/af-only-gnomad.hg19_chr.vcf.gz \
    --panel-of-normals ./GATK_result/pon.vcf.gz \
    --f1r2-tar-gz ./GATK_result/f1r2.tar.gz \
    -O ./GATK_result/somatic.vcf.gz "

# filter varient

"gatk --java-options '-XX:ParallelGCThreads=30 -Xmx30G -Djava.io.tmpdir=./tmp' LearnReadOrientationModel \
    -I ./GATK_result/f1r2.tar.gz \
    -O ./GATK_result/read-orientation-model.tar.gz"

"gatk --java-options '-XX:ParallelGCThreads=30 -Xmx30G -Djava.io.tmpdir=./tmp' GetPileupSummaries \
    -I ./wes_bam/tumor_mark_dup_fix_bqsr_chr1.bam \
    -V /sibcb1/wuliganglab1/changyunjian/annotation/human/gatk_snp_data/Mutect2/GetPileupSummaries/small_exac_common_3_hg19.vcf.gz \
    -L /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19_interval.list \
    -O ./GATK_result/getpileupsummaries.table"


"gatk --java-options '-XX:ParallelGCThreads=30 -Xmx30G -Djava.io.tmpdir=./tmp' CalculateContamination \
    -I ./GATK_result/getpileupsummaries.table \
    -tumor-segmentation ./GATK_result/segments.table \
    -O ./GATK_result/calculatecontamination.table"


"gatk --java-options '-XX:ParallelGCThreads=30 -Xmx30G -Djava.io.tmpdir=./tmp' FilterMutectCalls \
    -R /sibcb1/wuliganglab1/changyunjian/annotation/human/UCSC_hg19/hg19.fa \
    -V ./GATK_result/somatic.vcf.gz \
    --min-median-mapping-quality 20 \
    -tumor-segmentation ./GATK_result/segments.table \
    --contamination-table ./GATK_result/calculatecontamination.table \
    --ob-priors ./GATK_result/read-orientation-model.tar.gz \
    -O ./GATK_result/somatic_filtered.vcf"

#
"bcftools merge -m none ./strelka2/results/variants/somatic.snvs.vcf.gz \
    ./GATK_result/somatic_filtered.vcf.gz -Ov -o somatic_result/combine_snv_v.vcf""

##############################################################
# varscan
"java -jar VarScan.jar somatic {normal.pileup} {tumor.pileup} output "

##############################################################
# somaticsniper
"bam-somaticsniper -q 40 -G -L -f {refer.fa} {tumor.bam} {normal.bam} -F vcf {output} "

##############################################################
#strelka
'''
/sibcb1/wuliganglab1/changyunjian/software/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam {normal.bam} \
    --tumorBam {tumor.bam} \
    --referenceFasta {refer.fa} \
    --runDir strelka2 --exome

/sibcb2/wuliganglab2/changyunjian/2022_neoantigen_hla/strelka2/runWorkflow.py -m local -j 40

'''


##############################################################
# optitype
"OptiTypePipeline.py -i {r1.fastq} {r2.fastq} --rna -v -o ./optitype "

"OptiTypePipeline.py -i {r1.fastq} {r2.fastq} --dna -v -o ./optitype "

###############################################################
# seq2HLA
#run under the dict /seq2HLA
"conda activate py2.7"
"seq2HLA -1 raw_data/tumor_RNA_R1.fastq.gz -2 raw_data/tumor_RNA_R2.fastq.gz -r seq2hla -p 40"

###############################################################
#PHLAT
#ssh snode046 -q docker.q

"docker commit -m 'phlat with vim' -a 'changyunjian' 8926734a5599 changyunjian/phlat"

"docker run -P --name rna \
    -v /sibcb2/wuliganglab2/changyunjian/2022_neoantigen_hla/raw_data:/opt/phlat-release/example \
    -v /sibcb2/wuliganglab2/changyunjian/2022_neoantigen_hla/PHLAT:/opt/phlat-release/results \
        changyunjian/phlat bash /opt/phlat-release/example.sh"

"docker cp rna:/opt/phlat-release/results/tumor_RNA_HLA.sum ./PHLAT/"



# pVACseq
#conda env gatk

"vcf-genotype-annotator strelka2/results/variants/somatic.snvs.vcf.gz test 0/1 -o pVACseq/genotype.vcf"




