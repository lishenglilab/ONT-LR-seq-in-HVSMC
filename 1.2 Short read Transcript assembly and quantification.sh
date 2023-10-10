## Transcript assembly and quantification from RNA-seq data.

#!/bin/bash
# trimmomatic version=0.39
# HISAT2 version=2.7.6a
# Stringtie version=2.1.8
# genome version GRCh38
# gtf version gencode.v38.annotation.gtf


sample=$1

## quality control
dir_raw='/data/raw'
dir_clean='data_clean'
trim='/softwares/trimmomatic-0.39.jar'
mkidr -p ${dir_clean}
fq1=${dir_raw}/$i/$i'_1.fq.gz'
fq2=${dir_raw}/$i/$i'_2.fq.gz'
filt_fn_r1=${dir_clean}/${sample}/${sample}_filtered_R1.fastq.gz
unp_fn_r1=${dir_clean}/${sample}/${sample}_unpaired_R1.fastq.gz
filt_fn_r2=${dir_clean}/${sample}/${sample}_filtered_R2.fastq.gz
unp_fn_r2=${dir_clean}/${sample}/${sample}_unpaired_R2.fastq.gz
java -jar ${trim} PE ${fq1} ${fq2} ${dir_trim}/$i/${filt_fn_r1} ${dir_trim}/$i/${unp_fn_r1} \
    ${dir_trim}/$i/${filt_fn_r2} ${dir_trim}/$i/${unp_fn_r2} \
    'ILLUMINACLIP:'${adapter}':2:30:10' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    -threads 12 2> ${dir_clean}/${sample}/read_surviving_stat.txt

##  read alignment

hisat='/softwares/HISAT2'
gdir='/ref/hisat_gencodev38_index'
file_gtf='/ref/gencode.v38.annotation.gtf'
dir_alignment='/data/star_alignments'
mkdir -p ${dir_alignment}/${sample}
fq1=${dir_clean}/${sample}/${sample}_filtered_R1.fastq.gz
fq2=${dir_clean}/${sample}/${sample}_filtered_R2.fastq.gz
 
${hisat} -x ${gdir} -1 ${fq1} -2 ${fq2} -S ${dir_alignment}/${sample}.sam -p 12
samtools view -bS ${dir_alignment}/${sample}.sam --threads 6 > ${dir_alignment}/${sample}.bam
samtools sort ${dir_alignment}/${sample}.bam -o ${dir_alignment}/${sample}.sort.bam --threads 6
samtools index${dir_alignment}/${sample}.sort.bam -@ 6
samtools flagstat ${dir_alignment}/${sample}.sort.bam --threads 6 > $i/flagstat.txt




## Stringtie reference-based transcript assembly

stringtie='/sotfwares/stringtie'
file_bam=${dir_align}/${sample}/${sample}.sort.bam
dir_stringtie='/data/stringtie_assembly'
mkdir -p ${dir_stringtie}'/'${sample}
out_gtf=${dir_stringtie}/${sample}/${sample}_stringtie.gtf
${stringtie} ${file_bam} -p 5 -o ${out_gtf} -G ${file_gtf}

## merge transcript assembly from individual samples

gtf_files='/data/gtf_files.txt' # list of transcript assembly of all samples
${stringtie} --merge -p 10 -o ${out_gtf} -G ${ref_gtf} ${gtf_files}

## Stringtie transcript quantification
dir_stringtie='/mnt/data132Tp3/public/201212/stringtie_quantification'
mkdir -p ${dir_stringtie}/${sample}
ref_gtf='/data/merged.gtf'
out_gtf=${dir_stringtie}/${sample}/${sample}'_stringtie.gtf'
file_gene=${dir_stringtie}/${sample}/${sample}'_gene_abundance.tab'
${stringtie} ${file_bam} -p 8 -o ${out_gtf} -G ${ref_gtf} -A ${file_gene} -e