## Transcript assembly and quantification from RNA-seq data.

#!/bin/bash
# porechop version=0.2.4
# FLAIR version=1.4
# genome version GRCh38
# gtf version gencode.v38.annotation.gtf

sample=$1

## remove adapters


porechop='/softwares/porechop'
dir_raw='/data/raw'
dir_clean='/data/clean'
# gdir='/mnt/data132Tp3/public/201212/ref/star_gencodev35_grch38.p13_index'
# file_gtf='/mnt/data132Tp3/public/201212/ref/gencode.v35.annotation.gtf'
mkdir -p ${dir_clean}/${sample}
fq=${dir_fq}/${sample}/pass.fq.gz
${porechop} -i ${fq} -o ${dir_clean}/${sample}/out.pass.fq.gz


## FLAIR for alignment and transcript assembly

flair='/softwares/flair'
dir_flair='/data/flair'
clean_fq=${dir_clean}/${sample}/pass.fq.gz
mkdir -p ${dir_flair}/${sample}
ref_gtf='/ref/gencode.v38.annotation.gtf'
ref_fa='/ref/GRCh38.primary_assembly.genome.fa'

### FLAIR align
${flair} align --reads ${clean_fq} --genome ${ref_fa} --output ${dir_flair}/${sample}/${sample} --threads 10

### FLAIR correct
query_bed=${dir_flair}/${sample}/${sample}.bed
flair correct --query ${query_bed} --genome ${ref_fa} --gtf ${gtf} --output ${dir_flair}/${sample}/${sample} --threads 10

## FLAIR collapse
all_bed=${dir_flair}/all_samples_corrected.bed
find ${dir_flair} -name *_all_corrected.bed | xargs cat > ${all_bed}
bam_list=`find ${dir_flair} -name *.bam`
reads=`find ${dir_raw} -name pass.fq.gz`
flair collapse --temp_dir ${dir_flair} -g ${ref_fa} -q ${all_bed} -r ${reads} -o ${dir_flair}/merged_collapse --stringent --threads ${threads_num} --gtf ${gtf}

## FLAIR quantify
flair quantify --temp_dir ${dir_flair} -r ${dir_flair}/reads_manifest.txt -i ${dir_flair}/merged_collapse.isoforms.fa --threads 10 --tpm --output ${dir_flair}/transcript
