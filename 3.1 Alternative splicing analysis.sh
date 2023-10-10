## Alternative splicing analysis

#!/bin/bash
# SUPPA2 version=2.3

suppa='/softwares/suppa.py'
ref_gtf='/data/Filair/filair.collapse.isoforms.gtf'

## generateEvents
events='/data/SUPPA/events'
${suppa} generateEvents -i ${ref_gtf} -o ${events} -e SE SS MX RI FL -f ioe

## merge all events
mer_events='/data/SUPPA/merged.all.events.ioe'
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > ${mer_events}

## calculate PSI
exp='/data/Filair/counts_matrix.tsv'
psi_mat='data/SUPPA/project_event'
${suppa} -i ${mer_events} -e ${exp} -o ${psi_mat}
