## Find ORF

#!/bin/bash
# ORFfinder version=0.4.3
tr_fa='/Data/Flair/flair.collapse.isoforms.fa'
orf_fa='/Data/ORFfind/transcript_orf.fa'
/sotfwares/ORFfinder -in ${tr_fa} -strand plus -g 1 -ml 100 -out ${orf_fa} -outfmt 1