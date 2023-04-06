#!/bin/bash


DL=$(ls -d *AD*)

for d in $DL
    do
   fold=$(echo ${d})
   sample=$(ls $fold | awk -F '_S' '{print $1}' | uniq)

    echo $fold
    echo $sample
    
   
    cellranger count --id ${fold}_transcriptome --transcriptome ~/data/ref_genome/refdata-gex-GRCh38-2020-A  --fastqs ${fold}/  --localcores 40  --localmem 500 --include-introns true --sample=$sample
    
    mv ${fold}_transcriptome ${fold}
 done
