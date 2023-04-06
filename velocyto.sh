#!/bin/bash

RUN="YES"

date


SAMPLE=$( ls -d * | grep "AD"|sed 's!.*/!!' )


for fold in $( echo $SAMPLE );
do   
    printf "\n"
    echo "Folder##################"
    echo $fold
    
    if [ "$RUN" = "YES" ]
        then

            printf "SAMPLE = $fold \n"
         

        ~/miniconda3/envs/RNAvelo/bin/velocyto run10x -m ~/data/ref_genome/GRCh38_rmsk.gtf ${fold}/${fold}_transcriptome ~/data/ref_genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf
    
    fi

done
