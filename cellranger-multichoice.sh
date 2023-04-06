 #!/bin/bash
 
 #samples=$(ls SRR/ | grep SRR)
 samples=$(ls -d */)

 samp=$(sed 's/ /-/g' metadata_Bobadilla_2023.csv| tail -n +2)
 spATAC='spatial-ATAC'
 snATAC='single-nuclei-ATAC'
 visium='visium'
 jpgs=$(ls -f *jpg)
 
 ref_ATAC='/home/tbou/data/ref_genome/refdata-cellranger-arc-mm10-2020-A-2.0.0'
 ref_visium='/home/tbou/data/ref_genome/refdata-gex-mm10-2020-A'
 
 for folds in $samp;
    do 
        STR=$(echo $folds)
        fastqs=$(echo $folds | awk -F ',' '{print $5}'| awk -F '_L' '{print $1}')
        sample=$(echo $fastqs | awk -F '_S' '{print $1}')
        
        AGE=$(echo $folds | awk -F ',' '{print $2}'| sed 's/\./-/g')
        rep=$(echo $folds | awk -F ',' '{print $4}')
        #echo $sample
        #echo $sample
        #RUN SPATIAL ATAC
        
    if [[ "$STR" =~ .*"$spATAC".* ]] && [[ ! -d ${sample}_atacome ]]; then
        
        echo  "prout \n"
      
        
        cellranger-atac count --id=${sample}_atacome --reference=${ref_ATAC} --fastqs=${fastqs} --localcores=40 --localmem=500 --sample=${sample}
        
        #RUN sn ATAC
   elif [[ "$STR" =~ .*"$snATAC".* ]] && [[ ! -d ${sample}_atacome ]]; then
       cellranger-atac count --id=${sample}_atacome --reference=${ref_ATAC} --fastqs=${fastqs} --localcores=40 --localmem=500 --sample=${sample}
       
       
    elif [[ "$STR" =~ .*"$visium".* ]] && [[ ! -d ${sample}_visiome ]] ; then
            
            
            slide=$(echo ${fastqs} | awk -F '_' '{print $1}' )
            capArea=$(echo ${fastqs} | awk -F '_' '{print $2}' )
            
            e12="220421"
            
                if [[ "$slide" = "$e12" ]];
                then 
                
                    slide="V10B01-031"
                    capArea="B1"
                             
                fi
            
            image=$(echo $jpgs | tr ' ' '\n'| grep $slide )
            
        
        spaceranger count --id=${sample}_visiome --transcriptome=${ref_visium} --fastqs=${fastqs} --sample=${sample} --image=${image} --slide=${slide} --area=${capArea} --localcores=40 --localmem=500
    

       
    fi
    
    

        
    done
 