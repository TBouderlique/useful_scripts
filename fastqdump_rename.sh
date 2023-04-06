 #!/bin/bash
 
 #samples=$(ls SRR/ | grep SRR)
 samples=$(ls SRR/ |grep .sra)
 N=60
 RUN="YES"
 J=3
 
 #CREATE METADATA FILE
 touch metadata_Bobadilla_2023.csv
 echo SRR,AGE,EXPERIMENT,REPLICATE,READ1,READ2,READ3,READ4 >> metadata_Bobadilla_2023.csv
 
 
 for i in $samples 
     do
                  
         echo SRR/$i

         file=$(echo $i | sed 's!.*/!!' | rev | cut -c5- | rev  | sed 's!.*/!!')
         echo $file
         base_search=$(esearch -db sra -query $file | efetch -format native  -mode text) 
           
         
                    
        
#####WRITE METADATA
        AGE=$(echo $base_search  |  sed 's/.*sample_title="\(.*"\)/\1/'| cut -d'"' -f1 | tail -1 | sed 's/, /,/g' | cut -d',' -f2)
        replicate=$(echo $base_search | sed 's/.*sample_title="\(.*"\)/\1/'| cut -d'"' -f1| tail -1 | sed 's/, /,/g' | cut -d',' -f3)
        EXPERIMENT=$(echo $base_search | sed 's/.*sample_title="\(.*"\)/\1/'| cut -d'"' -f1| tail -1 | sed 's/, /,/g' | awk -F, '{print $NF}')

        echo replicate is equal to $replicate
        
        
        
        #get reads for Visium
         if [ "$EXPERIMENT" ==  "visium" ];
        then
                    
        
         read1=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep R1|head -1| sed 's/ //g')
         read2=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep R2| head -1| sed 's/ //g')
         read3=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep I1| head -1| sed 's/ //g')
         read4=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep I2|head -1|  sed 's/ //g')
         
         else
         read1=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep R1|head -1|  sed 's/ //g')
         read2=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep R2| head -1| sed 's/ //g')
         read3=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep R3| head -1| sed 's/ //g')
         read4=$(echo $base_search  |tr '><' '\n'| grep -E "read|filename"|tr '\"' '\n'| grep fastq.gz | grep -v read | grep I1| head -1| sed 's/ //g')
        fi
        
        
         #Create folder for fastqs
         folder=$(echo $read1 | awk -F'_L' '{print $1}'| sed 's/ //g')
         
         if [ ! -d  SRR/${folder} ]
            then
                mkdir -p SRR/${folder} 
            fi
            
        
        
        ##change value of read if it already exist somewhere
        
                next=$(echo ${read1} | sed 's/L001/L002/g')
                nexter=$(echo ${read1} | sed 's/L001/L003/g')
                nexterer=$(echo ${read1} | sed 's/L001/L004/g')
                nextererer=$(echo ${read1} | sed 's/L001/L005/g')
                nexterererer=$(echo ${read1} | sed 's/L001/L006/g')

            if [ -f  "SRR/${folder}/${read1}" ] && [ ! -f  "SRR/${folder}/${next}" ];
            then
                    read1=$(echo ${read1}| sed 's/L001/L002/g')
                    read2=$(echo ${read2}| sed 's/L001/L002/g')
                    read3=$(echo ${read3}| sed 's/L001/L002/g')
                    read4=$(echo ${read4}| sed 's/L001/L002/g')
                
                
                elif [ -f  "SRR/${folder}/${read1}" ] && [ ! -f  "SRR/${folder}/${nexter}" ];
                    then
                    read1=$(echo ${read1}| sed 's/L001/L003/g')
                    read2=$(echo ${read2}| sed 's/L001/L003/g')
                    read3=$(echo ${read3}| sed 's/L001/L003/g')
                    read4=$(echo ${read4}| sed 's/L001/L003/g')
                    
                       
                    
                 elif [ -f  "SRR/${folder}/${read1}" ] && [ ! -f  "SRR/${folder}/${nexterer}" ];
                    then
                    read1=$(echo ${read1}| sed 's/L001/L004/g')
                    read2=$(echo ${read2}| sed 's/L001/L004/g')
                    read3=$(echo ${read3}| sed 's/L001/L004/g')
                    read4=$(echo ${read4}| sed 's/L001/L004/g')
                    
                      
                    
                 elif [ -f  "SRR/${folder}/${read1}" ] && [ ! -f  "SRR/${folder}/${nextererer}" ];
                    then
                    next1=$(echo ${read1}| sed 's/L001/L005/g')
                    next2=$(echo ${read2}| sed 's/L001/L005/g')
                    next3=$(echo ${read3}| sed 's/L001/L005/g')
                    next4=$(echo ${read4}| sed 's/L001/L005/g')
                    
                    
                 elif [ -f  "SRR/${folder}/${read1}" ] && [ ! -f  "SRR/${folder}/${nexterererer}" ];
                    then
                    next1=$(echo ${read1}| sed 's/L001/L006/g')
                    next2=$(echo ${read2}| sed 's/L001/L006/g')
                    next3=$(echo ${read3}| sed 's/L001/L006/g')
                    next4=$(echo ${read4}| sed 's/L001/L006/g')
                    
                      
            
            fi
        
        
        
        
        
        
        
        
        if [ "$replicate" ==  "$EXPERIMENT" ];
        then
            replicate="replicate 1"
        fi



         
         #echo SRR/${i} to SRR/${folder}/${i}
         
         
         #touch SRR/${folder}/${read4}

echo ${file},${AGE},${EXPERIMENT},${replicate},${read1},${read2},${read3},${read4} >> metadata_Bobadilla_2023.csv

##########DO THE FASTQDUMP ONLY IF IT HASNT BEEN DONE BEFORE
        if [ "$RUN" = "YES" ];  
        then
            if [ ! -f  "~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${read1}"  ] ;
            then
            
             printf "\n"
                echo "SAMPLE##################"
                echo ${i}
                
            mv SRR/${i} SRR/${folder}/${i}


             ~/sratoolkit.3.0.1-ubuntu64/bin/fastq-dump ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${i} --split-files --split-spot --disable-multithreading -O ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/  ;
            fi
            
                pigz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_1.fastq                   
                pigz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_2.fastq 
                pigz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_3.fastq 
                pigz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_4.fastq 
                   
                   
                   
                mv  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_1.fastq.gz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${read1}    
                mv  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_2.fastq.gz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${read2}
                mv  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_3.fastq.gz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${read3}
                mv  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_4.fastq.gz  ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${read4}
            

   if [ -f  "SRR/${folder}/${read1}" ] && [ -f  "SRR/${folder}/${read2}" ] && [ -f  "SRR/${folder}/${read3}" ] && [ -f  "SRR/${folder}/${read4}" ]
       then
        #rm ${p}
        #rm ~/data/mouse_dev_jaw/spatial_data/2023_Bobadilla/SRR/${folder}/${file}_?.fastq.gz
        echo YEAH IT WORKED FOR ${file}
       fi
       
         fi
         
done
wait