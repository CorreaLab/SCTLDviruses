#This script references paired end reads rogether and runs them through fastp with nextflow interfacing

#!/bin/bash

threads="25"    ##set number of threads to use
avqual="25"  ##Average read quality - forward or reverse reads will be discarded if average base quality across the read is below the number set below (25 is a good start)
maxn="3"     ##Maximum number of "N"s acceptable in the forward or reverse reads (default for fastp is 5)
trimq="15"  ##  Minmum base quality to be trimmed


echo "Starting fastp read processing for flsctld"

mkdir fastpreports
for x in *R1*
do      echo "Starting with $x at $(date)" #prints progress as updated for each read pair in real time
        name="$(echo $x | awk -F "_R1" '{print $1}')" #sets variable name to create a string combining temporary loop variable x and blah_R1.fastq using awk. Awk finds the first line of the file with {print $1}.
        R2="$(echo $x | sed 's/R1/R2/g')" #sets variable R2 to the R2 filename by simply replacing R1 with R2, leveraging paired end read duplicity. sed substitutes R1 with R2 globally.
        fastp -i $x -I $R2 -o "$name"_clean_R1.fq -O "$name"_clean_R2.fq --detect_adapter_for_pe --average_qual $avqual -n $maxn -c --overrepresentation_analysis --html "$name".fastp.html --json "$name".fastp.json --thread "$threads" --report_title "$name" #executes fastp, see next line for more annotation info
        #the fastp program executes quality control, adapter trimming, filtering, splitting, and merging
        echo "Finished with $x  at $(date)" #prints progress as updated for each read pair in real time
done

mv *html ./fastpreports/
mv *json ./fastpreports/
