##Download refseq genomes using manually created list of groups and unique ID
## Authors: Alex Veglia and Samantha Coy 05-11-22 version 1

#!/bin/bash


typelist="/home/ajv5/Databases/refseq/type.list" ##list of refseq groups we will loop through
listdir="/home/ajv5/Databases/refseq/uniquefilesnames" #You must create this list with the unique ID for each refseq group. We can write this into script later... used regex manually.
#amt="" ##number of files  ##this does not work if ends in anything but .1 ##bacteria.3012.1.genomic.fna.gz

URL="https://ftp.ncbi.nlm.nih.gov/refseq/release" #set the variable URL to the following string

# This is a loop for downloading refseq genomes
for i in $( cat $typelist ) #i variable takes the value of each line in type.list, downloads, then moves to next line in type.list to repeat
do
        #echo to command what is happening in real time to know progress
        echo "Downloading refseq file: ${i}" #Alex, I think we should remove curly brackets in lieu of $ capture, so it's not redundant?
        # Downloading refseq entries
        mkdir $i #make directory for the refseq group you are on in the loop from type.list
        cd ./"$i" #cd to the directory you just created
        for z in $(cat "$listdir"/"$i".txt) #Now we are going to start the download loop, referencing the uniquefilenames.txt file to set z to a new unique ID for each iteration of the loop
        do      wget "$URL"/"$i"/"$i"."$z".genomic.fna.gz # uses wget to download from url, that combines i ifor the refseq group and z for the unique ID .genomic.fna.gz file
                echo "Done downloading "$i"."$z".genomic.fna.gz" #echo tells you the progress in real time
                # Construct concatenated file of all
                zcat "$i"."$z".genomic.fna.gz >> "$i".fna #use zcat instead of cat to combine one file, cat can't handle compressed .gz file but zcat can
                # We don't need individual copy of each file that was created from NCBI, so remove it so all we have left i the concatenated file
                rm "$i"."$z".genomic.fna.gz #rm command to remove the file
        done # end of the download loop, returns to beginning of this loop to see if any other values are in uniquefilenames to apply
        pigz -p 33 "$i".fna #use pigz to compress file and outputs compressed version 
        cd .. ##changes directory to parent so we are prepared to start back at the first loop
done #end of first loop.
