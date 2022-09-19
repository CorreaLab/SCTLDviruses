##This script is a backbone for conducting iterative processes.
##This version is specific for bbsplit.sh so all sample lists, echo strings, etc. are specific for this example.

#!/bin/bash


samplelist="cnat_sample.list" #Set variable samplelist to the unique sequencing IDs contained in this file
ref="rdrp.fasta" # Set variable ref to the fasta file containing the coding sequences

echo "bbsplit loop starting...." #echo progress printed to the terminal
for x in $(cat $samplelist) #start a loop file where x takes the value of each unique sequencing ID one at a time
do      echo "Starting with $x at $(date)" #echo progress printed to the terminal, indicating which sequencing file you are on
        R1=$(ls $x* | grep  "R1") #Define R1 variable. List both R1 and R2 for unique file ID, pipe these two to grep and use grep to select only the one that has R1 in it to assign R1 variable.
        R2=$(ls $x* | grep  "R2") #Repeat for R2 variable definition
        echo "input files = $R1 $R2" #echo progress printed to the terminal, check to see R1 and R2 are correct
        echo "reference file(s) = $ref" #echo prints variable ref to let user know correct reference file is used
        #seal.sh in1="$R1" in2="$R2" ref="$ref" stats="$x"_sealstats.txt rpkm="$x"_sealrpkm.txt out="$x"_chfv_R1.fq.gz out2="$x"_chfv_R2.fq.gz ambig=random
        bbsplit.sh in1="$R1" in2="$R2" ref="$ref" scafstats="$x"_scafstats.out basename="$x"_match_%_#.fq #Call bbsplit tool located in /home/ajv5/miniconda3/envs/vAMPirus/bin, unfortunately this is in Alex's miniconda installation so it's only accessible using his conda enviornment right now.
        #bbmerge-auto.sh in1="$R1" in2="$R2" out="$x"_merged.fq outu="$x"_umerged.fq extend2=20 iterations=5
        #cat "$x"_merged.fq "$x"_umerged.fq >> "$x"_all.fq
        #diamond blastx -q "$x"_all.fq -d $ref -p 30 --al "$x"_rdrp.fq --min-orf 15 --min-score 50 --ultra-sensitive -o "$x"_diamond_rdrp.out -f 6 qseqid qlen sseqid qstart qend qseq sseq length qframe evalue bitscore pident --max-targe$
        #pigz -p 30 "$x"_merged.fq
        #pigz -p 30 "$x"_umerged.fq
        #pigz -p 30 "$x"_all.fq
        echo "Finished with $x at $(date)"
        echo "    "
        echo "---------"
        echo "    "
done

cat *match_rdrp_1* >> all_rdrp_match_R1.fq
cat *match_rdrp_2* >> all_rdrp_match_R2.fq

rnaspades.py -1 all_rdrp_match_R1.fq -2 all_rdrp_match_R2.fq -o rdrp_spades -t 34

mv ./rdrp_spades/transcripts.fasta ./flsctld_rdrp.fasta
