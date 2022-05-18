#!/bin/bash




samplelist="cnat_sample.list"
ref="RdRP_viruses_900-12000_cds.fasta"

for x in $(head -1 $samplelist)
do      echo "Starting with $x at $(date)"
        R1=$(ls $x* | grep  "R1")
        R2=$(ls $x* | grep  "R2")
        echo "input files = $R1 $R2"
        echo "reference file(s) = $ref"
        #seal.sh in1="$R1" in2="$R2" ref="$ref" stats="$x"_sealstats.txt rpkm="$x"_sealrpkm.txt out="$x"_chfv_R1.fq.gz out2="$x"_chfv_R2.fq.gz ambig=random
        bbsplit.sh in1="$R1" in2="$R2" ref="$ref" pattern="$x"_chfv_%.fq scafstats="$x"_scafstats.out #out="$x"_chfv_R1.fq.gz out2="$x"_chfv_R2.fq.gz
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
