#/bin/#!/usr/bin/env bash


type="" ##bacter,arch,vert
amt="" ##number of files
 ##bacteria.3012.1.genomic.fna.gz

for z in $(seq "$amt");
do  curl -o "$type"."$z".1.genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria."$z".1.genomic.fna.gz    
done
