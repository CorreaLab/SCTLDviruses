##this is the script used to rename flsctld rna seq lib
type="" ##dsRNA or polyA
infofile="" ##csvfile from seq facility

for x in *fastq.gz;
do      echo "$x"
        grep -w "$(echo "$x" | awk -F "_" '{print $2}')" $infofile | awk -F "," '{print $2}'
        samp="$(grep -w "$(echo "$x" | awk -F "_" '{print $2}')" $infofile | awk -F "," '{print $2}')"
        echo "samp = $samp"
        dir=$(echo "$x" | awk -F "_" '{print $5}')
        mv $x "$samp"_"$type"_"$dir".fastq.gz
        echo ""$x","$samp"_"$type"_"$dir".fastq.gz" >> new_names_"$type"libs_5_6_22.csv
done
