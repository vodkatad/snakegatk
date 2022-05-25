#!/bin/bash

while read p; do
        url1=$(echo $p | tr -s " " "\t"| cut -f 11)
        url2=$(echo $p | tr -s " " "\t"| cut -f 12)
        f=$(echo $p | tr -s " " "\t"| cut -f 2)
        f1=${f}"_SA_L001_R1_001.fastq.gz"
        f2=${f}"_SA_L001_R2_001.fastq.gz"
        echo "#downloading $f"
        date
        echo -e "$url1\t$f1"
        echo -e "$url2\t$f2"
        wget -O $f1 $url1
        wget -O $f2 $url2
        date
done < $1
