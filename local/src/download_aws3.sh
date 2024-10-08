#!/bin/bash

while read p; do
        url1=$(echo $p | tr -s " " "\t"| cut -f 9)
        url2=$(echo $p | tr -s " " "\t"| cut -f 10)
        f=$(echo $p | tr -s " " "\t"| cut -f 1)
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
