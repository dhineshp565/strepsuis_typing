#!/bin/env bash

# Remove the first line from the file and process the remaining lines
sed '1d' $1 | while IFS=, read -r sample _; do
    if grep -E -wq "cps-2|cps-1" "${sample}_sero.csv" && cut -f 2 "${sample}_vcf.csv" | grep -wq "483" ; then
        sed -i 's,serotype-2,serotype-1/2,g' "${sample}_sero.csv"
        sed -i 's,serotype-14,serotype-1,g' "${sample}_sero.csv"
        sed -i 's/_flye.fasta//g' "${sample}_sero.csv"
        cut -f 1,6,10,11,15 "${sample}_sero.csv" > "${sample}_serotype.csv"
    else
        sed -i 's,_flye.fasta,,g' "${sample}_sero.csv"
        cut -f 1,6,10,11,15 "${sample}_sero.csv" > "${sample}_serotype.csv"
    fi
done
