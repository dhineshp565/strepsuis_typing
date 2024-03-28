
# Remove the first line from the file and process the remaining lines
sed '1d' $1 | while IFS=, read -r sample _; do
# Check if the serotype is cps-2 or cps-1 and the SNP is at position 483 
    if grep -E -wq "cps-2|cps-1" "${sample}_sero.csv" && cut -f 2 "${sample}_vcf.csv" | grep -wq "483" ; then
        # Replace the serotype-2 with serotype-1/2 
        sed -i 's,serotype-2,serotype-1/2,g' "${sample}_sero.csv"
        # Replace the serotype-14 with serotype-1
        sed -i 's,serotype-14,serotype-1,g' "${sample}_sero.csv"
        sed -i 's/_flye.fasta//g' "${sample}_sero.csv"
        cut -f 1,6,10,11,15 "${sample}_sero.csv" > "${sample}_serotype.csv"
    else
        # No snp at position 483
        sed -i 's,_flye.fasta,,g' "${sample}_sero.csv"
        cut -f 1,6,10,11,15 "${sample}_sero.csv" > "${sample}_serotype.csv"
    fi
done
