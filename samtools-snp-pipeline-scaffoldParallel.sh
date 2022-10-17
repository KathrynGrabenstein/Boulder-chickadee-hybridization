#!/bin/bash

# Call snps in samtools

ref="BCCH_sorted.nuc.fixed.fasta"
bamdir="/data1/kagr.wgs.chickadee/sorted_bam_files/"
parallel_bin="/home/tchafin/local/bin/parallel"
threads=20
scaffold_vcf_dir="called_scaffolds_raw"

ID="called" # This will be used as a prefix for the output file

echo "making a pileup file for" $ID >> log

mkdir $scaffold_vcf_dir

call_scaffold(){
   scaffold=$1
   bamdir=$2
   ref=$3
   scaffold_vcf_dir=$4
        	
   bcftools mpileup -d 8000 -Ou -f $ref --ignore-RG -a AD,ADF,DP,SP,INFO/AD,INFO/ADF -r $1 "$bamdir"*.bam | bcftools call -mv -O z -o $scaffold_vcf_dir'/'$scaffold'.vcf.gz'
   bcftools index $scaffold_vcf_dir'/'$scaffold'.vcf.gz'
}
export -f call_scaffold

#can also add the -R flag joined with a scaffold list to subset and parralel
echo "Running samples from $bamdir in parallel for each scaffold, across $threads threads" >>log
cat $ref | grep ">" | sed "s/>//g" | $parallel_bin -j $threads call_scaffold {} $bamdir $ref $scaffold_vcf_dir

#merge sample VCFs
echo "Merging sample VCF files for $ID" >> log
bcftools concat -O v -o $ID"_raw_variants.vcf" $scaffold_vcf_dir/*.vcf.gz

#echo "removing all lines with two comment marks" >> log
#grep -v "##" "$ID"_snps_indels.vcf > "$ID"_snps_indels_short.vcf
echo "filtering low quality snps (<100)" >> log
awk '$1~/^#/ || $6 > 100 {print $0}' > \
"$ID"_snps_indels_filtered.vcf "$ID"_raw_variants.vcf
echo "add the header and check length of column 4 and 5 to make sure they are snp type variants" >> log
awk '$1~/^#/ || length($4)==1 && length($5)==1 {print $0}'> \
"$ID"_snps_filtered.vcf "$ID"_snps_indels_filtered.vcf
