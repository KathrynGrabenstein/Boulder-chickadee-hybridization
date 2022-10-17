# Boulder Chickadee Study
### Whole genome analysis of 2018, 2019, 2020 samples
#### Kathryn Grabenstein
#### 17 Dec 2021

## Goal: Quantify hybridization between BCCH and MOCH in Boulder, CO using whole genome samples from 2018, 2019, & 2020 breeding season.


## Overview of workflow

### Samples have all already been demultiplexed, begin with trimming & qc
 
### Step 1: Trim & qc
### Step 2: Assembly and prep for variant calling
### Step 3: Call variants
### Step 4: merge, genotype, and filter
### Step 5: Calculate hybrid indexes for all individuals using fixed SNPs

#############################

## Step 1: Trim & qc

# Run trim script
bash trim-and-QC.sh -i 2020_Chickadee_WGS_Samples.txt -p /data1/kagr.wgs.chickadee/fastq/ -f R1_001.fastq.gz -r R2_001.fastq.gz -a TruSeq3-PE.fa -t 18

#check pre-trimmed and post-trimmed QC scores
#use multiqc

#download multiqc
#into working directory
pip install multiqc --user

#run multiqc
multiqc .


#############################
## Step 2: Align and Sort


##Open screen for this command so runs in the background
bash align-and-sort.sh -i chickadee.wgs.MASTER.txt -r /data1/kagr.wgs.chickadee/reference.BCCH.genome.2021/BCCH.reference.final.fasta -t 18 -p /data1/kagr.wgs.chickadee/trimmed_fastq/


#############################
## Need to merge bams for 2018 data
##Use merge_bams.sh


#copy 2018 bam files to directory
for file in $(cat chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/bam_files/"$file".bam /data1/kagr.wgs.chickadee/2018_bams_to_merge/; done

#copy 2018 file list to edit to remove lane info
cp chickadee.wgs.2018.txt chickadee.wgs.2018.no.lane.txt


#rename file
mv chickadee.wgs.2018.no.lane_1.txt chickadee.wgs.2018.no.lane.txt

#edit merge script for file paths for txt file and output directory in vim

#navigate to /data1/kagr.wgs.chickadee/bam_files
#run merge_bam.sh

#reconnect to screen to merge

bash merge_bams.sh

## rename merged bam files to not have "merged"
rename 's/merged//' *.bam


#re-loop through bams, and mark RGs and duplicates

### edit align-sort script to include only sorting, bc bam files already exist
cp align-and-sort.sh sort-bams-only.sh

bash sort-bams-only.sh -i bams.to.sort.2018.2020.txt -r ma -t 10 -b /data1/kagr.wgs.chickadee/bam_files/ -s /data1/kagr.wgs.chickadee/broken_pipe_sorted_bams/


## Run mpileup in parallel
### tyler chafin made script, move to working directory

## run mpileup script
bash samtools-snp-pipeline-scaffoldParallel.sh 




#count number of SNPs

#also use bcftools
bcftools stats called_snps_filtered.vcf > vcf.stats
cat vcf.stats

number of SNPs: 56172492    
number of indels:    0
number of multiallelic SNP sites:  0

bcftools stats called_raw_variants.vcf > vcf.stats.raw
cat vcf.stats.raw

number of SNPs: 79525688   
number of indels:   6784542
number of multiallelic SNP sites: 3664216 

#count number of SNPs
grep -v "^#" called_snps_filtered.vcf | wc -l
#56172492 lines


##########################
### Step three: Variant filtering

##calculate summary statistics for vcf table to get an idea of how to set filtering requirements

### compress vcf table bc MASSIVE
bgzip -c called_snps_filtered.vcf > called_snps_filtered.vcf.gz


#####################

## randomly subset VCF to speed up analyses
bcftools view called_snps_filtered.vcf.gz | vcfrandomsample -r 0.0012 > chickadee_subset.vcf

## compress subset vcf
bgzip chickadee_subset.vcf

# index vcf
bcftools index called_snps_filtered.vcf.gz

## set variables to save typing below
SUBSET_VCF=/data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz
OUT=/data1/kagr.wgs.chickadee/vcftools/chickadee_wgs

### Calculate allele frequency
vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2

### Calculate mean depth per individual
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT

### Calculate mean depth per variant
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT

### Calculate site quality score for each variant
vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT

### Calculate proportion missing data per inidividual
vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT

### Calculate missing data per site
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

#####################

### Step 4: merge, genotype, and filter
#### Decide how to filter SNPs 

 #what individuals have missing data?
mawk '$5 > 0.3' chickadee_wgs.imiss | cut -f1 > lowDP.indv.txt

VCF_IN=/data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz
VCF_OUT=/data1/kagr.wgs.chickadee/chickadee_wgs_FILTERED.vcf.gz

##set filters
###filter: 30 phred quality score, 0.1 minor allele freq, 90% missingness, max depth 10x
###no min depth bc all shallow depth

QUAL=80
MAF=0.05
MISS=0.50
MAX_DEPTH=10
MIN_DEPTH=2

#run vcftools filter command
#Filter #1
vcftools --gzvcf /data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz \
--remove-indels --maf 0.05 --max-missing 0.50 --minQ 80 \
--min-meanDP 2 \
--max-meanDP 10 \
--minDP 2 \
--maxDP 10 --recode --stdout | gzip -c > chickadee_wgs_FILTERED.vcf.gz 

#### Filter #2 
#run vcftools filter command
vcftools --gzvcf /data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz \
--remove-indels --maf 0.1 --max-missing 0.50 --minQ 80 \
--min-meanDP 2 \
--max-meanDP 10 \
--minDP 2 \
--maxDP 10 --recode --stdout | gzip -c > chickadee_wgs_FILTERED_2.vcf.gz 

#### Filter #3 
#run vcftools filter command
vcftools --gzvcf /data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz \
--remove-indels --maf 0.1 --max-missing 0.90 --minQ 80 \
--min-meanDP 2 \
--max-meanDP 10 \
--minDP 2 \
--maxDP 10 --recode --stdout | gzip -c > chickadee_wgs_FILTERED_3.vcf.gz 

#### Filter #4 (selected filter)
#run vcftools filter command
vcftools --gzvcf /data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz \
--remove-indels --maf 0.1 --max-missing 0.80 --minQ 80 \
--min-meanDP 2 \
--max-meanDP 10 \
--minDP 2 \
--maxDP 10 --recode --stdout | gzip -c > chickadee_wgs_FILTERED_4.vcf.gz 
### SNPs: 412351


## Calc Fst for alleles
##Calculate Fst per SNP to use informative SNPs to calculate hybrid indeces


##Calc by locus Fst
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz --weir-fst-pop allo_BCCH.txt --weir-fst-pop allo_MOCH.txt --out Fst_chickadee_wgs_filter_4

##calculate Fst between populations 
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz --weir-fst-pop allo_BCCH.txt --weir-fst-pop allo_MOCH.txt --out pop1_vs_pop2
##mean fst between allopatric pops: 0.32


##Output txt with loci with only alleles with 0.65 or greater
awk '$3 >= 0.65' Fst_chickadee_wgs_filter_4.weir.fst | awk '!/-nan/' > Loci_Fst_0.65_chickadee_wgs_filter_4.txt

#how many loci with fst 0.65 or greater?
wc -l Loci_Fst_0.65_chickadee_wgs_filter_4.txt
#70423 loci Fst >0.65

##output txt with loci with only alleles with 0.80 or greater, ignore lines with -nan (missing values)
awk '$3 >= 0.80' Fst_chickadee_wgs_filter_4.weir.fst | awk '!/-nan/' > Loci_Fst_0.80_chickadee_wgs_filter_4.txt

#how many loci with fst 0.65 or greater?
wc -l Loci_Fst_0.80_chickadee_wgs_filter_4.txt
#46443 loci >0.80

##output txt with loci with only fixed alleles
awk '$3 >= 1' Fst_chickadee_wgs_filter_4.weir.fst | awk '!/-nan/' > Loci_Fst_1_chickadee_wgs_filter_4.txt

#how many fixed alleles?
wc -l Loci_Fst_1_chickadee_wgs_filter_4.txt
#22506 loci Fst = 1


### Weir and Cockerham mean Fst estimate: 0.37572
### Weir and Cockerham weighted Fst estimate: 0.55881


#Filter VCF.4 with position file for Fst 0.65
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz \
--positions Loci_Fst_0.65_chickadee_wgs_filter_4.txt \
--recode --stdout | gzip -c > \
VCF_Loci_0.65_chickadee_wgs_FILTERED_2.vcf.gz


#Filter VCF.4 with position file for Fst 0.8
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz \
--positions Loci_Fst_0.8_chickadee_wgs_filter_4.txt \
--recode --stdout | gzip -c > \
VCF_Loci_0.8_chickadee_wgs_FILTERED_2.vcf.gz

#Filter VCF.4 with position file for fixed loci
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz \
--positions Loci_Fst_1_chickadee_wgs_filter_4.txt \
--recode --stdout | gzip -c > \
VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf.gz



# Need to make triangle plots for whole genome 
## Calculate heterozygosity and combine with gghyhbrid hybrid index

###To calculate hybrid index as an average genotype and to compare to heterozygosity 
###and to gghybrid index, you first need to exclude non-reference heterozygote 
###in your “zero” allopatric population. If you open your vcf table (the one you created for 
###loci fixed to alternative states (Fst=1) in allopatric populations), 
###you’ll see that one of your allopatric populations will have almost all genotypes as “0”, 
###while only few will be “2”. So you need to trim those “2” away.


# output your allopatric populations and look which has the majority of refence genotypes (“0”):
vcftools --gzvcf VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf.gz --keep allo_BCCH.txt --recode --out bcch.only
vcftools --gzvcf VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf.gz --keep allo_MOCH.txt --recode --out moch.only


# output loci with non-reference allele frequency above 0.9
vcftools --vcf bcch.only.recode.vcf --non-ref-af 0.9 --recode --out bcch.only.non.ref.loci
#After filtering, kept 386 out of a possible 22505 Sites


# extract locus names
vcftools --vcf bcch.only.non.ref.loci.recode.vcf --012 --out bcch.only.non.ref.loci.names


# recode the main file and exclude those non-ref loci
vcftools --gzvcf VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf.gz --exclude-positions bcch.only.non.ref.loci.names.012.pos --recode --out complete.non.ref.excl
#kept 22119 loci

# output file as 012
vcftools --vcf complete.non.ref.excl.recode.vcf --012 --out complete.non.ref.excl

