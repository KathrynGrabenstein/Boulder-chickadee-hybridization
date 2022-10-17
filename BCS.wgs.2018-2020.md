# Boulder Chickadee Study
# Whole genome analysis of 2018, 2019, 2020 samples
# Kathryn Grabenstein
# 17 Dec 2021

#Goal: Quantify hybridization between BCCH and MOCH in Boulder, CO using whole genome samples from 2018, 2019, & 2020 breeding season.



## Overview of workflow

# Samples have all already been demultiplexed, begin with trimming & qc
 
# Step 1: Trim & qc
# Step 2: Assembly and prep for variant calling
# Step 3: Call variants
# Step 4: merge, genotype, and filter
# Step 5: Calculate hybrid indexes for all individuals using fixed SNPs

#General how tos
########################
#open up new screen
screen

#escape screen
CTRL + AD

#see running screens
screen -list

#reconnect to screen
screen -r [INSERT SCREEN ID]

#kill screen
#reconnect to screen and then type "exit"

#kill screen
screen -X -S [session # you want to kill] quit


#Check data drive space
df -h /data1/
df -h /nas/

#############################

##### Step 1: Trim & qc

# Move all scripts from local computer to working directory
#navigate to local computer on terminal
scp /Users/kathryngrabenstein/Documents/Data Analysis/Chickadee Whole Genome/ *.sh kagr3234@chickadee:/data1/kagr.wgs.chickadee

# Move adapter sequence file from local computer to working directory
#navigate to local computer on terminal
scp /Users/kathryngrabenstein/Documents/Data Analysis/Chickadee \Whole \Genome/TruSeq3-PE.fa kagr3234@chickadee:/data1/kagr.wgs.chickadee

# Move samples from nas to data1
## move fastq files to working directory
## move into directory where files are located

#move 2018 samples
cp -r *.gz /data1/kagr.wgs.chickadee/fastq

# dont move 2019 files bc already have trimmed files
#move trimmed 2019 files to trimmed_fastq folder

#move 2020 raw files
cp **/*.gz /data1/kagr.wgs.chickadee/fastq

#get rid of A prefix for 2020 samples
rename 's/A*//' *.fq.gz

#make sure naming system consistent for 2018 and 2020
rename 's/_1.fq.gz/_R1_001.fastq.gz/g' *_1.fq.gz

rename 's/_2.fq.gz/_R2_001.fastq.gz/g' *_2.fq.gz

# Generate list of samples to trim
#from fastq directory
ls -1 > 2020_Chickadee_WGS_Samples.txt

#move txt file to working directory
mv 2020_Chickadee_WGS_Samples.txt /data1/kagr.wgs.chickadee/

#remove prefixes from sample txt file for trimmomatic
sed -i 's/_R1_001.fastq.gz//g' 2020_Chickadee_WGS_Samples.txt
sed -i 's/_R2_001.fastq.gz//g' 2020_Chickadee_WGS_Samples.txt

#remove duplicates
uniq 2020_Chickadee_WGS_Samples.txt > 2020_Chickadee_WGS_Samples_1.txt

#rename file
mv 2020_Chickadee_WGS_Samples_1.txt 2020_Chickadee_WGS_Samples.txt

# Run trim script
bash trim-and-QC.sh -i 2020_Chickadee_WGS_Samples.txt -p /data1/kagr.wgs.chickadee/fastq/ -f R1_001.fastq.gz -r R2_001.fastq.gz -a TruSeq3-PE.fa -t 18

#check pre-trimmed and post-trimmed QC scores
#use multiqc

#download multiqc
#into working directory
pip install multiqc --user

#run multiqc
multiqc .

#move multiqc file to local computer
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/post_trim_QC_files/multiqc_report.html /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/

scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/pre_trim_QC_files/multiqc_report.html /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/


#############################
#Step 2: Align and Sort

##prep work for align and sort; ie manage all the new trimmed files
#copy trimmed 2018 & 2020 back to raw folders to save for previous use (ie if updating ref in future)
#move trimmed to master trimmed directory

##### copy 2018 trimmed to master 2018 folder
#make file with all sample names for 2018
#navigate to 2018 directory
ls > chickadee.wgs.2018.txt

#change suffixes to match trimmed files
sed -i 's/_R2_001.fastq.gz//g' chickadee.wgs.2018.txt
sed -i 's/_R1_001.fastq.gz//g' chickadee.wgs.2018.txt

#remove duplicates (shouldn't be any?)
uniq chickadee.wgs.2018.txt > chickadee.wgs.2018_1.txt

#rename file
mv chickadee.wgs.2018_1.txt chickadee.wgs.2018.txt

#copy ALL trimmed files to 2018/trimmed master directory for posterity 
#four files for every sample
for file in $(cat /nas/chickadee.wgs/2018/chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_2P.fq.gz /nas/chickadee.wgs/2018/trimmed/; done

for file in $(cat /nas/chickadee.wgs/2018/chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_1P.fq.gz /nas/chickadee.wgs/2018/trimmed/; done

for file in $(cat /nas/chickadee.wgs/2018/chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_2U.fq.gz /nas/chickadee.wgs/2018/trimmed/; done

for file in $(cat /nas/chickadee.wgs/2018/chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_1U.fq.gz /nas/chickadee.wgs/2018/trimmed/; done



###### copy 2020 trimmed
#make sample list for 2020
ls > chickadee.wgs.2020.txt

#get rid of A prefix for 2020 samples
sed -i 's/A*//' chickadee.wgs.2020.txt

#copy ALL trimmed 2020 files to master 2020/trimmed directory for posterity
#four files for every sample

for file in $(cat /nas/chickadee.wgs/2020/chickadee.wgs.2020.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_2P.fq.gz /nas/chickadee.wgs/2020/trimmed/; done

for file in $(cat /nas/chickadee.wgs/2020/chickadee.wgs.2020.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_1P.fq.gz /nas/chickadee.wgs/2020/trimmed/; done

for file in $(cat /nas/chickadee.wgs/2020/chickadee.wgs.2020.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_2U.fq.gz /nas/chickadee.wgs/2020/trimmed/; done

for file in $(cat /nas/chickadee.wgs/2020/chickadee.wgs.2020.txt); do scp /data1/kagr.wgs.chickadee/fastq/"$file"_trimmed_1U.fq.gz /nas/chickadee.wgs/2020/trimmed/; done


##now that everything has been "backed up", move trimmed fastq to trimmed_fastq with 2019 data

#move trimmed 2018 & 2020 from fastq to trimmed_fastq
mv *fq.gz /data1/kagr.wgs.chickadee/trimmed_fastq/


# make MASTER sample file with 2018, 2019 & 2020 data
#navigate to trimmmed_fastq
ls > chickadee.wgs.MASTER.txt

#remove suffixes 1P, 2P, 1U, 2U so only sample name 
sed -i 's/_trimmed_*1P.fq.gz//g' chickadee.wgs.MASTER.txt
sed -i 's/_trimmed_*2P.fq.gz//g' chickadee.wgs.MASTER.txt
sed -i 's/_trimmed_*1U.fq.gz//g' chickadee.wgs.MASTER.txt
sed -i 's/_trimmed_*2U.fq.gz//g' chickadee.wgs.MASTER.txt


#remove duplicates (shouldn't be any?)
uniq chickadee.wgs.MASTER.txt > chickadee.wgs.MASTER_1.txt

#rename file
mv chickadee.wgs.MASTER_1.txt chickadee.wgs.MASTER.txt

#move sample file to working directory
mv chickadee.wgs.MASTER.txt /data1/kagr.wgs.chickadee/

##Open screen for this command so runs in the background
bash align-and-sort.sh -i chickadee.wgs.MASTER.txt -r /data1/kagr.wgs.chickadee/reference.BCCH.genome.2021/BCCH.reference.final.fasta -t 18 -p /data1/kagr.wgs.chickadee/trimmed_fastq/


#############################
## Need to merge bams for 2018 data
##Use merge_bams.sh


#copy 2018 bam files to new directory in case make mistake, so don't have to remake bam files
#make directory for bam to merge
mkdir 2018_bams_to_merge

#copy 2018 bam files to directory
for file in $(cat chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/bam_files/"$file".bam /data1/kagr.wgs.chickadee/2018_bams_to_merge/; done

#copy 2018 file list to edit to remove lane info
cp chickadee.wgs.2018.txt chickadee.wgs.2018.no.lane.txt

#edit 2018 file list to use with merge_bams.sh to include ONLY sample ID, no lane info
sed -i 's/_S.*//g' chickadee.wgs.2018.no.lane.txt

#remove duplicates (shouldn't be any?)
uniq chickadee.wgs.2018.no.lane.txt > chickadee.wgs.2018.no.lane_1.txt

#rename file
mv chickadee.wgs.2018.no.lane_1.txt chickadee.wgs.2018.no.lane.txt

#edit merge script for file paths for txt file and output directory in vim

#navigate to /data1/kagr.wgs.chickadee/bam_files
#run merge_bam.sh

#reconnect to screen to merge

bash merge_bams.sh

## rename merged bam files to not have "merged"
rename 's/merged//' *.bam

## move merged bams to bam_files directory
for file in $(cat chickadee.wgs.2018.no.lane.txt); do mv /data1/kagr.wgs.chickadee/merged_2018_bams/"$file".bam /data1/kagr.wgs.chickadee/bam_files/; done


## Move allopatric BCCH files to sorted_bams 
for file in $(cat allo_BCCH.txt); do scp /nas/usftp21.novogene.com/raw_data/MO_chickadee.transect.raw/consolidated/sorted_bam_files/"$file"* /data1/kagr.wgs.chickadee/sorted_bam_files/; done



#################################
# need to make more space on data drive
### remove sorted bams from 2018 bc merged bams will need to be resorted
rm  *L004_sorted_RGadded_dupmarked.bam.bai
rm *L004_dupmarked_metrics.txt
rm *L004_sorted_RGadded_dupmarked.bam
rm  *L003_sorted_RGadded_dupmarked.bam.bai
rm *L003_dupmarked_metrics.txt
rm *L003_sorted_RGadded_dupmarked.bam

##make txt file with 2019 samples
ls -1 > chickadee.wgs.2019.txt

##remove prefixes from sample txt file for trimmomatic
sed -i 's/_trimmed.*//g' chickadee.wgs.2019.txt

##remove duplicates
uniq chickadee.wgs.2019.txt > chickadee.wgs.2019_1.txt

##rename file
mv chickadee.wgs.2019_1.txt chickadee.wgs.2019.txt

## move bam files from 2019 to nas to make more space on data1
for file in $(cat /nas/chickadee.wgs/2019/trimmed/chickadee.wgs.2019.txt); do scp /data1/kagr.wgs.chickadee/bam_files/"$file".bam /nas/chickadee.wgs/2019/bam_files; done

## remove 2019 bam files from data1 after moved
for file in $(cat /nas/chickadee.wgs/2019/trimmed/chickadee.wgs.2019.txt); do rm /data1/kagr.wgs.chickadee/bam_files/"$file".bam; done

##move 2018 bam files to nas for more space
for file in $(cat /nas/chickadee.wgs/2018/chickadee.wgs.2018.txt); do scp /data1/kagr.wgs.chickadee/bam_files/"$file".bam /nas/chickadee.wgs/2018/bam_files; done

## remove 2018 bam files from data1 after moved
for file in $(cat /nas/chickadee.wgs/2018/chickadee.wgs.2018.txt); do rm /data1/kagr.wgs.chickadee/bam_files/"$file".bam; done

## move merged bams to bam_files
mv *.bam /data1/kagr.wgs.chickadee/bam_files/


## some 2019 files skipped (the poxes); copy bam to bam_files on data1
for file in $(cat /data1/kagr.wgs.chickadee/bams.to.sort.2018.2020.txt); do scp /nas/chickadee.wgs/2019/bam_files/"$file".bam /data1/kagr.wgs.chickadee/bam_files/; done


## move trimmed_fastq files to nas to make more space on data1
mv /data1/kagr.wgs.chickadee/fastq/*.fq.gz /nas/chickadee.wgs/2020/trimmed/ 


### Space issues resolved, for now

#re-loop through bams, and mark RGs and duplicates

## combine broken pipe & merged bam file & sort both at the same time
### only sorting bam files, bc bam files already exist 

### copy 2018 file and broken pipe file to combine
cp chickadee.wgs.2018.no.lane.txt chickadee.wgs.2018.no.lane_1.txt
cp broken.pipe.chickadee.wgs.txt broken.pipe.chickadee.wgs_1.txt

cat chickadee.wgs.2018.no.lane_1.txt broken.pipe.chickadee.wgs_1.txt > bams.to.sort.2018.2020.txt

rm chickadee.wgs.2018.no.lane_1.txt
rm broken.pipe.chickadee.wgs_1.txt

### edit align-sort script to include only sorting, bc bam files already exist
cp align-and-sort.sh sort-bams-only.sh

bash sort-bams-only.sh -i bams.to.sort.2018.2020.txt -r ma -t 10 -b /data1/kagr.wgs.chickadee/bam_files/ -s /data1/kagr.wgs.chickadee/broken_pipe_sorted_bams/



##aaaaand we're back on track

# Prep for allopatric MOCH files
## Download from google drive

## move fastq files from local computer to server
### navigate to local computer on terminal
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/Allo_MOCH_Fastq/*.fastq.gz kagr3234@chickadee:/data1/kagr.wgs.chickadee/allo_MOCH/


# Run trim script
bash trim-and-QC.sh -i allo_MOCH_samples.txt -p /data1/kagr.wgs.chickadee/allo_MOCH/ -f R1_001.fastq.gz -r R2_001.fastq.gz -a TruSeq3-PE.fa -t 18

## align and sort allopatric MOCH files
bash align-and-sort.sh -i allo_MOCH_samples.txt -r /data1/kagr.wgs.chickadee/reference.BCCH.genome.2021/BCCH.reference.final.fasta -t 18 -p /data1/kagr.wgs.chickadee/allo_MOCH/

## move all allo_moch contents to sorted bam
cp -a /data1/kagr.wgs.chickadee/allo_MOCH/sorted_bam_files/. /data1/kagr.wgs.chickadee/sorted_bam_files/

## move hybrid files over 
cp -a /home/anha2922/hybrid/. /data1/kagr.wgs.chickadee/sorted_bam_files/


## Run mpileup in parallel
### tyler chafin made script, move to working directory

### move script from local computer to working directory
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/samtools-snp-pipeline-scaffoldParallel.sh kagr3234@chickadee:/data1/kagr.wgs.chickadee/

### make script executable
chmod +x samtools-snp-pipeline-scaffoldParallel.sh 


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

# get sample names
bcftools query -l /data1/kagr.wgs.chickadee/called_snps_filtered.vcf > wgs_chickadee_sample_names.txt

## move sample name file to local computer to make population ID file
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/wgs_chickadee_sample_names.txt /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/

##########################
#Variant filtering
## the most tedious and least fulfilling step :(

##calculate summary statistics for vcf table to get an idea of how to set filtering requirements

### compress vcf table bc MASSIVE
bgzip -c called_snps_filtered.vcf > called_snps_filtered.vcf.gz


#####################
# the following is in a script bc I didn't want to deal

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

# End of script
#####################

## switch to R and plot to figure out filtering parameters
###first, transfer vcftools files to local desktop
###run command from local terminal, not remote into server
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/vcftools/chickadee_wgs.* /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/

## now switch to R to plot, come back after figuring out how to filter


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

#how many variants remain after filtering?
bcftools view -H chickadee_wgs_FILTERED.vcf.gz | wc -l
 
##re check stats after filtering
### run bash script
bash vcf_stats.sh

## switch to R and plot to figure out filtering parameters
###first, transfer vcftools files to local desktop
###run command from local terminal, not remote into server
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/vcftools/chickadee_wgs_Filter_1* /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/

## now switch to R to plot, come back after figuring out how to filter

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

#### Filter #4 
#run vcftools filter command
vcftools --gzvcf /data1/kagr.wgs.chickadee/called_snps_filtered.vcf.gz \
--remove-indels --maf 0.1 --max-missing 0.80 --minQ 80 \
--min-meanDP 2 \
--max-meanDP 10 \
--minDP 2 \
--maxDP 10 --recode --stdout | gzip -c > chickadee_wgs_FILTERED_4.vcf.gz 
#SNPs: 412351

##re check stats after filtering
### run bash script
<!--#### VCFIN = chickadee_wgs_FILTERED_2.vcf.gz -->
#### VCF OUT = chickadee_wgs_Filter_2
bash vcf_stats.sh


## move population file from local computer to data1
#navigate to local computer on terminal
scp wgs_chickadee_pop_ID.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/


## move pdf from data1 to local computer
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/PCA.Whole.Genome.nolabels.pdf /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/

## plot only parentals
### exclude nestlings

#remove nestlings from analysis
#move txt file with nestlings to exclude from local desktop to remote
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/Nestlings_to_Exclude.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/

#move to Whole Genome working directory

VCF_IN=/data1/kagr.wgs.chickadee/chickadee_wgs_FILTERED.vcf
VCF_OUT=/data1/kagr.wgs.chickadee/chickadee_wgs_FILTERED.nestlings.excluded.vcf.gz

#run vcftools filter command
vcftools --vcf $VCF_IN \
--remove Nestlings_to_Exclude.txt \
--recode --stdout | gzip -c > \
$VCF_OUT





## Calc Fst for alleles
##Calculate Fst per SNP to use informative SNPs to calculate hybrid indeces


##Make list of allo BCCH individuals
#Used regular expressions
#did this in excel bc I don't know how to do in terminal

##Make list of MOCH individuals
#same, did this in excel


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


#Weir and Cockerham mean Fst estimate: 0.37572
#Weir and Cockerham weighted Fst estimate: 0.55881


#move txt files to local computer
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/Loci_Fst_*_chickadee_wgs_filter_4.txt /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/

#create txt files with chromosome and position info for 0.65, 0.80, 1
#edit in text editor to remove Fst Values
#use regular expression find, replace to delete Fst column
\s0.\d+
#current file structre is chromosome (unchanged), TAB, and position
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/Loci_Fst_*_chickadee_wgs_filter_4.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/


#Filter VCF.4 with position file for 0.65
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz \
--positions Loci_Fst_0.65_chickadee_wgs_filter_4.txt \
--recode --stdout | gzip -c > \
VCF_Loci_0.65_chickadee_wgs_FILTERED_2.vcf.gz


#Filter VCF.4 with position file for 0.8
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz \
--positions Loci_Fst_0.8_chickadee_wgs_filter_4.txt \
--recode --stdout | gzip -c > \
VCF_Loci_0.8_chickadee_wgs_FILTERED_2.vcf.gz

#Filter VCF.4 with position file for fixed
vcftools --gzvcf chickadee_wgs_FILTERED_4.vcf.gz \
--positions Loci_Fst_1_chickadee_wgs_filter_4.txt \
--recode --stdout | gzip -c > \
VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf.gz



#save with .STR ending so gghybrid recognizes them
#add .STR to end of file name
#click to open with TextWrangler

#Use Find & Replace to get rid of long file name from directory stuff
#need to create populations file where each pop is duplicated
#bc structure files have each line repeated twice bc diploid
#Use code below in excel
=OFFSET(Sheet1!$A$1,INT((ROWS($A$1:A1)-1)/2),)

#then need to add pop file as second column
#probs easiest to do in excel where can copy& paste whole column
#save STR as .csv to facilitate this

#rerun gghybrid with these STR files
#switch to R script now





## Need to convert vcf to STR file in R on server
library(vcfR)
vcf <- read.vcfR(VCF_Loci_1_chickadee_wgs_FILTERED_2.vcf, verbose = FALSE )
loci <- extract.gt(vcf, element = "GT", as.numeric = TRUE, )

pop.id <- read.delim(wgs_chickadee_pop_ID.txt, header = TRUE, sep = "\t", dec = ".", ...)


## recode vcf table to genotype table 
vcftools --vcf VCF_Loci_1_chickadee_wgs_FILTERED_2.vcf --012 

vcftools --gzvcf VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf.gz --012 


#### need to create diploid input file from 012 file made with vcftools


# remove first column with sample order
cut -f2- out.012 > temp1.txt

cut -f3- ready.input.chickadee.txt > temp2.txt

tail -n +2 temp2.txt > temp3.txt


# replace 012 with diploid genotypes

sed -i -e 's/-1/M M/g' temp3.txt
sed -i -e 's/0/A A/g' temp3.txt
sed -i -e 's/1/A T/g' temp3.txt
sed -i -e 's/2/T T/g' temp3.txt
sed -i -e 's/A/1/g' temp3.txt
sed -i -e 's/T/2/g' temp3.txt
sed -i -e 's/M/-9/g' temp3.txt
awk  '{gsub(" ","\t",$0); print;}' temp3.txt > chickadee.wgs.diploid.3.txt


## add scaffold names back to first row of datasheet

cat scaffolds.txt chickadee.wgs.diploid.txt > INPUT.chickadee.wgs.diploid.txt


scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/out.012 /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/out.012.pos /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/wgs_gghybrid_loci_1.csv /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/Rplots.pdf /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/big.test.ploidy2..txt /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/


scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/wgs_chickadee_duplicated.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/wgs_chickadee_pop_ID.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/out.012.pos kagr3234@chickadee:/data1/kagr.wgs.chickadee/
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/gghybrid_chickadee.R kagr3234@chickadee:/data1/kagr.wgs.chickadee/
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/ready.input.chickadee.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/

Rscript gghybrid_chickadee.R 




pdf(file = paste0(outDir,"/",name,"gghybrid_wgs_chickadee_loci_1.pdf"), useDingbats=FALSE)

setkey(hindlabel$hi,POPID)    #function from data.table, for rapid sorting and subsetting#


abc = plot_h(data=hindlabel$hi[c("Allo_BCCH","BCCH","Hybrid","MOCH","Allo_MOCH")],
             test.subject=hindlabel$test.subject,
             mean.h.by="POPID",            #Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),    #Order test subjects along the x axis by the mean hybrid 
                                                                #index calculated above and also by individual hybrid index
             col.group="POPID",
             group.sep="POPID",
             fill.source=TRUE,
             basic.lines=FALSE,
             source.col=c("red","blue"),
             source.limits=c("red","blue"),
             cex=1,pch=16,
             cex.lab=1.5,cex.main=1.5,ylim=c(0,1))

#Add a legend using the 'plot_h' object 'abc'#

setkey(abc,rn)        #Order data by row number#
legend("topleft",    #Place the legend in the top left of the figure#
       abc[,POPID],         #Name of the field by which data point colours are grouped#
       bg="white",            #Background colour#
       text.col=c("black"), #Text colour#
       pch=22,                 #Text size#
       col=abc[,col.Dark2], #Name of the field containing colour information#
       pt.bg=abc[,col.Dark2],    #Name of the field containing colour information#
       ncol=2,                #Number of columns for the group names#
       cex=1, pt.cex=1)

dev.off()


##ahhhhh!! 


9 total hybrids
#two are nestlings from the same nest (out of 5 chicks, in a BCCH nest in Louisville; sex unknown)
#one BCCH caught at SGR (ie overlap); had a BP on 10 june (ie breeding female)
#One M BCCH at watershed school, provisioned chicks equally as female
#One male BCCH at Stacey wolffe (ie out east), reared chicks
#4 nestlings out of 6 chicks, nest at watershed school, chicks of the male above^^






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

#move file to local computer
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/complete.non.ref.excl.012 /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome


##get sample names for meta merging
bcftools query -l /data1/kagr.wgs.chickadee/complete.non.ref.excl.recode.vcf > wgs_sample_names_triangle.txt

#move file to local computer
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/wgs_sample_names_triangle.txt /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome


#Run NewHybrids to assign hybrid classes
## Need vcf table with allo only

##filter VCF to only include allopatric individuals
VCF_IN=/data1/kagr.wgs.chickadee/VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf
VCF_OUT=/data1/kagr.wgs.chickadee/allopatric_1_chickadee_wgs.vcf

#run vcftools filter command
vcftools --vcf $VCF_IN \
--keep allo_wgs.txt \
--recode --stdout | gzip -c > \
$VCF_OUT

#check that filtered properly
##get sample names for meta merging
bcftools query -l /data1/kagr.wgs.chickadee/allopatric_1_chickadee_wgs.vcf > wgs_sample_names_allo.txt

scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/allopatric_1_chickadee_wgs.vcf /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome

## Need to filter the vcf bc too many loci for newhybrids to run quickly
#### Filter #4 
#run vcftools filter command
vcftools --vcf /data1/kagr.wgs.chickadee/VCF_Loci_1_chickadee_wgs_FILTERED_4.vcf \
--remove-indels --maf 0.1 --max-missing 0.85 --minQ 80 \
--min-meanDP 2 \
--max-meanDP 10 \
--minDP 2 \
--maxDP 10 --recode --stdout | gzip -c > Loci_1_chickadee_wgs_FILTERED_5.vcf.gz
 
 ## count number snps now
 grep -v "^#" Loci_1_chickadee_wgs_FILTERED_5.vcf | wc -l

#SNPs: 143

#move to local computer
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/sex.zchr.recode.vcf /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome



#### Need to sex chickadees using genomics

## filter by scaffold
vcftools —vcf filename --chr --recode --out

##move Z marker file from local comp to data1/kagr
scp /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome/scaffolds.on.z.txt kagr3234@chickadee:/data1/kagr.wgs.chickadee/


##take VCF, filter to keep markers on Z chromosome
vcftools --vcf called_raw_variants.vcf --positions scaffolds.on.z.txt --recode --out sex.zchr

vcftools --vcf called_snps_indels_filtered.vcf --positions scaffolds.on.z.txt --recode --out sex.zchr.2


# output file as 012
vcftools --vcf sex.zchr.recode.vcf --012 --out sex.zchr

##move file to local computer and switch to R
scp kagr3234@chickadee:/data1/kagr.wgs.chickadee/sex.zchr.012 /Users/kathryngrabenstein/Documents/Data\ Analysis/Chickadee\ Whole\ Genome

## Figure out what scaffolds have in the big unfiltered called SNPs vcf 

##zip vcf 
bgzip called_raw_variants.vcf

##create list of scaffold names in VCF to see if Z chromosome scaffolds included...
tabix called_raw_variants.vcf.gz
zcat called_raw_variants.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names​





## figure out what scaffolds are in vcf table 
samtools faidx BCCH_pseudochr_rename.fasta 'Pseudogi|224381666|ref|NC_011462.1|' 'Pseudogi|224381667|ref|NC_011474.1|' 'Pseudogi|224381668|ref|NC_011475.1|' 'Pseudogi|224381671|ref|NC_011478.1|' 'Pseudogi|224381672|ref|NC_011479.1|' 'Pseudogi|224381674|ref|NC_011481.1|' 'Pseudogi|224381679|ref|NC_011465.1|' 'Pseudogi|224381689|ref|NC_011466.1|' 'Pseudogi|224381692|ref|NC_011469.1|' 'Pseudogi|224381694|ref|NC_011471.1|' > scaffolds.to.anontate.fasta
