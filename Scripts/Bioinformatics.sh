# 1. Reads QC with process radtags
/vol/animalbehaviour/davidlee/bin/stacks-2.60/process_radtags -p /prj/mar-in-gen/Pecten_M/originals -o /prj/mar-in-gen/Pecten_M/Clean -e pstI -c -q -r
# Run with:
qsub -P fair_share -cwd -l idle=1 -e /prj/mar-in-gen/Pecten_M/scripts -o /prj/mar-in-gen/Pecten_M/scripts -l vf=4G process_radtags.sh

# 2. Now we can align the reads to the reference genome
# First, index the reference
bwa index /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa
# Then align. We will have to use BWA aln/samse because the reads are really short. Pipe the output to samtools to convert to bam and sort.
for i in /prj/mar-in-gen/Pecten_M/Clean/*.fq.gz
do
gzip -d $i
/vol/biotools/bin/bwa aln -t 8 /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa ${i%.gz} > ${i%.fq.gz}.sai
/vol/biotools/bin/ bwa samse /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa ${i%.fq.gz}.sai ${i%.gz} | \
/vol/biotools/bin/samtools view -b -h | /vol/biotools/bin/samtools sort -o ${i%.fq.gz}_sorted.bam
gzip ${i%.gz}
done
# Run with:
qsub -P fair_share -cwd -l idle=1 -e /prj/mar-in-gen/Pecten_M/scripts -o /prj/mar-in-gen/Pecten_M/scripts -pe multislot 8 -l vf=8G Align_and_sort.sh


# 3. Now we can use stacks ref_map.pl to call SNPs.
/grp/animalbehaviour/davidlee/bin/bin/ref_map.pl --samples /prj/mar-in-gen/Pecten_M/Clean/ --popmap /prj/mar-in-gen/Pecten_M/Clean/popmap.txt -o ref_map_output -T 12

qsub -P fair_share -cwd -l idle=1 -e /prj/mar-in-gen/Pecten_M/scripts -o /prj/mar-in-gen/Pecten_M/scripts -pe multislot 12 -l vf=8G Ref_map.sh

# 4. Filtering
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 5 --minGQ 5 --recode --out Pec_Biall_DP5
# 323 x 1,051,046

# Let's output missingness per individual to see if we clearly have to remove someone.
vcftools --vcf Pec_Biall_DP5.recode.vcf --missing-indv --out Pec_Biall_DP5_MissInd

# Let's remove the individuals with more than 800000 missing genotypes, see "Pec_miss_data_perInd.pdf"
vcftools --vcf Pec_Biall_DP5.recode.vcf --remove Inds_to_remove.txt --recode --out Pec_Biall_DP5_GoodSam
# 284 x 1,051,046

# Let's explor the effect of --max-missing
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.9 --out Pec_Biall_DP5_MissInd_MD90 # 129,260
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.8 --out Pec_Biall_DP5_MissInd_MD80 # 337,210
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.7 --out Pec_Biall_DP5_MissInd_MD70 # 466,933

# Let's try both 80 and 90
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.8 --recode --out Pec_Biall_DP5_MissInd_MD80
# 284 x 337,210
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.9 --recode --out Pec_Biall_DP5_MissInd_MD90
# 284 x 129.260

# Too high cov 
vcftools --vcf Pec_Biall_DP5_MissInd_MD80.recode.vcf --site-mean-depth --out Pec_Biall_DP5_MissInd_MD80_depth
# Let's filter for this later if needed, no bimodal distrib.

# Let's make plink files
/vol/animalbehaviour/davidlee/bin/plink/plink --vcf Pec_Biall_DP5_MissInd_MD80.recode.vcf --out Pec_Biall_DP5_MissInd_MD80 --double-id --aec
/vol/animalbehaviour/davidlee/bin/plink/plink --vcf Pec_Biall_DP5_MissInd_MD90.recode.vcf --out Pec_Biall_DP5_MissInd_MD90 --double-id --aec

# Remove replicates,
vcftools --vcf Pec_Biall_DP5_MissInd_MD80.recode.vcf --missing-indv --out Pec_Biall_DP5_MissInd_MD80_SNPs
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD80 --aec --remove Remove_inds.txt --make-bed --out Pec_Biall_DP5_MissInd_MD80_NoRep
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD90 --aec --remove Remove_inds.txt --make-bed --out Pec_Biall_DP5_MissInd_MD90_NoRep
# 262 samples

# add FID in R

# make ped/map and raw
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD80_NoRep --aec --recode --out Pec_Biall_DP5_MissInd_MD80_NoRep
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD90_NoRep --aec --recode --out Pec_Biall_DP5_MissInd_MD90_NoRep
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD80_NoRep --aec --recodeA --out Pec_Biall_DP5_MissInd_MD80_NoRep
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD90_NoRep --aec --recodeA --out Pec_Biall_DP5_MissInd_MD90_NoRep
### 262 x (337,210 - MD80 | 129,260 - MD90)
### 262 = 27 Pj + 235 Pm




# make ped/map, raw also after excluding BC pops (found P. jacobeus there?, PEB-PEB, RDB-CAM, TND-T, PLY-C) 
# and artificail pops (MUL-A, BAN-R)
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD80_NoRep --aec --remove Remove_BC_and_artificial.txt --make-bed --out Pec_Biall_DP5_MissInd_MD80_NoRep_NoBC
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD90_NoRep --aec --remove Remove_BC_and_artificial.txt --make-bed --out Pec_Biall_DP5_MissInd_MD90_NoRep_NoBC
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD80_NoRep_NoBC --aec --recode --out Pec_Biall_DP5_MissInd_MD80_NoRep_NoBC
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD90_NoRep_NoBC --aec --recode --out Pec_Biall_DP5_MissInd_MD90_NoRep_NoBC
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD80_NoRep_NoBC --aec --recodeA --out Pec_Biall_DP5_MissInd_MD80_NoRep_NoBC
/vol/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_MissInd_MD90_NoRep_NoBC --aec --recodeA --out Pec_Biall_DP5_MissInd_MD90_NoRep_NoBC
### 150 x (337,210 - MD80 | 129,260 - MD90)
### 150 = 27 Pj + 123 Pm

