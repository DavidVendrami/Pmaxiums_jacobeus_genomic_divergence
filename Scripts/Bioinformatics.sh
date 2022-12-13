##########################################################
### This files provide the commands to process the raw ###
### reads all the way to SNP genotyping and fitlering. ###
##########################################################

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
/vol/biotools/bin/bwa samse /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa ${i%.fq.gz}.sai ${i%.gz} | \
/vol/biotools/bin/samtools view -b -h | /vol/biotools/bin/samtools sort -o ${i%.fq.gz}_sorted.bam
gzip ${i%.gz}
done
# Run with:
qsub -P fair_share -cwd -l idle=1 -e /prj/mar-in-gen/Pecten_M/scripts -o /prj/mar-in-gen/Pecten_M/scripts -pe multislot 8 -l vf=8G Align_and_sort.sh


# 3. Now we can use stacks ref_map.pl to call SNPs.
/grp/animalbehaviour/davidlee/bin/bin/ref_map.pl --samples /prj/mar-in-gen/Pecten_M/Clean/ --popmap /prj/mar-in-gen/Pecten_M/Clean/popmap.txt -o /prj/mar-in-gen/Pecten_M/Clean/ref_map_output -T 12 -X "populations: --vcf"

qsub -P fair_share -cwd -l idle=1 -e /prj/mar-in-gen/Pecten_M/scripts -o /prj/mar-in-gen/Pecten_M/scripts -pe multislot 12 -l vf=8G Ref_map.sh

# 4. Filtering
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 5 --minGQ 5 --recode --out Pec_Biall_DP5
# 280 x 1,009,368

# Let's output missingness per individual to see if we clearly have to remove someone.
vcftools --vcf Pec_Biall_DP5.recode.vcf --missing-indv --out Pec_Biall_DP5_MissInd

# Let's remove the individuals with more than 800000 missing genotypes, see "Pec_miss_data_perInd.pdf"
vcftools --vcf Pec_Biall_DP5.recode.vcf --remove Inds_to_remove.txt --recode --out Pec_Biall_DP5_GoodSam
# 280 x 1,009,368

# Let's explor the effect of --max-missing
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.9 --out Pec_Biall_DP5_MissInd_MD90 # 129,260
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.8 --out Pec_Biall_DP5_MissInd_MD80 # 337,210
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.7 --out Pec_Biall_DP5_MissInd_MD70 # 466,933

# Let's try both 80 and 90
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.8 --recode --out Pec_Biall_DP5_MissInd_MD80
# 280 x 321,973
vcftools --vcf Pec_Biall_DP5_GoodSam.recode.vcf --max-missing 0.9 --recode --out Pec_Biall_DP5_MissInd_MD90
# 280 x 129,260

# Too high cov 
vcftools --vcf Pec_Biall_DP5_MissInd_MD80.recode.vcf --site-mean-depth --out Pec_Biall_DP5_MissInd_MD80_depth
# Let's filter for this later if needed, no bimodal distrib.

# Let's make plink files
/vol/animalbehaviour/davidlee/bin/plink/plink --vcf Pec_Biall_DP5_MissInd_MD80.recode.vcf --out Pec_Biall_DP5_MissInd_MD80 --double-id --aec
/vol/animalbehaviour/davidlee/bin/plink/plink --vcf Pec_Biall_DP5_MissInd_MD90.recode.vcf --out Pec_Biall_DP5_MissInd_MD90 --double-id --aec

# Remove samples with too many missing genotypes.
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
### 262 x (321,973 - MD80 | 129,260 - MD90)
### 262 = 27 Pj + 235 Pm

###############################################################
### The reviewer asked to check whether different fitlering ###
### thresholds affect downstream results.                   ###
###############################################################
# 1. Let's see if there are SNPs genotyped only in P. jacobaeus (starting from a not --max-missing and --maf filtered dataset).

/grp/animalbehaviour/davidlee/bin/plink/plink --vcf Pec_Biall_DP5_GoodSam.recode.vcf --out Pec_Biall_DP5_GoodSam --double-id --aec # Convert to plink

R
data<-read.table("Pec_Biall_DP5_GoodSam.fam",h=F)
fid<-c()
fid[1:23]<-"Pmax"
fid[24:31]<-"Pjac"
fid[32:118]<-"Pmax"
fid[119:139]<-"Pjac"
fid[140:284]<-"Pmax"
data[,1]<-fid
write.table(data,"Pec_Biall_DP5_GoodSam.fam",quote=F,col.names=F,row.names=F)

/grp/animalbehaviour/davidlee/bin/plink/plink --bfile Pec_Biall_DP5_GoodSam --recodeA --aec --out Pec_Biall_DP5_GoodSam

R
library(data.table)
data<-fread("Pec_Biall_DP5_GoodSam.raw",h=T)
# Split dataset by species
pm_ind<-which(data$FID=="Pmax")
pmd<-data[pm_ind,]
pjd<-data[-pm_ind,]
pmd<-pmd[,-c(1:6)] # 255 samples
pjd<-pjd[,-c(1:6)] # 29 samples

# Total number of SNPs: 1009368

# Count missing genotypes per locus for the two species
pm_md<-apply(pmd,2, function(x) length(which(is.na(x))))
pj_md<-apply(pjd,2, function(x) length(which(is.na(x))))

# Amount of loci with only missing genotypes
length(which(pm_md==255))
# 5368: 0.5%
length(which(pj_md==29))
# 42519: 4.2%

# Loci with genotypes only in Pjac:
pmmd_ind<-which(pm_md==255)
pj_md_u<-pj_md[pmmd_ind]
length(which(pj_md_u!=29))
# 2520, 0.24%
# Of which with at least 80% genotypes:
length(which(pj_md_u<5.8))
# 33, 1.3% of pj_md_u, 0.003% of total SNPs and 0.07% of final SNPs

# Quite useless. But list the names anyways and do the FST thing.

ind<-which(pj_md_u<5.8)
pju<-names(pj_md_u[ind])
write.table(pju,"SNPs_uniquely_genotyped_in_Pj.txt",quote=F,col.names=F,row.names=F,sep="\t")

sed -i 's/_[ACTG]//g' SNPs_uniquely_genotyped_in_Pj.txt

# 2. No MAF filtering: already available

# 3. MAF 1%
vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep.recode.vcf --maf 0.01 --recode --out ../For_rev_comments/MAF/MAF_01

# 4. HWE fitlerings. Calculate HWE stats separately for each pop and then exclude SNPs out of HWE in (i) half and (ii) all pops.
# Subset file by pop to calculate HWE separately
sed -n '1,20p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Mulroy_sub.txt
sed -n '21,28p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Adriatic_sub.txt
sed -n '29,49p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Oban_sub.txt
sed -n '50,69p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Plymouth_sub.txt
sed -n '70,87p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > CAM_sub.txt
sed -n '88,107p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Toralla_sub.txt
sed -n '108,126p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Gdl_sub.txt
sed -n '127,146p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Shetland_sub.txt
sed -n '147,166p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Norway1_sub.txt
sed -n '167,187p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Norway2_sub.txt
sed -n '188,208p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > NPA_sub.txt
sed -n '209,226p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > PEB_sub.txt
sed -n '227,246p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Bantry_sub.txt
sed -n '247,262p' ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.fam | cut -d' ' -f2 > Tinduff_sub.txt

for i in *_sub.txt
do
vcftools --vcf ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --keep $i --recode --out ${i%_sub.txt}
/grp/animalbehaviour/davidlee/bin/plink/plink --vcf ${i%_sub.txt}.recode.vcf --out ${i%_sub.txt} --double-id --aec # Let's do a plink file before MD filtering
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile ${i%_sub.txt} --aec --hardy --out ${i%_sub.txt}_HWE
done

ls -1 *HWE.hwe > file_list.txt

R
files<-read.table("file_list.txt",h=F)
snps01<-c()
snps05<-c()
for (i in 1:length(files[,1])){
data<-read.table(paste(files[i,1]),h=T)
ind05<-which(data$P<0.05)
ind01<-which(data$P<0.01)
snps05<-c(snps05,as.character(data$SNP[ind05]))
snps01<-c(snps01,as.character(data$SNP[ind01]))
}
table(table(snps05))
#     1     2     3     4     5     6     7     8     9    10    11    12    13	14
# 12312  3364   723   162    62    58    28    25    26    38    30    33    54	44

table(table(snps01))
#    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 6102  446   75   55   40   28   15   13   23   24   14   28   51 
hwe05_all_ind<-names(which(table(snps05)==14))
hwe05_half_ind<-names(which(table(snps05)>=7))
hwe01_all_ind<-names(which(table(snps01)>=7))

write.table(hwe05_all_ind,"hwe05_all.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(hwe05_half_ind,"hwe05_half.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(hwe01_all_ind,"hwe01_all.txt",quote=F,col.names=F,row.names=F,sep="\t")

grep -f hwe05_all.txt ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf | cut -d$'\t' -f1,2 > hwe05_all.pos
grep -f hwe05_half.txt ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf | cut -d$'\t' -f1,2 > hwe05_half.pos
grep -f hwe01_all.txt ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf | cut -d$'\t' -f1,2 > hwe01_all.pos

vcftools --vcf ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --exclude-positions hwe05_all.pos --recode --out HWE_fitlered_files/hwe05_all
vcftools --vcf ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --exclude-positions hwe05_half.pos --recode --out HWE_fitlered_files/hwe05_half
vcftools --vcf ../../ref_map_output/Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --exclude-positions hwe01_all.pos --recode --out HWE_fitlered_files/hwe01_all

