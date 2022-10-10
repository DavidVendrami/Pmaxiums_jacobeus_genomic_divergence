# Divide SNPs by chromosome
grep $'HiC_scaffold_1\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr1.txt
grep $'HiC_scaffold_1\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr1.txt
                                                                             
grep $'HiC_scaffold_2\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr2.txt
grep $'HiC_scaffold_2\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr2.txt
                                                                              
grep $'HiC_scaffold_3\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr3.txt
grep $'HiC_scaffold_3\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr3.txt
                                                                              
grep $'HiC_scaffold_7\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr7.txt
grep $'HiC_scaffold_7\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr7.txt
                                                                             
grep $'HiC_scaffold_9\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr9.txt
grep $'HiC_scaffold_9\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr9.txt
	 
grep $'HiC_scaffold_10\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr10.txt
grep $'HiC_scaffold_10\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr10.txt
                                                                               
grep $'HiC_scaffold_11\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr11.txt
grep $'HiC_scaffold_11\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr11.txt
                                                                               
grep $'HiC_scaffold_12\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr12.txt
grep $'HiC_scaffold_12\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr12.txt
                                                                               
grep $'HiC_scaffold_13\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr13.txt
grep $'HiC_scaffold_13\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr13.txt
                                                                              
grep $'HiC_scaffold_14\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr14.txt
grep $'HiC_scaffold_14\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr14.txt
                                                                              
grep $'HiC_scaffold_15\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr15.txt
grep $'HiC_scaffold_15\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr15.txt
                                                                               
grep $'HiC_scaffold_17\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr17.txt
grep $'HiC_scaffold_17\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr17.txt
                                                                               
grep $'HiC_scaffold_18\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr18.txt
grep $'HiC_scaffold_18\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr18.txt
                                                                               
grep $'HiC_scaffold_19\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr19.txt
grep $'HiC_scaffold_19\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr19.txt

grep $'HiC_scaffold_4\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr4.txt
grep $'HiC_scaffold_4\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr4.txt
                                                                              
grep $'HiC_scaffold_5\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr5.txt
grep $'HiC_scaffold_5\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr5.txt
                                                                              
grep $'HiC_scaffold_6\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr6.txt
grep $'HiC_scaffold_6\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr6.txt
                                                                              
grep $'HiC_scaffold_8\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr8.txt
grep $'HiC_scaffold_8\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr8.txt
	 
grep $'HiC_scaffold_16\t' Final_onlyPmax.map | cut -d$'\t' -f2 > LD_stuff/Pmax_chr16.txt
grep $'HiC_scaffold_16\t' Final_onlyPjac.map | cut -d$'\t' -f2 > LD_stuff/Pjac_chr16.txt
 
################
# Make bed files
for i in {1..19}
do
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile ../Final_onlyPmax --aec --extract Pmax_chr$i.txt --make-bed --out Pmax_chr$i
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile ../Final_onlyPjac --aec --extract Pjac_chr$i.txt --make-bed --out Pjac_chr$i
done

# Calculate pairwise LD among all loci within the same chromosome
for i in {1..19}
do
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile Pmax_chr$i --aec --r2 square --out Pmax_$i
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile Pjac_chr$i --aec --r2 square --out Pjac_$i
done

#########
# Plot LD maps

R

max<-read.table("Pmax_1.ld",h=F)
jac<-read.table("Pjac_1.ld",h=F)

pdf("Pmax_Pjac_chr1.pdf",width=5,height=5)
LDheatmap(as.matrix(max),add.map=F,title=substitute(paste(italic('P. maximus'), " - Chromosome 1")))
LDheatmap(as.matrix(jac),add.map=F,title=substitute(paste(italic('P. jacobeus'), " - Chromosome 1")),newpage=T)
dev.off()

tiff("Pmax_chr1tiff",width=5,height=5, units= 'in', res=600, pointsize=1/600)
LDheatmap(as.matrix(max),add.map=F,title=substitute(paste(italic('P. maximus'), " - Chromosome 1")))
dev.off()

tiff("Pjac_chr1.tiff",width=5,height=5, units= 'in', res=600, pointsize=1/600)
LDheatmap(as.matrix(jac),add.map=F,title=substitute(paste(italic('P. jacobeus'), " - Chromosome 1")))
dev.off()

# Repeat for all othe chromosomes
