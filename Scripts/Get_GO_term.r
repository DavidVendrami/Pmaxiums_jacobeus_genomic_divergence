### Define regions of interest:
# 1. Diagnostic SNPs
# 2. Highly differentiated windows
#	a) windows with Fst > 0.95 quantile
#	b) windows with Fst > 0.99 quantile
# 3. Highly similar chromosomes (mean Fst significanly lower than average)

# 1. Find their locations in: /prj/mar-in-gen/Pecten_M/Enrichment_new/Diagnostic_loci/Diag_loc_rep.txt

# 2.
# Generate a "map" file:
# grep -v '#' Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf | cut -d '   ' -f1,2,3,4,5 > All_loci_report.txt

map<-read.table("/prj/mar-in-gen/Pecten_M/Clean/ref_map_output/All_loci_report.txt",h=F)
colnames(map)<-c("CHR","POS","ID","REF","ALT")

data<-read.table("Window_Fst_MAF05.windowed.weir.fst",h=T) # Read in stepping window fst file
data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside

q95<-quantile(data$WEIGHTED_FST,0.95) # Define q95 threshold
q99<-quantile(data$WEIGHTED_FST,0.99) # Define q99 threshold

w95<-data[which(data$WEIGHTED_FST>=q95),] # Extract relevant windows
w99<-data[which(data$WEIGHTED_FST>=q99),] # Extract relevant windows

w95_out<-data.frame()
for (i in 1:dim(w95)[1]){
w95_out_t<-map[map$CHR==w95[i,1] & map$POS>=w95[i,2] &  map$POS<w95[i,3],]
w95_out<-rbind(w95_out,w95_out_t)
}

w99_out<-data.frame()
for (i in 1:dim(w99)[1]){
w99_out_t<-map[map$CHR==w99[i,1] & map$POS>=w99[i,2] &  map$POS<w99[i,3],]
w99_out<-rbind(w99_out,w99_out_t)
}

write.table(w95_out,"/prj/mar-in-gen/Pecten_M/Enrichment_new/Highly_diff/q95/q95_loc_rep.txt",quote=F,col.names=T,row.names=F,sep="\t")
write.table(w99_out,"/prj/mar-in-gen/Pecten_M/Enrichment_new/Highly_diff/q99/q99_loc_rep.txt",quote=F,col.names=T,row.names=F,sep="\t")

# a) Find their locations in: "/prj/mar-in-gen/Pecten_M/Enrichment_new/Highly_diff/q95/q95_loc_rep.txt"
# b) Find their locations in: "/prj/mar-in-gen/Pecten_M/Enrichment_new/Highly_diff/q99/q99_loc_rep.txt"

# 3.
# Generate stats to define chromosomes that are highly similar
t_test<-data.frame(CHR=c(),p_value=c())
w_test<-data.frame(CHR=c(),p_value=c())

for (i in 1:19){
t_test[i,1]<-paste("HiC_scaffold_",i,sep="")
t_test[i,2]<-t.test(data$WEIGHTED_FST,data$WEIGHTED_FST[data$CHROM==paste("HiC_scaffold_",i,sep="")],alternative="greater")$p.value
w_test[i,1]<-paste("HiC_scaffold_",i,sep="")
w_test[i,2]<-wilcox.test(data$WEIGHTED_FST,data$WEIGHTED_FST[data$CHROM==paste("HiC_scaffold_",i,sep="")],alternative="greater")$p.value
}
w_test[,3]<-p.adjust(w_test$V2,"bonferroni")
w_test[,4]<-p.adjust(w_test$V2,"BH")
t_test[,3]<-p.adjust(t_test$V2,"bonferroni")
t_test[,4]<-p.adjust(t_test$V2,"BH")
colnames(t_test)<-c("CHR","p_value","Bonferroni","FDR")
colnames(w_test)<-c("CHR","p_value","Bonferroni","FDR")

# Let's do this for the moment: let's only look at CHR 4, 5, 6, 8 and 16
# because they have no SNP whose Fst is > q95 and because they are on average significantly undifferentiated as
# compared to genome-wide Fst

ld<-data[data$CHROM=="HiC_scaffold_4" | data$CHROM=="HiC_scaffold_5" | data$CHROM=="HiC_scaffold_6" | data$CHROM=="HiC_scaffold_8" | data$CHROM=="HiC_scaffold_16",]

ld_out<-data.frame()
for (i in 1:dim(ld)[1]){
ld_out_t<-map[map$CHR==ld[i,1] & map$POS>=ld[i,2] &  map$POS<ld[i,3],]
ld_out<-rbind(ld_out,ld_out_t)
}
write.table(ld_out,"/prj/mar-in-gen/Pecten_M/Enrichment_new/Low_diff/ld_loc_rep.txt",quote=F,col.names=T,row.names=F,sep="\t")
# Find their locations in: "/prj/mar-in-gen/Pecten_M/Enrichment_new/Low_diff/ld_loc_rep.txt"


###########################################

### Extract genes (where SNPs are located in or the closest one).
### Do the following for all regions of interest (change input for 'snps' accordingly).
snps<-read.table("/prj/mar-in-gen/Pecten_M/Enrichment_new/Diagnostic_loci/Diag_loc_rep.txt",h=T) # Read in file with SNP locations
data<-read.table("/prj/mar-in-gen/Pecten_M/Enrichment_new/Annotation_files/Pecten_maximus_mRNA.gff",h=F) # Read in annotation file.
colnames(data)<-c("CHR","Start","Stop","Annotation")
data[,1]<-as.character(data[,1])
snps[,1]<-as.character(snps[,1])

out<-data.frame(CHR=c(),Start=c(),Stop=c(),Within_gene=c(),Annotation=c())

for (i in 1:length(snps[,1])){
  
  dat<-data[data$CHR==snps[i,1],]
  ind<-which(dat$Start<=snps[i,2] & dat$Stop>=snps[i,2])
  
  if (length(ind)!=0){
    if (length(ind)>1){
      ind<-ind[1]}
    bla<-dat[ind,]
    bla$Within_gene<-"YES"
    out<-rbind(out,bla)
  } else {
    a<-dat[dat$Start>snps[i,2],]
    b<-dat[dat$Start<=snps[i,2],]
    am<-min(a$Start)
    bm<-max(b$Stop)
    
    if(am<bm){
      ind<-which(dat$Start==am)
      bla<-dat[ind,]
      bla$Within_gene<-"NO"
      out<-rbind(out,bla)
    } else {
      ind<-which(dat$Stop==bm)
      bla<-dat[ind,]
      bla$Within_gene<-"NO"
      out<-rbind(out,bla)
    }
  }
}

rem<-which(duplicated(out$Annotation))
out_r<-out[-rem,]

yes<-out[out$Within_gene=="YES",]
rem_y<-which(duplicated(yes$Annotation))
yes_r<-yes[-rem_y,]

write.table(out_r,"SNP_all_rem.gff",quote=F,col.names=T,row.names=F,sep="\t") # only for GO without duplicates
#write.table(out,"SNP_all_Annotations.gff",quote=F,col.names=T,row.names=F,sep="\t") # to be fully annotated and for GO enrich.anal. to see if having multiple entries of the same GO term change things
write.table(yes_r,"SNP_YES_rem.gff",quote=F,col.names=T,row.names=F,sep="\t") # only for GO without duplicates
#write.table(yes,"SNP_YES_Annotations.gff",quote=F,col.names=T,row.names=F,sep="\t") # only for GO with duplicates

### Switch to BASH and 
# 1. extract gene names for GO enrichment analysis
#cut -d$'\t' -f4 SNP_all_rem.gff | cut -d ';' -f1 | sed 's/ID=//g' | sed '1d' > AllRem_Gene_IDs_for_GOenrich.txt
cut -d$'\t' -f4 SNP_all_Annotations.gff | cut -d ';' -f1 | sed 's/ID=//g' | sed '1d' > All_Gene_IDs_for_GOenrich.txt
#cut -d$'\t' -f4 SNP_YES_rem.gff | cut -d ';' -f1 | sed 's/ID=//g' | sed '1d' > YesRem_Gene_IDs_for_GOenrich.txt
cut -d$'\t' -f4 SNP_YES_Annotations.gff | cut -d ';' -f1 | sed 's/ID=//g' | sed '1d' > Yes_Gene_IDs_for_GOenrich.txt

# 2. Annotate genes with NR, KEGG, SWISSPROT and GO
sed 's/$/\t/g' AllRem_Gene_IDs_for_GOenrich.txt > patterns.txt
grep -f patterns.txt /prj/mar-in-gen/Pecten_M/Enrichment_new/Annotation_files/NR.result | cut -d $'\t' -f1,5,14 | tr ' ' '_' > DB.txt
grep -f patterns.txt /prj/mar-in-gen/Pecten_M/Enrichment_new/Annotation_files/Swissprot.result | cut -d $'\t' -f1,5,16 | tr ' ' '_' > DB_swiss.txt
grep -f patterns.txt /prj/mar-in-gen/Pecten_M/Enrichment_new/Annotation_files/KEGG.result | cut -d $'\t' -f1,5,16 | tr ' ' '_' > DB_kegg.txt


R
data<-read.table("SNP_all_Annotations.gff",h=T)
db<-read.table("DB.txt",h=F)
dbs<-read.table("DB_swiss.txt",h=F)
dbk<-read.table("DB_kegg.txt",h=F)
newd<-read.table("All_Gene_IDs_for_GOenrich.txt",h=F)
data$Annotation<-newd[,1]

for (i in 1:length(data[,1])){
ind<-which(as.character(db[,1])==as.character(data$Annotation[i]))
if (length(ind)!=0){
data$NR_ID[i]<-as.character(db[ind,2])
data$NR_Description[i]<-as.character(db[ind,3])
} else {
data$NR_ID[i]<-NA
data$NR_Description[i]<-NA
}
}

for (i in 1:length(data[,1])){
ind<-which(as.character(dbs[,1])==as.character(data$Annotation[i]))
if (length(ind)!=0){
data$SwissP_ID[i]<-as.character(dbs[ind,2])
data$SwissP_Description[i]<-as.character(dbs[ind,3])
} else {
data$SwissP_ID[i]<-NA
data$SwissP_Description[i]<-NA
}
}

for (i in 1:length(data[,1])){
ind<-which(as.character(dbk[,1])==as.character(data$Annotation[i]))
if (length(ind)!=0){
data$KEGG_ID[i]<-as.character(dbk[ind,2])
data$KEGG_Description[i]<-as.character(dbk[ind,3])
} else {
data$KEGG_ID[i]<-NA
data$KEGG_Description[i]<-NA
}
}

write.table(data,"All_Diag_SNPs_Annotations.txt",quote=F,col.names=T,row.names=F,sep="\t")

# 3. Extract gene sequences, in case you want to annotate them
# bash
cut -d$'\t' -f1,2,3 SNP_all_rem.gff > Genes.bed
bedtools getfasta -fi /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa -bed Genes.bed -fo ./Gene_seqs.fasta

#################################################

# Ready to run enrichment analysis




