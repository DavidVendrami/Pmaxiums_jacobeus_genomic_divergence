#####################
# Calculate fst (per-SNP and within sliding window) using vcftools
#

# Make list of Pjac and Pmax sample IDs:
R

data<-read.table("Pec_Biall_DP5_MissInd_MD90_NoRep_NoBC.fam",h=F)
pjac<-as.character(data$V2[data$V1=="P_jacobeus"])
pmax<-as.character(data$V2[data$V1=="P_maximus"])

write.table(pjac,"P_jacobeus.txt",quote=F,col.names=F,row.names=F,sep="\t")
write.table(pmax,"P_maximus.txt",quote=F,col.names=F,row.names=F,sep="\t")

q()
n

### Fst
##### Need to filter for MAF
vcftools ... --maf 0.05
################## All_MD_80_MAF05: 262 x 47944
vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --weir-fst-pop P_jacobeus.txt --weir-fst-pop P_maximus.txt --out SNP_Fst_MAF05
vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --weir-fst-pop P_jacobeus.txt --weir-fst-pop P_maximus.txt --fst-window-size 100000 --out Window_Fst_MAF05
vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --weir-fst-pop P_jacobeus.txt --weir-fst-pop P_maximus.txt --fst-window-size 100000 --fst-window-step 10000 --out Window_Step_Fst_MAF05

# Plot single-SNP FST:
data<-read.table("SNP_Fst_MAF05.weir.fst",h=T) # Read in single-SNP fst file

col<-c()
cols<-c("grey25","grey85")
col[1]<-"grey25"
for (i in 2:length(data$CHROM)){
  if(data$CHROM[i]==data$CHROM[i-1]){
    col[i]<-col[i-1]
  } else {
    col[i]<-cols[which(cols!=col[i-1])]
  }
}
col[47397:length(col)]<-"grey85" 
data$col<-col

labs<-table(data$CHROM)[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
bla[20]<-47750

pdf("SingleSNP_Fst_MD80_MAF05.pdf",width=16,height=6)

plot(1:length(data[,1]),data$WEIR_AND_COCKERHAM_FST,xaxt='n',col=transp(data$col,.6),ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Single SNP Fst",cex.lab=1.25, cex.axis=1.25)
axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
abline(h=c(quantile(data$WEIR_AND_COCKERHAM_FST,probs=c(0.95,0.99,0.999),na.rm=T)),lty=2,col=c("#feb24c","#fc9272","#de2d26"))

dev.off()

############# Plot stepping window
data<-read.table("Window_Fst_MAF05.windowed.weir.fst",h=T) # Read in stepping window fst file
data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside

col<-c()
cols<-c("grey25","grey85")
col[1]<-"grey25"
for (i in 2:length(data$CHROM)){
  if(data$CHROM[i]==data$CHROM[i-1]){
    col[i]<-col[i-1]
  } else {
    col[i]<-cols[which(cols!=col[i-1])]
  }
}
col[4737:length(col)]<-"grey85" 
data$col<-col

labs<-table(data$CHROM)[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
bla[20]<-4750

pdf("Stepping_Windowed_SNP_Fst_MD80_MAF05.pdf",width=16,height=6)

plot(1:length(data[,1]),data$WEIGHTED_FST,xaxt='n',col=transp(data$col,.6),ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Stepping window Fst",cex.lab=1.25, cex.axis=1.25)
axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
abline(h=c(quantile(data$WEIGHTED_FST,probs=c(0.95,0.99,0.999))),lty=2,col=c("#feb24c","#fc9272","#de2d26"))

dev.off()

### Sliding window
data<-read.table("Window_Step_Fst_MAF05.windowed.weir.fst",h=T) # Read in sliding window fst file
data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside

col<-c()
cols<-c("grey25","grey85")
col[1]<-"grey25"
for (i in 2:length(data$CHROM)){
  if(data$CHROM[i]==data$CHROM[i-1]){
    col[i]<-col[i-1]
  } else {
    col[i]<-cols[which(cols!=col[i-1])]
  }
}
col[47240:length(col)]<-"grey85" 
data$col<-col


labs<-table(data$CHROM)[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
bla[20]<-bla[19]+labs[19]/2+17

pdf("Sliding_Windowed_SNP_Fst_MD80_MAF05.pdf",width=16,height=6)

plot(1:length(data[,1]),data$WEIGHTED_FST,xaxt='n',col=transp(data$col,.3),ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Sliding window Fst",cex.lab=1.25, cex.axis=1.25)
axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
abline(h=c(quantile(data$WEIGHTED_FST,probs=c(0.95,0.99,0.999))),lty=2,col=c("#feb24c","#fc9272","#de2d26"))

dev.off()

