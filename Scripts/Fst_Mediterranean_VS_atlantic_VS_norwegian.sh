# Let's get windowed fst stats for med-atl-nor

# get IDs
cp P_jacobeus.txt Med_Atl_Nor/
grep -f Norvegian_IDs.txt P_maximus.txt > Med_Atl_Nor/P_norvegicus.txt # XD
grep -f Norvegian_IDs.txt -v P_maximus.txt > Med_Atl_Nor/P_maximus.txt # these are the atlantic ones

vcftools --vcf ../Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --weir-fst-pop P_norvegicus.txt --weir-fst-pop P_maximus.txt --fst-window-size 100000 --out Window_Fst_MaxNor
vcftools --vcf ../Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --weir-fst-pop P_jacobeus.txt --weir-fst-pop P_maximus.txt --fst-window-size 100000 --out Window_Fst_MaxJac
vcftools --vcf ../Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --weir-fst-pop P_jacobeus.txt --weir-fst-pop P_norvegicus.txt --fst-window-size 100000 --out Window_Fst_JacNor

# NorMax
# Weir and Cockerham mean Fst estimate: 0.03455
# Weir and Cockerham weighted Fst estimate: 0.039194
# 
# MaxJac
# Weir and Cockerham mean Fst estimate: 0.085529
# Weir and Cockerham weighted Fst estimate: 0.17089
# 
# JacNor
# Weir and Cockerham mean Fst estimate: 0.10728
# Weir and Cockerham weighted Fst estimate: 0.20085

pdf("/homes/davidlee/Desktop/WinFst_NorMax_2.pdf",width=15,height=10)
#layout(matrix(c(1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8), 4, 4, byrow = TRUE))
par(mfrow=c(2,1))

# MaxJac
data<-read.table("../Window_Fst_MAF05.windowed.weir.fst",h=T) # Read in stepping window fst file
data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside
data<-data[-c(4458:length(data$CHROM)),] # Remove unplaced scaffolds (only for the figure!)

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
#col[4458:length(col)]<-"grey85" 
data$col<-col

labs<-table(factor(data$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
#bla[20]<-4457

plot(1:length(data[,1]),data$WEIGHTED_FST,xaxt='n',yaxt='n',col=transp(data$col,.6),main="",ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Fst",cex.lab=1.25, cex.axis=1.25)
axis(1,at=bla,labels=c(seq(1,19,by=1)),cex.lab=0.4, cex.axis=1.25)
axis(2,at=c(0,0.5,1),cex.lab=0.4, cex.axis=1.25,las=2)
abline(h=c(quantile(data$WEIGHTED_FST,probs=c(0.95))),lty=2,col=c("#de2d26"))

#d <- density(data$WEIGHTED_FST)
#plot(d,main="",xlab="Fst",xlim=c(-0.1,1),yaxt='n',ylab="")
#polygon(d, col=transp("#7570b3",.6), border="#7570b3")

# Plot MaxatlJac
# data<-read.table("Window_Fst_MaxJac.windowed.weir.fst",h=T) # Read in stepping window fst file
# data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside
# 
# col<-c()
# cols<-c("grey25","grey85")
# col[1]<-"grey25"
# for (i in 2:length(data$CHROM)){
#   if(data$CHROM[i]==data$CHROM[i-1]){
#     col[i]<-col[i-1]
#   } else {
#     col[i]<-cols[which(cols!=col[i-1])]
#   }
# }
# col[4457:length(col)]<-"grey85" 
# data$col<-col
# 
# labs<-table(factor(data$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
# bla<-c()
# bla[1]<-labs[1]/2
# for (i in 2:19){
# bla[i]<-sum(labs[1:i-1])+labs[i]/2
# }
# bla[20]<-bla[19]+labs[19]/2+17
# 
# # layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
# 
# plot(1:length(data[,1]),data$WEIGHTED_FST,xaxt='n',col=transp(data$col,.6),main="P. maximus (atlantic population) vs P. jacobeus, N = 4459",ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Stepping window Fst",cex.lab=1.25, cex.axis=1.25)
# axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
# abline(h=c(quantile(data$WEIGHTED_FST,probs=c(0.95))),lty=2,col=c("#de2d26"))
# 
# d <- density(data$WEIGHTED_FST)
# plot(d,main="",xlab="Fst",xlim=c(-0.1,1),yaxt='n',ylab="")
# polygon(d, col=transp("#7570b3",.6), border="#7570b3")

# Plot JacNor
# data<-read.table("Window_Fst_JacNor.windowed.weir.fst",h=T) # Read in stepping window fst file
# data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside
# 
# col<-c()
# cols<-c("grey25","grey85")
# col[1]<-"grey25"
# for (i in 2:length(data$CHROM)){
#   if(data$CHROM[i]==data$CHROM[i-1]){
#     col[i]<-col[i-1]
#   } else {
#     col[i]<-cols[which(cols!=col[i-1])]
#   }
# }
# col[4439:length(col)]<-"grey85" 
# data$col<-col
# 
# labs<-table(factor(data$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
# bla<-c()
# bla[1]<-labs[1]/2
# for (i in 2:19){
# bla[i]<-sum(labs[1:i-1])+labs[i]/2
# }
# bla[20]<-bla[19]+labs[19]/2+17

# layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))

# plot(1:length(data[,1]),data$WEIGHTED_FST,xaxt='n',col=transp(data$col,.6),main="P. maximus (norvegian population) vs P. jacobeus, N = 4441",ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Stepping window Fst",cex.lab=1.25, cex.axis=1.25)
# axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
# abline(h=c(quantile(data$WEIGHTED_FST,probs=c(0.95))),lty=2,col=c("#de2d26"))
# 
# d <- density(data$WEIGHTED_FST)
# plot(d,main="",xlab="Fst",xlim=c(-0.1,1),yaxt='n',ylab="")
# polygon(d, col=transp("#7570b3",.6), border="#7570b3")

# Plot NorMax
data<-read.table("Window_Fst_MaxNor.windowed.weir.fst",h=T) # Read in stepping window fst file
data<-data[data$N_VARIANTS>=5,] # Retain only windows with at least 10 SNPs inside
data<-data[-c(4446:length(data$CHROM)),] # Remove unplaced scaffolds (only for the figure!)

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
#col[4446:length(col)]<-"grey85" 
data$col<-col

labs<-table(factor(data$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
#bla[20]<-bla[19]+labs[19]/2+17


# layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))

plot(1:length(data[,1]),data$WEIGHTED_FST,xaxt='n',yaxt='n',col=transp(data$col,.6),main="",ylim=c(0,1),pch=16,xlab="Chromosomes",ylab="Fst",cex.lab=1.25, cex.axis=1.25)
axis(1,at=bla,labels=c(seq(1,19,by=1)),cex.lab=0.4, cex.axis=1.25)
axis(2,at=c(0,0.5,1),cex.lab=0.4, cex.axis=1.25,las=2)
abline(h=c(quantile(data$WEIGHTED_FST,probs=c(0.95))),lty=2,col=c("#de2d26"))

#d <- density(data$WEIGHTED_FST)
#plot(d,main="",xlab="Fst",xlim=c(-0.1,1),yaxt='n',ylab="")
#polygon(d, col=transp("#7570b3",.6), border="#7570b3")

dev.off()




