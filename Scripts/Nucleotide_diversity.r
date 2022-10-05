### Calculate Pi for all samples and for the 2 species separately

vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --keep P_jacobeus.txt --recode --out Final_Pjca # 40317 SNPs
vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --remove P_jacobeus.txt --recode --out Final_Pmax # 45540 SNPs

vcftools --vcf Pec_Biall_DP5_MissInd_MD80_NoRep_MAF05.recode.vcf --window-pi 100000 --out Pi_All # 45704 SNPs
vcftools --vcf Final_Pmax.recode.vcf --window-pi 100000 --out Pi_Pmax
vcftools --vcf Final_Pjca.recode.vcf --window-pi 100000 --out Pi_Pjac

### R
### All
# Repeat for p. maximus and p. jacobeus
pi<-read.table("Pi_Pjac.windowed.pi",h=T)
pi<-pi[pi$N_VARIANTS>=5,]
labs<-table(factor(pi$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
    bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
bla[20]<-bla[19]+labs[19]/2+17

col<-c()
cols<-c("grey25","grey85")
col[1]<-"grey25"
for (i in 2:length(pi$CHROM)){
    if(pi$CHROM[i]==pi$CHROM[i-1]){
        col[i]<-col[i-1]
    } else {
        col[i]<-cols[which(cols!=col[i-1])]
    }
}
#col[4446:length(col)]<-"grey85" #All
#col[4446:length(col)]<-"grey85" #Pmax
col[3879:length(col)]<-"grey85" #Pjac
pi$col<-col

# Define outliers in red
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
col[4442:length(col)]<-"grey85" 
data$col<-col
q95<-quantile(data$WEIGHTED_FST,prob=(0.95))
data$col[data$WEIGHTED_FST>=q95]<-"red"

t1<-paste(data$CHROM,data$BIN_START,sep="_")
t2<-paste(pi$CHROM,pi$BIN_START,sep="_")
indd<-which(t1%in%t2)
indp<-which(t2%in%t1)

#for pmax
datam<-data[indd,]
pim<-pi[indp,]
ind<-which(datam$col=="red")
pim$col[ind]<-"red"

#for pjac
dataj<-data[indd,]
pij<-pi[indp,]
ind<-which(dataj$col=="red")
pij$col[ind]<-"red"

# Plot pi while highlighting highly divergent windows
pdf("Pi_fst_red.pdf",width=10,height=8)
par(mfrow=c(2,1))
par(xpd=T)
#plot(1:length(pia[,1]),pia$PI,xlim=c(0,length(pia[,1])),xaxt='n',col=transp(pia$col,.6),ylim=c(0,6e-05),pch=16,xlab="Chromosomes",ylab=expression(pi),cex.lab=1.25, cex.axis=1.25)
#pia$x<-1:length(pia[,1])
#par(new=T)
#piar<-pia[pia$col=="red",]
#plot(piar$x,piar$PI,xlim=c(0,length(pia[,1])),xaxt='n',col=transp(piar$col,.6),ylim=c(0,6e-05),pch=16,xlab="Chromosomes",ylab=expression(pi),cex.lab=1.25, cex.axis=1.25)
#labs<-table(factor(pia$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
#bla<-c()
#bla[1]<-labs[1]/2
#for (i in 2:19){
#bla[i]<-sum(labs[1:i-1])+labs[i]/2
#}
#bla[20]<-bla[19]+labs[19]/2+17
#axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
#text(-160,8e-5,"(a) Both species",cex=1.4)

plot(1:length(pim[,1]),pim$PI,xlim=c(0,length(pim[,1])),xaxt='n',yaxt='n',col=transp(pim$col,.6),ylim=c(0,6e-05),pch=16,xlab="Chromosomes",ylab=expression(pi),cex.lab=1.25, cex.axis=1.25)
pim$x<-1:length(pim[,1])
par(new=T)
pimr<-pim[pim$col=="red",]
plot(pimr$x,pimr$PI,xlim=c(0,length(pim[,1])),xaxt='n',yaxt='n',col=transp(pimr$col,.6),ylim=c(0,6e-05),pch=16,xlab="",ylab="",cex.lab=1.25, cex.axis=1.25)
labs<-table(factor(pim$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
bla[20]<-bla[19]+labs[19]/2+17
axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
axis(2,at=c(0,0.00003,0.00006),labels=c("0","3e-05","6e-05"))
text(-170,8e-5,substitute(paste('(a) ',italic('P. maximus'))),cex=1.4)

plot(1:length(pij[,1]),pij$PI,xlim=c(0,length(pij[,1])),xaxt='n',yaxt='n',col=transp(pij$col,.6),ylim=c(0,6e-05),pch=16,xlab="Chromosomes",ylab=expression(pi),cex.lab=1.25, cex.axis=1.25)
pij$x<-1:length(pij[,1])
par(new=T)
pijr<-pij[pij$col=="red",]
plot(pijr$x,pijr$PI,xlim=c(0,length(pij[,1])),xaxt='n',yaxt='n',col=transp(pijr$col,.6),ylim=c(0,6e-05),pch=16,xlab="",ylab="",cex.lab=1.25, cex.axis=1.25)
labs<-table(factor(pij$CHROM,levels=c("HiC_scaffold_1","HiC_scaffold_2","HiC_scaffold_3","HiC_scaffold_4","HiC_scaffold_5","HiC_scaffold_6","HiC_scaffold_7","HiC_scaffold_8","HiC_scaffold_9","HiC_scaffold_10","HiC_scaffold_11","HiC_scaffold_12","HiC_scaffold_13","HiC_scaffold_14","HiC_scaffold_15","HiC_scaffold_16","HiC_scaffold_17","HiC_scaffold_18","HiC_scaffold_19")))[1:19]
bla<-c()
bla[1]<-labs[1]/2
for (i in 2:19){
bla[i]<-sum(labs[1:i-1])+labs[i]/2
}
bla[20]<-bla[19]+labs[19]/2+17
axis(1,at=bla,labels=c(seq(1,19,by=1),"U"),cex.lab=0.4, cex.axis=1.25)
axis(2,at=c(0,0.00003,0.00006),labels=c("0","3e-05","6e-05"))
text(-160,8e-5,substitute(paste('(b) ',italic('P. jacobeus'))),cex=1.4)
dev.off()

# Alternative way of visualizing relationship between Fst and pi
pdf("/homes/davidlee/Desktop/Fst_pi_spec.pdf",width=8,height=6)
datam$col[datam$col=="grey85"]<-"grey25"
par(fig=c(0,0.5,0,0.8))
plot(pim$PI,datam$WEIGHTED_FST,xlab=expression(pi ~ (x ~ 10^-5)),ylab=expression(italic(F)[ST]),col=transp(datam$col,.2),main="",xaxt='n',pch=16)
axis(1,at=c(0,0.00001,0.00002,0.00003,0.00004,0.00005),labels=c("0","1","2","3","4","5"))

par(fig=c(0,0.5,0.5,1), new=TRUE)
hist(pim$PI,breaks=100,border=NA,col="grey85",xaxt='n',yaxt='n',xlab="",ylab="",main="")
par(xpd = TRUE)
text(-0.000004,150,substitute(paste('(a) ',italic('P. maximus'))),cex=1.2)

par(new=T)
dataj$col[dataj$col=="grey85"]<-"grey25"
par(fig=c(0.5,1,0,0.8),new=T)
plot(pij$PI,dataj$WEIGHTED_FST,xlab=expression(pi ~ (x ~ 10^-5)),ylab=expression(italic(F)[ST]),col=transp(dataj$col,.2),main="",xaxt='n',pch=16)
axis(1,at=c(0,0.00001,0.00002,0.00003,0.00004,0.00005),labels=c("0","1","2","3","4","5"))

par(fig=c(0.5,1,0.5,1), new=TRUE)
hist(pij$PI,breaks=100,border=NA,col="grey85",xaxt='n',yaxt='n',xlab="",ylab="",main="")
par(xpd = TRUE)
text(-0.0000044,150,substitute(paste('(b) ',italic('P. jacobeus'))),cex=1.2)
dev.off()

