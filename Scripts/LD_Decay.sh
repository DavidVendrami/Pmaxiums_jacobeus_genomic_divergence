plink --file FILE --maf 0.05 --out FILE_MAF05 --aec --make-bed
plink --bfile FILE_MAF05 --r2 --out FILE_LD_measurments --aec --ld-window-r2 0 --ld-window-kb 1000 --ld-window 99999


# Now process it in Rstudio:
library(data.table)
ld_t<-fread("FILE_LD_measurments.ld",header=T)
ld_t$dist<-(ld_t$BP_B-ld_t$BP_A)/1000
r2<-ld_t$R2
n <- length(unique(c(ld_t$CHR_A[ld_t$dist<10], ld_t$CHR_B[ld_t$dist<10])))
LD.data<-ld_t$R2[ld_t$dist<10]
distance<-ld_t$dist[ld_t$dist<10]
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
ld.df<-data.frame(distance,fpoints)
ld.df<-ld.df[order(ld.df$distance),]
png("LD_Decay_MAF05_final.tiff",width=6, height=7, units= 'in', res=600)
plot(distance,LD.data,pch=19,cex=0.9,xlim=c(0,4),ylim=c(0,0.5),col=transp("black",.2))
lines(ld.df$distance,ld.df$fpoints,lwd=2,col="blue")
dev.off()

#################

# Calculate LD only for highly divergent SNPs and plot it together with previous estimates based on all loci

# Extract highly differentiated SNPs
# To do so, use the positions.txt file used for Enrichment analysis

pmax<-read.table("Final_onlyPmax.map",h=F)
bla<-read.table("positions.txt",h=T)
bla$check<-paste(bla[,1],bla[,2],sep="_")
pmax$check<-paste(pmax[,1],pmax[,4],sep="_")
ind<-which(pmax$check%in%bla$check) #Retain all SNPs located within highly divergent windows
maxout<-pmax[ind,]
outm<-maxout$V2
write.table(outm,"Pmax_HD.txt",quote=F,row.names=F,col.names=F,sep="/t") #This can be used also for P.jacobeus as these are the very same loci (obviously).

plink --bfile Final_onlyPmax --aec --extract Pmax_HD.txt --make-bed --out Pmax_HD
plink --bfile Final_onlyPjac --aec --extract Pmax_HD.txt --make-bed --out Pjac_HD

plink --bfile Pmax_HD --maf 0.05 --out Pmax_HD --aec --make-bed
plink --bfile Pjac_HD --maf 0.05 --out Pjac_HD --aec --make-bed

plink --bfile Pmax_HD --r2 --out Pmax_HD_LD_measurments --aec --ld-window-r2 0 --ld-window-kb 1000 --ld-window 99999
plink --bfile Pjac_HD --r2 --out Pjac_HD_LD_measurments --aec --ld-window-r2 0 --ld-window-kb 1000 --ld-window 99999


library(data.table)
ld_t<-fread("Pmax_LD_measurments.ld",header=T)
ld_t$dist<-(ld_t$BP_B-ld_t$BP_A)/1000
r2<-ld_t$R2
n <- length(unique(c(ld_t$CHR_A[ld_t$dist<100], ld_t$CHR_B[ld_t$dist<100])))
LD.data<-ld_t$R2[ld_t$dist<100]
distance<-ld_t$dist[ld_t$dist<100]
HW.st<-c(C=0.1)
HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
tt<-summary(HW.nonlinear)
new.rho<-tt$parameters[1]
fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
ld.df<-data.frame(distance,fpoints)
ld.df<-ld.df[order(ld.df$distance),]

# Repeat using the other input files and storing the relevant variables as below:
# P. maximus all loci (file: Pmax_LD_measurments.ld):
pmd<-distance
pmLD<-LD.data
pmld.df<-ld.df

# P. jacobeus all loci (file: Pjac_LD_measurments.ld):
pjd<-distance
pjLD<-LD.data
pjld.df<-ld.df

# P. maximus highly divergent loci (file: Pmax_HD_LD_measurments.ld):
pmhdd<-distance
pmhdLD<-LD.data
pmhdld.df<-ld.df

# P. jacobeus highly divergent loci (file: Pjac_HD_LD_measurments.ld):
pjhdd<-distance
pjhdLD<-LD.data
pjhdld.df<-ld.df

tiff("LD_Decay_plots.tiff",width=12, height=10, units= 'in', res=600, pointsize=1/600)
par(mfrow=c(2,2))

plot(pmd,pmLD,ylim=c(0,0.5),xlab="Distance (kb)", ylab="r2", pch=20, col=transp("black",.2),cex=2)
lines(pmld.df$distance,pmld.df$fpoints,lwd=2,col="blue")

plot(pjd,pjLD,ylim=c(0,0.5),xlab="Distance (kb)", ylab="r2", pch=20, col=transp("black",.2),cex=2)
lines(pjld.df$distance,pjld.df$fpoints,lwd=2,col="blue")

plot(pmhdd,pmhdLD,ylim=c(0,0.5),xlab="Distance (kb)", ylab="r2", pch=20, col=transp("black",.2),cex=2)
lines(pmhdld.df$distance,pmhdld.df$fpoints,lwd=2,col="blue")

plot(pjhdd,pjhdLD,ylim=c(0,0.5),xlab="Distance (kb)", ylab="r2", pch=20, col=transp("black",.2),cex=2)
lines(pjhdld.df$distance,pjhdld.df$fpoints,lwd=2,col="blue")

dev.off()
