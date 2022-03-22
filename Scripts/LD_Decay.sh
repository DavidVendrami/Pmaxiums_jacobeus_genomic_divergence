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
