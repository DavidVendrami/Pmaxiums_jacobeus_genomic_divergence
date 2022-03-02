# 71 Diagnostic loci

data<-read.table("Final_All_MAF05.raw",h=T)
geno<-data[,-c(1:6)]
genom<-geno[data$FID=="P_maximus",]
genoj<-geno[data$FID=="P_jacobeus",]
lm<-apply(genom,2,function(x) dim(table(x)))
lj<-apply(genoj,2,function(x) dim(table(x)))
jna<-colnames(genoj[,which(lj==1)])
mna<-colnames(genom[,which(lm==1)])
length(which(jna%in%mna))

Diag1<-(genoj[,jna[which(jna%in%mna)]])
Diag2<-(genom[,mna[which(mna%in%jna)]])
Fin<-rbind(Diag1,Diag2)
ind<-apply(Fin,2,function(x) dim(table(x)))
length(which(ind==2))
[1] 71
indd<-which(ind==2)

mne<-which(colnames(data)%in%names(indd))
mne<-mne-6
map<-read.table("Final_All_MAF05.map",h=F)
Floc<-map[mne,]
View(Floc)
write.table(Floc,"Diagnostic_loci.txt",quote=F,col.names=F,row.names=F,sep="\t")
