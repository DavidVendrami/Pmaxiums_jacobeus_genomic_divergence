### Get % of mapped reads per sample

for i in *_sorted.bam
do
samtools view $i | cut -d$'\t' -f2,5 > bam_stats/${i%_1_sorted.bam}.txt
done

cd bam_stats/

ls -1 *.txt > files.txt

R

files<-read.table("files.txt",h=F)

Sample<-c()
Tot_reads<-c()
Mapped_reads<-c()
HQ_mapped_reads<-c()

for (i in 1:length(files[,1])){
data<-read.table(paste(files[i,1],sep=""),h=F)

Sample[i]<-gsub(".txt","",files[i,1])
Tot_reads[i]<-dim(data)[1]
Mapped_reads[i]<-dim(data[data[,1]!=4,])[1]/dim(data)[1]
HQ_mapped_reads[i]<-dim(data[data[,1]!=4 & data[,2]>=10,])[1]/dim(data)[1]
}

out<-cbind(Sample,Tot_reads,Mapped_reads,HQ_mapped_reads)
colnames(out)<-c("Sample","Tot_reads","Mapped_reads","HQ_mapped_reads")
write.table(out,"Bam_stats.txt",quote=F,col.names=T,row.names=F,sep="\t")
q()
n


