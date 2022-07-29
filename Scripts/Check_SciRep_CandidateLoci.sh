### Check Srep candidates ###
# Start from the .map file containing locations of the Srep candidate loci and extract the contigs sequences from the de novo catalog
# build for the Srep paper

cut -d ' ' -f 1 Candidate.map | sed 's/contig_/>/g' | sed 's/$/ /g' > Candidates_fasta_headers.txt
grep -A 1 -f Candidates_fasta_headers.txt /prj/mar-in-gen/Pecten_M/All_denovo_10Adri_344/catalog.fa | grep -v '-' > Candidate_contiga.fa
# 271 sequences because 8 contigs contain 2 SNPs

# Map the contigs to the genome
/vol/biotools/bin/bwa aln -t 8 /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa Candidate_contiga.fa > Candidate_contiga.sai
/vol/biotools/bin/bwa samse /prj/mar-in-gen/Pecten_M/New_Reference/Pecten_maximus.genomic.fa Candidate_contiga.sai Candidate_contiga.fa > Candidate_mapping.sam

sed -i '1,107d' Candidate_mapping.sam

# Keep interesting columns
cut -d$'\t' -f 1,2,3,4,5 Candidate_mapping.sam > Candidate_mapping.txt
# and add column names

# R
# Windows
data<-read.table("Window_Fst_MAF05.windowed.weir.fst",h=T)
quantile(data$WEIGHTED_FST,prob=0.95) # Based only on Windows with at least 5 SNPs
#      95% 
#0.6264734

cand<-read.table("Candidate_mapping.txt",h=T)
table(cand[,2])
# 34 did not map
cand<-cand[cand$SAM!=4,]

CHROM<-as.character(data$CHROM)
data$CHROM<-as.character(CHROM)

Chr<-as.character(cand$Chr)
cand$Chr<-as.character(Chr)

fst<-c()
map<-c()
N_V<-c()
for (i in 1:dim(cand)[1]){
subs<-data[data$CHROM==cand$Chr[i],]
if (length(which(subs$BIN_START<=cand$Pos[i] & subs$BIN_END>=cand$Pos[i]))!=0){
map[i]<-which(subs$BIN_START<=cand$Pos[i] & subs$BIN_END>=cand$Pos[i])
fst[i]<-subs$WEIGHTED_FST[which(subs$BIN_START<=cand$Pos[i] & subs$BIN_END>=cand$Pos[i])]
N_V[i]<-subs$N_VARIANTS[which(subs$BIN_START<=cand$Pos[i] & subs$BIN_END>=cand$Pos[i])]
} else {
map[i]<-NA
fst[i]<-NA
N_V[i]<-NA
}
}

cand$Window_Fst<-fst
cand$N_V<-N_V

length(which(is.na(cand$Window_Fst)))
# the re are 14 candidate for which there is no window
# these probably are candidate that mapped to multiple locations and we did not called any SNP from there because the 
# corresponding reads also map to multiple locations probably

cand<-cand[!is.na(cand$N_V),]

length(which(cand$N_V<5))
# 60 SNPs are located within windows with < 5 SNPs

length(which(cand$Window_Fst>=0.62))
# 13 are located within >95% windows

Win_95<-c() # is the SNP located within a >95 percentile window?
for (i in 1:length(cand$Chr)){
if (cand$Window_Fst[i]>=0.62){
Win_95[i]<-"Y"
} else {
Win_95[i]<-"N"
}
}

cand$Win_95<-Win_95
write.table(cand,"CandidateSrep_SNPs_Windows_report.txt",quote=F,col.names=T,row.names=F,sep="\t")

# These candidate loci are only associated with O2 and T! 

### Single Fst
# Locating exactly the SNP may not be doable because of differences (indels for example) between genome and denovo assembly
# Let's then find the closest reference-based SNP to each denovo-based candidate SNP and use those Fst values
data<-read.table("SNP_Fst_MAF05.weir.fst",h=T)
cand<-read.table("Candidate_mapping.txt",h=T)
cand<-cand[cand$SAM!=4,]
CHROM<-as.character(data$CHROM)
data$CHROM<-as.character(CHROM)
Chr<-as.character(cand$Chr)
cand$Chr<-as.character(Chr)

dif<-c()
fst<-c()
for (i in 1:dim(cand)[1]){
subs<-data[data$CHROM==cand$Chr[i],]
diff<-abs(cand$Pos[i]-subs$POS)
dif[i]<-min(diff)
fst[i]<-subs$WEIR_AND_COCKERHAM_FST[which(diff==min(diff))]
}

cand$dif<-dif
cand$fst<-fst

sig<-c()
for (i in 1:length(cand$Chr)){
if (cand$fst[i]>=0.65){
sig[i]<-"Y"
} else {
sig[i]<-"N"
}
}

# 14 candidates with high fst
# among them 6 are new (queries: 10709, 16316, 60808, 77928, 126116, 159363).
write.table(cand,"CandidateSrep_SNPs_SingleFst_report.txt",quote=F,col.names=T,row.names=F,sep="\t")

