# Perform enrichment analysis on regions of interest

library(topGO)
# library(qvalue)
geneID2GO <- readMappings(file = choose.files()) #GO annotations of all genes
geneUniverse <- names(geneID2GO)
genesOfInterest <- read.table(choose.files(),header=FALSE) #Significant loci
genesOfInterest <- as.character(genesOfInterest[,1])
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdataMF <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO) #Choose "BP","MF" (or "CC")
myGOdataBP <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO) #Choose "BP","MF" (or "CC")
# resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher") # 'classic' does not take into account gene hierarchy: false positives
resultFisherMFw <- runTest(myGOdataMF, algorithm="weight01", statistic="fisher") # 'weight01' takes into account gene hierarchy.
resultFisherMFe <- runTest(myGOdataMF, algorithm="elim", statistic="fisher") # 'elim' starts with the most specific term. Most stringent.
resultFisherBPw <- runTest(myGOdataBP, algorithm="weight01", statistic="fisher") # 'weight01' takes into account gene hierarchy.
resultFisherBPe <- runTest(myGOdataBP, algorithm="elim", statistic="fisher") # 'elim' starts with the most specific term. Most stringent.

# allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 1386) #N is the number of scored GO terms (you can see it when you type resultFisher)
allResMFw <- GenTable(myGOdataMF, classicFisher = resultFisherMFw, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(resultFisherMFw@score)) #N is the number of scored GO terms (you can see it when you type resultFisher)
allResMFe <- GenTable(myGOdataMF, classicFisher = resultFisherMFe, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(resultFisherMFe@score)) #N is the number of scored GO terms (you can see it when you type resultFisher)
allResBPw <- GenTable(myGOdataBP, classicFisher = resultFisherBPw, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(resultFisherBPw@score)) #N is the number of scored GO terms (you can see it when you type resultFisher)
allResBPe <- GenTable(myGOdataBP, classicFisher = resultFisherBPe, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = length(resultFisherBPe@score)) #N is the number of scored GO terms (you can see it when you type resultFisher)

write.table(allResMFw,"MF_Diag_Weight.txt",quote=F,col.names=T,row.names=F,sep="\t")
write.table(allResMFe,"MF_Diag_Elim.txt",quote=F,col.names=T,row.names=F,sep="\t")
write.table(allResBPw,"BP_Diag_Weight.txt",quote=F,col.names=T,row.names=F,sep="\t")
write.table(allResBPe,"BP_Diag_Elim.txt",quote=F,col.names=T,row.names=F,sep="\t")



# Print position of significant GOs
# showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 3, useInfo ='all')
# printGraph(myGOdata, resultFisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE) # pdf




