#Reduced stress defense responses contribute to the higher toxicity of a pesticide under warming. 
#Delnat V., Swaegers J., Asselman A., and Stoks R. (2020). 
#Molecular Ecology, accepted.
#R code tested on 01/06/2020

######Packages#######
#install packages
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("tximport")     
BiocManager::install("jsonlite")     
BiocManager::install("readr")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
install.packages("car")     
install.packages("lme4")     
install.packages("lsmeans")     
install.packages("effects")     
install.packages("afex")     

#load packages
library(tximport)
library(jsonlite)
library(readr)
library(DESeq2)
library(ggplot2)
library(car)     
library(lme4)
library(emmeans)
library(effects)     
library(afex)

#package to export graphs to svg or ppt
install.packages("digest")
library(digest)
install.packages("devtools")
library(devtools)
devtools::install_github("tomwenseleers/export", force=TRUE)
library(export)


######Read in count files#######
#read in count files using tximport and a table with sample names and treatments
sampleTable=read.table("sampleTableText.txt", sep="\t", header=TRUE, na.strings=c(""))
sampleTable$SampleName=as.factor(sampleTable$SampleName)
sampleTable$Temperature=as.factor(sampleTable$Temperature)
sampleTable$Temperature=factor(sampleTable$Temperature, levels=c("ambient","warming"))
sampleTable$Chlorpyrifos=as.factor(sampleTable$Chlorpyrifos)
sampleTable$Chlorpyrifos=factor(sampleTable$Chlorpyrifos, levels=c("control","low","high"))
#use count files received from Salmon (with --validateMappings and --gcBias flags to improve quality)
files <- file.path("./CountFiles", sampleTable$SampleName,"quant.sf")
txi <- tximport(files, type = c("salmon"), txOut = TRUE)

#make the dds (dataset of DESeq) using txi
dataDESeq = DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ Temperature * Chlorpyrifos)
dataDESeq$group <- factor(paste0(dataDESeq$Temperature, ".", dataDESeq$Chlorpyrifos))
design(dataDESeq) <- ~ group
dds <- DESeq(dataDESeq)

#remove low counts
keep <- rowSums(counts(dds)) >= 10
dds10 <- dds[keep,]
#Reference dataset of 62527 transcripts after filtering (count(dds) >=10)


######DESeq2#######
#Contrast 1: 24°C-solvent control vs 20°C-solvent control (single effect of warming)
DESeqWarmingSolventControl=results(dds10, alpha=0.05, lfcThreshold=1, contrast=c("group", "warming.control", "ambient.control"))
DESeqWarmingSolventControlALL = as.data.frame(DESeqWarmingSolventControl[order(DESeqWarmingSolventControl$log2FoldChange),])
write.table(DESeqWarmingSolventControlALL,file="DESeqWarmingSolventControlALL.txt", sep='\t', quote=F)
summary(DESeqWarmingSolventControl)
DESeqDEGsWarmingSolventControl=subset(DESeqWarmingSolventControl,padj < 0.05 & abs(log2FoldChange) > 1) 
DESeqDEGsWarmingSolventControl = as.data.frame(DESeqDEGsWarmingSolventControl[order(DESeqDEGsWarmingSolventControl$log2FoldChange),])
write.table(DESeqDEGsWarmingSolventControl,file="DESeqDEGsWarmingSolventControl.txt", sep='\t', quote=F)

#Contrast 2: 20°C-low-effect chlorpyrifos vs 20°C-solvent control (single effect of low-effect chlorpyrifos)
DESeqAmbientLowChlorpyrifos=results(dds10, alpha=0.05, lfcThreshold=1, contrast=c("group", "ambient.low", "ambient.control"))
DESeqAmbientLowChlorpyrifosALL = as.data.frame(DESeqAmbientLowChlorpyrifos[order(DESeqAmbientLowChlorpyrifos$log2FoldChange),])
write.table(DESeqAmbientLowChlorpyrifosALL,file="DESeqAmbientLowChlorpyrifosALL.txt", sep='\t', quote=F)
DESeqAmbientLowChlorpyrifos=results(dds10, alpha=0.05, contrast=c("group", "ambient.low", "ambient.control"))
summary(DESeqAmbientLowChlorpyrifos)
DESeqDEGsAmbientLowChlorpyrifos=subset(DESeqAmbientLowChlorpyrifos,padj < 0.05 & abs(log2FoldChange) > 1) 
DESeqDEGsAmbientLowChlorpyrifos = as.data.frame(DESeqDEGsAmbientLowChlorpyrifos[order(DESeqDEGsAmbientLowChlorpyrifos$log2FoldChange),])
write.table(DESeqDEGsAmbientLowChlorpyrifos,file="DESeqDEGsAmbientLowChlorpyrifos.txt", sep='\t', quote=F)

#Contrast 3: 20°C-high-effect chlorpyrifos vs 20°C-solvent control (single effect of high-effect chlorpyrifos)
DESeqAmbientHighChlorpyrifos=results(dds10, alpha=0.05, lfcThreshold=1, contrast=c("group", "ambient.high", "ambient.control"))
DESeqAmbientHighChlorpyrifosALL = as.data.frame(DESeqAmbientHighChlorpyrifos[order(DESeqAmbientHighChlorpyrifos$log2FoldChange),])
write.table(DESeqAmbientHighChlorpyrifosALL,file="DESeqAmbientHighChlorpyrifosALL.txt", sep='\t', quote=F)
DESeqAmbientHighChlorpyrifos=results(dds10, alpha=0.05, contrast=c("group", "ambient.high", "ambient.control"))
summary(DESeqAmbientHighChlorpyrifos)
DESeqDEGsAmbientHighChlorpyrifos=subset(DESeqAmbientHighChlorpyrifos,padj < 0.05 & abs(log2FoldChange) > 1) 
DESeqDEGsAmbientHighChlorpyrifos = as.data.frame(DESeqDEGsAmbientHighChlorpyrifos[order(DESeqDEGsAmbientHighChlorpyrifos$log2FoldChange),])
write.table(DESeqDEGsAmbientHighChlorpyrifos,file="DESeqDEGsAmbientHighChlorpyrifos.txt", sep='\t', quote=F)

#Contrast 4: 24°C-low-effect chlorpyrifos vs 20°C-solvent control (combined effect of low-effect chlorpyrifos and warming)
DESeqWarmingLowChlorpyrifos=results(dds10, alpha=0.05, lfcThreshold=1, contrast=c("group", "warming.low", "ambient.control"))
DESeqWarmingLowChlorpyrifosALL = as.data.frame(DESeqWarmingLowChlorpyrifos[order(DESeqWarmingLowChlorpyrifos$log2FoldChange),])
write.table(DESeqWarmingLowChlorpyrifosALL,file="DESeqWarmingLowChlorpyrifosALL.txt", sep='\t', quote=F)
summary(DESeqWarmingLowChlorpyrifos)
DESeqDEGsWarmingLowChlorpyrifos=subset(DESeqWarmingLowChlorpyrifos,padj < 0.05 & abs(log2FoldChange) > 1) 
DESeqDEGsWarmingLowChlorpyrifos = as.data.frame(DESeqDEGsWarmingLowChlorpyrifos[order(DESeqDEGsWarmingLowChlorpyrifos$log2FoldChange),])
write.table(DESeqDEGsWarmingLowChlorpyrifos,file="DESeqDEGsWarmingLowChlorpyrifos.txt", sep='\t', quote=F)

#Contrast 5: 24°C-high-effect chlorpyrifos vs 20°C-solvent control (combined effect of high-effect chlorpyrifos and warming)
DESeqWarmingHighChlorpyrifos=results(dds10, alpha=0.05, lfcThreshold=1, contrast=c("group", "warming.high", "ambient.control"))
DESeqWarmingHighChlorpyrifosALL = as.data.frame(DESeqWarmingHighChlorpyrifos[order(DESeqWarmingHighChlorpyrifos$log2FoldChange),])
write.table(DESeqWarmingHighChlorpyrifosALL,file="DESeqWarmingHighChlorpyrifosALL.txt", sep='\t', quote=F)
summary(DESeqWarmingHighChlorpyrifos)
DESeqDEGsWarmingHighChlorpyrifos=subset(DESeqWarmingHighChlorpyrifos,padj < 0.05 & abs(log2FoldChange) > 1) 
DESeqDEGsWarmingHighChlorpyrifos = as.data.frame(DESeqDEGsWarmingHighChlorpyrifos[order(DESeqDEGsWarmingHighChlorpyrifos$log2FoldChange),])
write.table(DESeqDEGsWarmingHighChlorpyrifos,file="DESeqDEGsWarmingHighChlorpyrifos.txt", sep='\t', quote=F)


######Venn diagram#######

AmbientLowChlorpyrifos=row.names(DESeqDEGsAmbientLowChlorpyrifos)
AmbientHighChlorpyrifos=row.names(DESeqDEGsAmbientHighChlorpyrifos)
WarmingSolventControl=row.names(DESeqDEGsWarmingSolventControl)

install.packages("VennDiagram")
library(VennDiagram)

venn.diagram(x=list(AmbientLowChlorpyrifos, AmbientHighChlorpyrifos, WarmingSolventControl), 
             category.names=c("(a) 20°C-low-effect chlorpyrifos", "(b) 20°C-high-effect chlorpyrifos", "(c) warming (24°C-solvent control)"),
             cat.default.pos = "outer", cat.pos = c(-18, 18, 170), cex= 4, cat.cex = 2, fill = c("grey53", "grey29","grey11"),
             resolution = 300, filename = 'Venn_diagram_General_Stress_Response.png', output=TRUE)

#Not used, only DEGs with known blast were used. Venn diagram was made in venny version 2.1 (Oliveros, 2015)             

######Reaction norm######

install.packages("reshape")     
install.packages("lattice")     
library(reshape)
library(lattice)

#Log2FC values of DEGs present in contrast 1, 2 and 3 (main effects, effect of single stressors) with a fully defined blast annotation
LFC_BlastGeneralStressDefence=read.table("BlastReactionNormMainEffects.txt", sep="\t", header=TRUE)
graph.data.blast<-melt(LFC_BlastGeneralStressDefence, id=c("replicate"))

#Plot graph and save it as svg (further editing can be done in e.g. InkScape)
svg("ReactionNorms_GeneralStressDefenceBlast.svg", width = 9, height = 9)
xyplot(value~  variable, group=as.factor(replicate), 
       cex=1.5, layout=(c(1,1)), pch=16, lty=1, col="black",
       type = c("p","l"), data=graph.data.blast, scales=list(tck=c(1,1),
       xlab =list(label="Treatment",fontsize=13), ylab = list(label=expression(paste("LogFC")),fontsize=13), 
       x=list(cex=1.0, labels=c("24°C Control","20°C Low CPF","20°C High CPF")), y=list(cex=1.0)))
dev.off()

#Log2FC values of all DEGs present in contrast 1, 2 and 3 (main effects, effect of single stressors)
LFC_AllGeneralStressDefence=read.table("AllReactionNormMainEffects.txt", sep="\t", header=TRUE)
graph.data.all<-melt(LFC_AllGeneralStressDefence, id=c("replicate"))

#Plot graph and save it as svg (further editing can be done in e.g. InkScape)
svg("ReactionNorms_GeneralStressDefenceAll.svg", width = 9, height = 9)
xyplot(value~  variable, group=as.factor(replicate), 
      cex=1.5, layout=(c(1,1)), pch=16, lty=1, col="grey",
      type = c("p","l"), data=graph.data.all, scales=list(tck=c(1,1),
      xlab =list(label="Treatment",fontsize=13), ylab = list(label=expression(paste("LogFC")),fontsize=13), 
      x=list(cex=1.0, labels=c("24°C Control","20°C Low CPF","20°C High CPF")), y=list(cex=1.0)))
dev.off()


######MA plot#######
#Supplementary G
#ylab="log2 Fold Change", xlab="Mean Expression"

##Main effects
#MA plots showing the mean expression against the log2 fold changes of the differentially expressed genes (DEGs) against 20°C-solvent control  
#for (A) 24°C-solvent control (green, warming), 20°C-low-effect chlorpyrifos (blue) and 20°C-high-effect chlorpyrifos (red).
plot(DESeqWarmingSolventControl$log2FoldChange ~ DESeqWarmingSolventControl$baseMean,log='x', cex.axis=2.0, 
     ylab="", xlab="", pch=16, cex=0.5, col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif")
abline(h=0, col='black',lwd=1.5)
points(DESeqDEGsAmbientHighChlorpyrifos$log2FoldChange ~ DESeqDEGsAmbientHighChlorpyrifos$baseMean, col = "red2", pch = 18, cex=1.2)
points(DESeqDEGsAmbientLowChlorpyrifos$log2FoldChange ~ DESeqDEGsAmbientLowChlorpyrifos$baseMean, col = "royalblue2", pch = 17, cex=0.9)
points(DESeqDEGsWarmingSolventControl$log2FoldChange ~ DESeqDEGsWarmingSolventControl$baseMean, col = "springgreen2", pch = 15, cex=0.9)
#Export plot to svg
graph2svg(file="MAplot_MainEffects", width=8, height=10)

##Interaction Temp x Low CPF
#MA plots showing the mean expression against the log2 fold changes of the differentially expressed genes (DEGs) against 20°C-solvent control  
#for (B) 24°C-solvent control (green, warming), 20°C-low-effect chlorpyrifos (blue) and 24°C-loweffect chlorpyrifos (dark blue).
plot(DESeqWarmingSolventControl$log2FoldChange ~ DESeqWarmingSolventControl$baseMean,log='x', cex.axis=2.0, 
     ylab="", xlab="", pch=16, cex=0.5, col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif")
abline(h=0, col='black',lwd=1.5)
points(DESeqDEGsWarmingLowChlorpyrifos$log2FoldChange ~ DESeqDEGsWarmingLowChlorpyrifos$baseMean, col = "darkblue", pch = 17, cex=0.9)
points(DESeqDEGsAmbientLowChlorpyrifos$log2FoldChange ~ DESeqDEGsAmbientLowChlorpyrifos$baseMean, col = "royalblue2", pch = 18, cex=1.2)
points(DESeqDEGsWarmingSolventControl$log2FoldChange ~ DESeqDEGsWarmingSolventControl$baseMean, col = "springgreen2", pch = 15, cex=0.9)
#Export plot to svg
graph2svg(file="MAplot_InteractionLow", width=8, height=10)

##Interaction Temp x High CPF
#MA plots showing the mean expression against the log2 fold changes of the differentially expressed genes (DEGs) against 20°C-solvent control  
#for (C) 24°C-solvent control (green, warming), 20°C-high-effect chlorpyrifos (red) and 24°C-high-effect chlorpyrifos (dark red). 
plot(DESeqWarmingSolventControl$log2FoldChange ~ DESeqWarmingSolventControl$baseMean,log='x', cex.axis=2.0, 
     ylab="", xlab="", pch=16, cex=0.5, col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif")
abline(h=0, col='black',lwd=1.5)
points(DESeqDEGsWarmingHighChlorpyrifos$log2FoldChange ~ DESeqDEGsWarmingHighChlorpyrifos$baseMean, col = "darkred", pch = 17, cex=0.9)
points(DESeqDEGsAmbientHighChlorpyrifos$log2FoldChange ~ DESeqDEGsAmbientHighChlorpyrifos$baseMean, col = "red2", pch = 18, cex=1.2)
points(DESeqDEGsWarmingSolventControl$log2FoldChange ~ DESeqDEGsWarmingSolventControl$baseMean, col = "springgreen2", pch = 15, cex=0.9)
#Export plot to svg
graph2svg(file="MAplot_InteractionHigh", width=8, height=10)


######Volcano plot#######
#Supplementary H
#xlab = "log2 fold change", ylab = "-log10 padj"

#Contrast 1: 24°C-solvent control vs 20°C-solvent control (single effect of warming)
tab = data.frame(logFC = DESeqWarmingSolventControl$log2FoldChange, negLogPval = -log10(DESeqWarmingSolventControl$padj))
plot(tab, pch = 16, cex = 0.6, cex.axis = 2, ylim=c(0,10), xlab = "", ylab = "",
     col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif", cex.lab=1.3)
xtick<-seq(-20, 20, by=5)
axis(side=1, at=xtick, labels = FALSE,fg="gray35")
#Highlight DEGs
lfc = 1
padj = 0.05
signGenesP = (tab$logFC > lfc & tab$negLogPval > -log10(padj))
signGenesN = (tab$logFC < -lfc & tab$negLogPval > -log10(padj))
points(tab[signGenesP, ], pch = 16, cex = 0.8, col = "springgreen2") #magenta3
points(tab[signGenesN, ], pch = 16, cex = 0.8, col = "red2")   #orange
#Visualize thresholds
abline(h = -log10(padj), col = "blue", lty = 2, lwd=2)  #blue
abline(v = c(-lfc, lfc), col = "blue4", lty = 2, lwd=2) #green
#Export plot to svg
graph2svg(file="VolcanoPlots_WarmingSolventControl", width=8, height=6)

#Contrast 2: 20°C-low-effect chlorpyrifos vs 20°C-solvent control (single effect of low-effect chlorpyrifos)
tab = data.frame(logFC = DESeqAmbientLowChlorpyrifos$log2FoldChange, negLogPval = -log10(DESeqAmbientLowChlorpyrifos$padj))
plot(tab, pch = 16, cex = 0.6, cex.axis = 2, xlab = "", ylab = "",
     col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif", cex.lab=1.3)
xtick<-seq(-30,20, by=5)
axis(side=1, at=xtick, labels = FALSE,fg="gray35")
#Highlight DEGs
lfc = 1
padj = 0.05
signGenesP = (tab$logFC > lfc & tab$negLogPval > -log10(padj))
signGenesN = (tab$logFC < -lfc & tab$negLogPval > -log10(padj))
points(tab[signGenesP, ], pch = 16, cex = 0.8, col = "springgreen2") 
points(tab[signGenesN, ], pch = 16, cex = 0.8, col = "red2")  
#Visualize thresholds
abline(h = -log10(padj), col = "blue", lty = 2, lwd=2) 
abline(v = c(-lfc, lfc), col = "blue4", lty = 2, lwd=2) 
#Export plot to svg
graph2svg(file="VolcanoPlots_AmbientLowChlorpyrifos", width=8, height=6)

#Contrast 3: 20°C-high-effect chlorpyrifos vs 20°C-solvent control (single effect of high-effect chlorpyrifos)
tab = data.frame(logFC = DESeqAmbientHighChlorpyrifos$log2FoldChange, negLogPval = -log10(DESeqAmbientHighChlorpyrifos$padj))
plot(tab, pch = 16, cex = 0.6, cex.axis = 2, xlab = "", ylab = "",
     col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif", cex.lab=1.3)
xtick<-seq(-30,20, by=5)
axis(side=1, at=xtick, labels = FALSE,fg="gray35")
#Highlight DEGs
lfc = 1
padj = 0.05
signGenesP = (tab$logFC > lfc & tab$negLogPval > -log10(padj))
signGenesN = (tab$logFC < -lfc & tab$negLogPval > -log10(padj))
points(tab[signGenesP, ], pch = 16, cex = 0.8, col = "springgreen2") 
points(tab[signGenesN, ], pch = 16, cex = 0.8, col = "red2")  
#Visualize thresholds
abline(h = -log10(padj), col = "blue", lty = 2, lwd=2) 
abline(v = c(-lfc, lfc), col = "blue4", lty = 2, lwd=2) 
#Export plot to svg
graph2svg(file="VolcanoPlots_AmbientHighChlorpyrifos", width=8, height=6)

#Contrast 4: 24°C-low-effect chlorpyrifos vs 20°C-solvent control (combined effect of low-effect chlorpyrifos and warming)
tab = data.frame(logFC = DESeqWarmingLowChlorpyrifos$log2FoldChange, negLogPval = -log10(DESeqWarmingLowChlorpyrifos$padj))
plot(tab, pch = 16, cex = 0.6, xlab = "", ylab = "",
     col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif", cex.lab=1.3)
xtick<-seq(-25, 20, by=5)
axis(side=1, at=xtick, labels = FALSE,fg="gray35")
#Highlight DEGs
lfc = 1
padj = 0.05
signGenesP = (tab$logFC > lfc & tab$negLogPval > -log10(padj))
signGenesN = (tab$logFC < -lfc & tab$negLogPval > -log10(padj))
points(tab[signGenesP, ], pch = 16, cex = 0.8, col = "springgreen2") 
points(tab[signGenesN, ], pch = 16, cex = 0.8, col = "red2") 
#Visualize thresholds
abline(h = -log10(padj), col = "blue", lty = 2, lwd=2) 
abline(v = c(-lfc, lfc), col = "blue4", lty = 2, lwd=2) 
#Export plot to svg
graph2svg(file="VolcanoPlots_WarmingLowChlorpyrifos", width=8, height=6)

#Contrast 5: 24°C-high-effect chlorpyrifos vs 20°C-solvent control (combined effect of high-effect chlorpyrifos and warming)
tab = data.frame(logFC = DESeqWarmingHighChlorpyrifos$log2FoldChange, negLogPval = -log10(DESeqWarmingHighChlorpyrifos$padj))
plot(tab, pch = 16, cex = 0.6, xlab = "", ylab = "",
     col="gray45",col.lab="gray35",col.axis="gray35",fg="gray35",family="serif", cex.lab=1.3)
xtick<-seq(-30, 20, by=5)
axis(side=1, at=xtick, labels = FALSE,fg="gray35")
#Highlight DEGs
lfc = 1
padj = 0.05
signGenesP = (tab$logFC > lfc & tab$negLogPval > -log10(padj))
signGenesN = (tab$logFC < -lfc & tab$negLogPval > -log10(padj))
points(tab[signGenesP, ], pch = 16, cex = 0.8, col = "springgreen2") 
points(tab[signGenesN, ], pch = 16, cex = 0.8, col = "red2")  
#Visualize thresholds
abline(h = -log10(padj), col = "blue", lty = 2, lwd=2) 
abline(v = c(-lfc, lfc), col = "blue4", lty = 2, lwd=2) 
#Export plot to svg
graph2svg(file="VolcanoPlots_WarmingHighChlorpyrifos", width=8, height=6)


######PCA plot#######
#Supplementary L
#All unique differentially expressed genes (DEGs) based on contrast 1-5 of DESeq2 (with unique = each DEG transcript ID occuring only once)
DEGsWarmingSolventControl=rownames(DESeqDEGsWarmingSolventControl)
DEGsAmbientLowChlorpyrifos=rownames(DESeqDEGsAmbientLowChlorpyrifos)
DEGsAmbientHighChlorpyrifos=rownames(DESeqDEGsAmbientHighChlorpyrifos)
DEGsWarmingLowChlorpyrifos=rownames(DESeqDEGsWarmingLowChlorpyrifos)
DEGsWarmingHighChlorpyrifos=rownames(DESeqDEGsWarmingHighChlorpyrifos)
allDEGs<-c(DEGsWarmingSolventControl,DEGsAmbientLowChlorpyrifos,DEGsAmbientHighChlorpyrifos,
           DEGsWarmingLowChlorpyrifos,DEGsWarmingHighChlorpyrifos)
allDEGsUnique=unique(allDEGs)
select=allDEGsUnique

#Apply a variance stabilizing transformation on dds
vsd=vst(dds10)
#Principal component analysis
pcaDataDEG = plotPCA(vsd[select,],intgroup=c("Temperature","Chlorpyrifos"),returnData=TRUE)
percentVar <- round(100 * attr(pcaDataDEG, "percentVar"))
#PCA plot of the 36 samples with temperature and chlorpyrifos treatments indicated
ggplot(pcaDataDEG, aes(x = PC1, y = PC2, color = Temperature, shape = Chlorpyrifos)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
#Export plot to ppt
graph2ppt(file="PCAplot_DEGs", width=9, height=6, append=TRUE)


######Interaction types: Frequency#######
#Numbers extracted from 'SupplementaryK_InteractionType.xlsx' using Pivot Table 

dataAntagonistic=read.table("./InteractionType_Frequency_Antagonistic.txt", sep="\t", header=TRUE)
row.names(dataAntagonistic)=c("Antagonistic","Non-Antagonistic")
Antagonistic=data.matrix(dataAntagonistic, rownames.force = NA)
Antagonistic
barplot(Antagonistic, beside=T, legend=T)
fisher.test(Antagonistic,alternative = "two.sided",  conf.int=T, conf.level=0.95)

dataSynergistic=read.table("./InteractionType_Frequency_Synergistic.txt", sep="\t", header=TRUE)
row.names(dataSynergistic)=c("Synergistic","Non-Synergistic")
Synergistic=data.matrix(dataSynergistic, rownames.force = NA)
Synergistic
barplot(Synergistic, beside=T, legend=T)
fisher.test(Synergistic,alternative = "two.sided",  conf.int=T, conf.level=0.95)

dataAdditive=read.table("./InteractionType_Frequency_Additive.txt", sep="\t", header=TRUE)
row.names(dataAdditive)=c("Additive","Non-Additive")
Additive=data.matrix(dataAdditive, rownames.force = NA)
Additive
barplot(Additive, beside=T, legend=T)
fisher.test(Additive,alternative = "two.sided",  conf.int=T, conf.level=0.95)


######Interaction types: Strength#######

##All DEGs of which the interaction type was determined

#read in the dataset
dataMDR=read.csv("./InteractionTypeStrength/InteractionType_MDR.csv", sep=",", na.strings=c(""))
#MDR values between 0 and 1 indicate the strength of antagonistic interaction effects, 
#while MDR values > 1 indicate the strength of synergistic interaction effects.
#Problem = yes indicates MDR values < 0 due to the nature of gene expression data (19/577 DEGs, ~3%)
#only use the transcripts where the interaction type is certain (see Supplementary B and C - Implications Ymax(down) in IA model)
dataMDRsubset=subset(dataMDR,Problem=="no" & InteractionTypeOK=="yes")
dataMDRsubset$InteractionSimplified=factor(dataMDRsubset$InteractionSimplified, 
                                        levels=c("additive","antagonistic downregulation","antagonistic upregulation",
                                                 "synergistic downregulation","synergistic upregulation"))
str(dataMDRsubset)

#model with interaction terms: effect coding
set_sum_contrasts()

#MDR is log10 transformed to meet the assumptions
#log10MDR < 0 indicate the strength of the antagonistic interaction effects
#log10MDR > 0 indicate the strength of the synergistic interaction effects
dataMDRsubset$MDRtransformed=log10(dataMDRsubset$MDR+0.1) 

MDRall=lm(MDRtransformed ~ CPF*InteractionSimplified, data=dataMDRsubset, na.action=na.omit) 
Anova(MDRall, type="III",white.adjust = TRUE) 
#Posthoc test - Contrasts with fdr correction
pairs(emmeans(MDRall, ~ CPF|InteractionSimplified, adjust="fdr"))
#emmeans and standard errors for figure
emmeans(MDRall, ~ CPF|InteractionSimplified, adjust="fdr")

#Assumption - Normality - OK after transformation
shapiro.test(resid(MDRall)) 
hist(resid(MDRall)) 
#Assumption - Homogeneity of variance - not OK --> white.adjust=TRUE in Anova
leveneTest(MDRtransformed ~ InteractionType, data = dataMDRsubset)
aggregate(MDRtransformed ~ InteractionType, data = dataMDRsubset, var)


##DEG subset based on Sulmon et al. (2015) of which the interaction type was determined

#read in the dataset
dataMDRstress=read.csv("./InteractionTypeStrength/InteractionType_MDRofGeneralStressResponse.csv", sep=",", na.strings=c(""))
#MDR values between 0 and 1 indicate the strength of antagonistic interaction effects, 
#while MDR values > 1 indicate the strength of synergistic interaction effects.
#Problem = yes indicates MDR values < 0 due to the nature of gene expression data (19/577 DEGs, ~3%)
#only use the transcripts where the interaction type is certain (see Supplementary B and C - Implications Ymax(down) in IA model)
dataMDRstressSubset=subset(dataMDRstress,Problem=="no" & InteractionTypeOK=="yes")
dataMDRstressSubset$InteractionSimplified=factor(dataMDRstressSubset$InteractionSimplified, 
                                           levels=c("additive","antagonistic downregulation","antagonistic upregulation",
                                                    "synergistic downregulation","synergistic upregulation"))
#Not enough DEGs to statistically test the MDR of synergistic interaction effects
dataMDRstressSubsetAddAnt=subset(dataMDRstressSubset,InteractionType!="Synergistic")
str(dataMDRstressSubsetAddAnt)

#MDR is log10 transformed to meet the assumptions
#log10MDR < 0 indicate the strength of the antagonistic interaction effects
#log10MDR > 0 indicate the strength of the synergistic interaction effects
dataMDRstressSubsetAddAnt$MDRtransformed=log10(dataMDRstressSubsetAddAnt$MDR+0.1) 

MDRstress=lm(MDRtransformed ~ CPF*InteractionSimplified, data=dataMDRstressSubsetAddAnt, na.action=na.omit) 
Anova(MDRstress, type="III",white.adjust = TRUE) 
#Posthoc test - Contrasts with fdr correction
pairs(emmeans(MDRstress, ~ CPF|InteractionSimplified, adjust="fdr"))
#emmeans and standard errors for figure
emmeans(MDRstress, ~ CPF|InteractionSimplified, adjust="fdr")

#Assumption - Normality - OK after transformation
shapiro.test(resid(MDRstress)) 
hist(resid(MDRstress)) 
#Assumption - Homogeneity of variance - not OK --> white.adjust=TRUE in Anova
leveneTest(MDRtransformed ~ InteractionType, data = dataMDRstressSubset)
aggregate(MDRtransformed ~ InteractionType, data = dataMDRstressSubset, var)


######Save Rdata######
save.image(file="Delnat et al_RNAseq-stress-responses_20200601.Rdata")
