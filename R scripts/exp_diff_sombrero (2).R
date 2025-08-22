####matriz datos sombrero
tabla=read.table("/Users/manuel/Documents/posdoc_berkeley/sombrero_proyecto/mat_sombrero.txt", header=TRUE, row.names = 1, sep="\t")
head(tabla)

tabla["col0C.R1_S"]=rowSums(tabla[,c(1,2)])
tabla["col0C.R2_S"]=rowSums(tabla[,c(3,4)])
tabla["col0P.R2_S"]=rowSums(tabla[,c(7,8)])
tabla["smbC.R1_S"]=rowSums(tabla[,c(10,11)])

head(tabla)

colnames(tabla)

tabla2=tabla[,c("col0C.R1_S", "col0C.R2_S", "col0C.R3_1", "col0P.R1_1", "col0P.R2_S", "col0P.R3_1", "smbC.R1_S", "smbC.R2_1", "smbC.R3_1", "smbP.R1_1", "smbP.R2_1", "smbP.R3_1")]
dim(tabla2)


library("edgeR")
library("limma")
library("statmod")
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)


tabla2 = tabla2[rowSums(cpm(tabla2) >= 3) >=5 ,]
class(tabla2)
matTab=as.matrix(tabla2)
row.names(matTab)=row.names(tabla2)
colnames(matTab)=colnames(tabla2)

batch <- c(1,2,3,1,2,3,1,2,3,1,2,3)
class(batch)

adjusted <- ComBat_seq(matTab, batch=batch, group=NULL)
head(tabla2)
head(adjusted)



outpathcount = "/Users/manuel/Documents/posdoc_berkeley/sombrero_proyecto/diff_totales/"
dir.create(outpathcount, showWarnings=FALSE)

counts=adjusted

grp = sub(".....$", "", colnames(counts))
grp
dge = DGEList(counts=counts, group=grp)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]




pdf(paste(outpathcount, "MDS_sombrero", ".pdf"))
plotMDS(dge, top=8000, col=c("black"), xlab="Dim1", ylab="Dim2", dim.plot=c(1,2))
#abline(h= c(1,-1), col="blue", lwd=1, lty=2)
dev.off()


#para calcular el factor de normalizacion.
dge = calcNormFactors(dge)

#dise??o experimental (matriz con todas las cepas) usando a las que tienen el mismo tratamiento como replica.

design = model.matrix(~0 + dge$samples$group)
colnames(design) = levels(dge$samples$group)
design


#calculando la dispersion de los datos.
dge = estimateGLMCommonDisp(dge, design= design)
dge = estimateGLMTrendedDisp(dge, design= design)
dge = estimateGLMTagwiseDisp(dge, design= design)

#plotBCV(dge)

#imprimir dispersion
print(dge$common.dispersion)
################################################################
#fit usando la dispersion calculada
#fit=glmFit(dge, design)
##cambie a esta función para estimate the QL dispersions around the dispersion trend
fit <- glmQLFit(dge, design, robust=TRUE)
plotQLDisp(fit)

# makeContrasts
contVector <- c(
  
  "smbPvssmbC"= "smbP - smbC",
  "col0Pvscol0C"= "col0P - col0C",
  "smbPvscol0P" = "smbP - col0P",
  "smbCvscol0C"= "smbC - col0C"

)


contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)

#differentail expression
contrast_FCFDR=list()
contrast=list()
id = c()
diff_genes=list()
expre_total=list()
list_up=list()
list_down=list()

for (comp in colnames(contMatrix)) { 
   
  print(comp)
  
  tr <- glmTreat(fit, contrast=contMatrix[, comp])
  topTab <- topTags(tr, n=Inf)$table
  
  
  
  deGenes<- rownames(topTab)[topTab$FDR < 0.05 & abs(topTab$logFC) > 1]
  deGenes_up <- rownames(topTab)[topTab$FDR < 0.05 & (topTab$logFC) >1]
  deGenes_down <- rownames(topTab)[topTab$FDR < 0.05 & (topTab$logFC) < -1]
  totaldiff<-length(deGenes)
  up<-length(deGenes_up)
  down<-length(deGenes_down)
  nodif<-(length(row.names(topTab))-totaldiff)
  genes<-c(totaldiff, up, down, nodif)
  names(genes)=c("DEG", "UP", "DOWN", "Non-DEG")
  contrast[[comp]]=genes
  cont=t(as.data.frame(contrast))
  
  list_up[[comp]]=deGenes_up
  list_down[[comp]]=deGenes_down
  print(length(deGenes))
  print(length(deGenes_up))
  print(length(deGenes_down))
  
  ###exporting figure
  pdf(paste(outpathcount, comp, ".pdf", sep=""))
  
  
  plotSmear (tr, de.tags=deGenes, ylab= paste(comp, "logFC"), cex= 0.3, col= "black")
  abline(h= c(1,-1), col="blue", lwd=1, lty=2)
  dev.off()
  
  diff=topTab[as.character(deGenes),]
  #TF=na.exclude(diff[TFgene$Gene_ID,])
  #Con=na.exclude(diff[congenes$genes_con,])
  #inter=na.exclude(topTab[inters$genes_con,])
  
  
  pdf(paste(outpathcount, comp, ".pdf"))
  plot(topTab$logCPM, topTab$logFC, col="darkmagenta", pch = 19, cex=0.3, main=comp, ylim = c(-10,10))
  points(diff$logCPM, diff$logFC, col="burlywood2", pch = 1, cex=1)
  #  points(TF$logCPM, TF$logFC, col="cornflowerblue",pch = 19, cex=1.5)
  # points(Con$logCPM, Con$logFC, col="brown1",pch = 19, cex=1.8 )
  #  text(inter$logCPM, inter$logFC, labels = inters$name,col="black", cex=1.2, pos=3 )
  dev.off()
  
  diff_genes[[comp]]=topTab[deGenes,]
  expre_total[[comp]]=topTab
  id = unique(c(id, deGenes))
  
  write.table(topTab, file=paste(outpathcount, comp, ".txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
  write.table(deGenes_up, file=paste(outpathcount, comp, "_up.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
  write.table(deGenes_down, file=paste(outpathcount, comp, "_down.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
  
}


##
x=list()
ids=list()
for (i in names(diff_genes)){
  print(paste(i, length(diff_genes[[i]]$logFC)))
  x[[i]]= length(diff_genes[[i]]$logFC)
  ids[[i]]=row.names(diff_genes[[i]])
  
}

ids2=as.data.frame(unlist(ids))
ids2=unique(ids2$`unlist(ids)`)
length(ids2)

#tab_coun=counts[ids2,]
###sustituyendo esto con la tabla de promedios de lecturas por condicion

tab_coun=counts[ids2,]
dim(tab_coun)
head(tab_coun)
# install if necessary
#install.packages("dendextend")
# load package
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendextend)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(tab_coun, 1, cal_z_score))
range(data_subset_norm)
my_hclust_gene <- hclust(dist(data_subset_norm), method = "complete")

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
names=read.csv("/Users/manuel/Google Drive/trichoseq/tricho_anot/anot_tricho.csv", row.names = 1)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 12)
clust=as.data.frame(my_gene_col)
l <- split(clust, clust$my_gene_col)
x11()
for(cl in names(l)) {
  tab=data_subset_norm[row.names(l[[cl]]),]
  row.names(tab)= sub('gene-TRIATDRAFT_', '', row.names(tab)) 
  pdf(paste("heatmap_fc3",cl,".pdf"))
  pheatmap(tab, cluster_cols = FALSE, show_rownames = FALSE)  
  dev.off() 
  tab2=tab_coun[row.names(l[[cl]]),]
  row.names(tab2)= sub('gene-TRIATDRAFT_', '', row.names(tab2)) 
  anot=names[row.names(tab),1]
  tab2=cbind(tab2, anot)
  tab2=as.data.frame(tab2)
  write.table(tab2, file=paste(cl, "fc3_cpm.txt", sep=""), sep="\t")
  print(length(l[[cl]]$my_gene_col))
}

names=read.csv("/Users/manuel/Google Drive/trichoseq/tricho_anot/anot_tricho.csv", row.names = 1)
head(names)

###cortando arbol y guardando pdf
pdf("heatmap_cpm_cutree_12_ncbi_FC2.pdf") 
x11()
pheatmap(data_subset_norm, cluster_cols = FALSE, show_rownames = FALSE, cutree_rows = 12)  
dev.off()

range(data_subset_norm)

name(data_subset_norm)
###generanting Upset plots (same idea of the Venn diagram)

library(UpSetR)

data_time_up=fromList(list_up)
data_time_down=fromList(list_down)
setwd("/Users/manuel/Documents/posdoc_berkeley/sombrero_proyecto/")
x11()
pdf("Up_regulated_sombrero.pdf") 
upset(data_time_up, main.bar.color = "SteelBlue", sets.bar.color = "red", nsets = 10,  
      mb.ratio = c(0.6, 0.4), nintersects = 15 , order.by = c("freq"), decreasing = c(TRUE,TRUE), matrix.color="black")
dev.off()

pdf("Down_regulated_sombrero.pdf") 
upset(data_time_down, main.bar.color = "SteelBlue", sets.bar.color = "green", nsets = 10, nintersects = 15, 
      mb.ratio = c(0.6, 0.4), order.by = c("freq"), decreasing = c(TRUE,TRUE), matrix.color="black")
dev.off()





