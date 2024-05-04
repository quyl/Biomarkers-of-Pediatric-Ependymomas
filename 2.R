
options(stringsAsFactors = F)
library(GEOquery)
library(tidyr)
library(dplyr)
library(limma)
library(WGCNA)
library(stringr)
library(dplyr)
#install.packages("AnnoProbe")
library(AnnoProbe)

gset<-eSet
gset[[1]]
a=gset[[1]]
dat=exprs(a)
dim(dat)
#dat[1:4,1:4]

#dat=normalizeBetweenArrays(dat)

boxplot(dat[,1:4],las=2)
metadata=pData(a) #矩阵的临床信息
library(stringr)
ids=idmap('GPL6947','soft')

idddd = idmap()
fd1=ids


load(file = "step2output.Rdata")
rownames(exp_sy) = exp_sy$symbol
exp_sy=exp_sy[,-1]
datExpr0 = as.data.frame(t(exp_sy))

gsg=goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0),method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
#plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h=1600000,col="red")

clust=cutreeStatic(sampleTree,cutHeight = 160000,minSize = 10)
table(clust)
keepSamples=(clust==1)
datExpr0=datExpr0[keepSamples,] #现在是去除了异常值之后放入输入数据

p=datExpr0
metadata <- metadata[rownames(metadata) %in% rownames(datExpr0), ]
p=p[rownames(p) %in% rownames(metadata),]
datExpr0 =p
row_names <- rownames(metadata)
col1 <- metadata[, 8] #取出临床信息里想要的那一列信息
traitData <- metadata[, -c(1:7)]
traitData <- traitData[, -c(3:31)]
traitData[, 2] <- traitData[, 1]

colnames(traitData)[1] <- "cancer"
colnames(traitData)[2] <- "normal"
traitData$cancer <- ifelse(grepl("tumour",traitData$cancer),1,0)
traitData$normal <- ifelse(grepl("BioChain",traitData$normal),1,0)
#traitData$cancer <- ifelse(traitData$cancer == "Paediartric brain tumour", 1, 0)
#traitData$normal <- ifelse(traitData$normal == "Non Tumoral Lung", 1, 0)

sampleTree2 = hclust(dist(datExpr0),method="average")
plot(sampleTree2)


sameSample = intersect(rownames(datExpr0),rownames(traitData))
sameSample
datExpr0=datExpr0[sameSample,]
datTraits=traitData[sameSample,]
traitColors = numbers2colors(datTraits, signed = FALSE)
##画图，红色表示阳性，白色表示阴性
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr0, datTraits, file = "01-dataInput.RData")

load(file="01-dataInput.RData")


#powers = c(c(1:10), seq(from = 12, to=20, by=1))
powers = c(1:20)
sft=pickSoftThreshold(datExpr0,powerVector = powers,verbose = 5)





par(mfrow = c(2,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red") #可以修改
# 保存图片
# 设置图片的保存路径和清晰度
png(file="Scale_independence.png", width=400*6, height=400*5, res=400)
# 重复绘图代码
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.9, col="red")
# 关闭图片保存
dev.off()


dev.off()
dev.new()
###平均连通性与power值散点图
png(file="Mean Connectivity.png", width=400*6, height=400*5, res=400)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 设置保存图形的参数，包括文件名、尺寸和分辨率
png(filename = "Mean_connectivity.png", width = 400*6, height = 400*5, res = 400)

# 绘制平均连接度图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",  # x轴标签
     ylab="Mean Connectivity",  # y轴标签
     type="n",  # 不绘制点
     main = paste("Mean connectivity"))  # 图片标题
# 在图上添加文本，显示不同的软阈值（power）
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels=powers,  # 显示的文本，powers应该是一个包含了所有软阈值的向量
     cex=cex1,  # 文本大小，cex1需要预先定义
     col="red")  # 文本颜色

# 关闭图形设备，完成图片保存
dev.off()



AsoftPower=sft$powerEstimate
sft
softPower = sft$powerEstimate
adjacency = adjacency(datExpr0,power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
save(dissTOM,file="step7output.Rdata")
rm(list = ls())
load("step7output.Rdata")

geneTree = hclust(as.dist(dissTOM),method = "average")
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


minModuleSize = 50      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#对每一个基因块附上颜色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
m=plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 保存基因聚类图
png(file = "Gene_clustering_TOM.png", width = 400*5, height = 400*3, res = 400)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# 保存带有动态剪切的树状图和模块颜色的图
png(file = "Gene_dendrogram_and_module_colors.png", width = 400*6, height = 400*3, res = 400)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# 设置全局图形参数使用中文字体
par(family = "SimSun")
pdf("julei.pdf")
plot1=plotDendroAndColors(geneTree, dynamicColors, "基因模块",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "聚类树状图")




ggsave("myplot.pdf", plot1, width = 8, height = 8, units = "cm")

load("01-dataInput.RData")
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")





nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue =
  
  
  


corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)


dev.off()
dev.new()
par(cex=1)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# 设置图像保存参数，包括文件名、尺寸和分辨率
png(file = "Module_trait_relationships1.png", width = 10*400, height = 9*400, res = 400)

# 设置图形参数
par(cex=1)
par(mar = c(10, 10, 3, 3))

# 绘制带标签的热图
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# 关闭图像设备
dev.off()




modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")



traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = dynamicColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}


for (mod in 1:nrow(table(dynamicColors)))
{  
  modules = names(table(dynamicColors))[mod]
  probes = colnames(datExpr0)
  inModule = (dynamicColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}


