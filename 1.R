#第一步：start GEO
#下载需要的数据
library(GEOquery)
gse = "GSE42656"
Sys.setenv("VROOM_CONNECTION_SIZE" = 2621440*2)
eSet<-getGEO(gse,destdir = '.',getGPL = F)

#提取表达矩阵exp
exp<-exprs(eSet[[1]])

#提取临床信息
pd <- pData(eSet[[1]])
load("Step1output.Rdata")
#调整exp的行名和pd的列名完全一致
p=identical(rownames(pd),colnames(exp));p
if(!p) exp=exp[,match(rownames(pd),colnames(exp))]

#提取芯片平台编号
gpl <- eSet[[1]]@annotation
save(eSet,gse,pd,exp,gpl,file = "Step1output.Rdata")


#第二步：group id
#group_list（实验分组）和ids（芯片注释）
rm(list = ls())
load("Step1output.Rdata")#加载过来之后上述的数据又回来啦
library(stringr)
library(dplyr)

#取出自己需要的胶质母细胞瘤的数据和正常的数据
#删除pd的前20行
pd = pd[-(1:25),]
#删除6:29
pd = pd[-(15:23),]
#删
pd = pd[-(15:23),]

gr = pd$source_name_ch1
k1=str_detect(gr,"tumour");table(k1)
k2=str_detect(gr,"BioChain");table(k2)

group=ifelse(k1,"cancer",'normal');table(group)
group
group=factor(group,levels = c("normal","cancer"))
table(group)
grouppp = ifelse(k1,"1","0");table(grouppp)
grouppp

#为了提取需要的部分：删除exp矩阵的列名不在pd矩阵的行名中的exp的所有列
# 获取exp矩阵中列名在pd矩阵行名中的列索引
valid_columns <- which(colnames(exp) %in% rownames(pd))
# 仅保留有效的列
exp_filtered <- exp[, valid_columns]
exp = exp_filtered



#install.packages("DT")
#探针注释
library(devtools)
library(GEOquery)
#install_github("jmzeng1314/idmap3")
#install_github("jmzeng1314/idmap2")
#install_github("jmzeng1314/idmap1")
library(idmap3)
library(idmap2)
library(idmap1)
library(dplyr)

temp=getIDs("GPL6947")
#ttttt = getIDs("GPL6947")

rowname=data.frame(probe_id=rownames(exp))
gene <- left_join(rowname,temp,by='probe_id')
gene<-na.omit(gene) #删除缺失值
exp<-exp[match(gene$probe_id,rownames(exp)),]
#对于多个探针对应一个基因的情况，取表达量最大的那个探针
result <- data.frame(ex=rowMeans(exp),symbol=gene$symbol)
result <- arrange(result,symbol,desc(ex))
de <- result[!duplicated(result$symbol),]
exp <- exp[match(rownames(de),rownames(exp)),]

exp = normalizeBetweenArrays(exp)
exp  = log2(exp)

exp_sy <- data.frame(symbol=de$symbol,exp)

save(eSet,pd,exp,exp_sy,group,gpl,file="step2output.Rdata")



rm(list = ls())
#load(file = "Step1output.Rdata")#step1的数据有误，不加载
load(file = "step2output.Rdata")
#exp = exp_sy[,-1]
group_list=group
table(group_list)

dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra)
#install.packages("FactoMineR")
#install.packages("htmlools")
#install.packages("factoextra")
#画pca
dat.pca<-PCA(dat,graph = FALSE)
pca_plot<-fviz_pca_ind(dat.pca,
                       geom.ind = "point",
                       col.ind = group_list,
                       addEllipses = TRUE,
                       legend.title="Groups")
pca_plot
save(pca_plot,file = "pca_plot.Rdata")



#第四步 ：DEG
rm(list = ls())
load(file = "step2output.Rdata")
#差异分析，用limma
#需要表达矩阵和group_list
exp=exp_sy[,-1]

#exp_sy = 

rownames(exp)= exp_sy$symbol
group_list=group
group
library(limma)
design=model.matrix(~group_list)
design
fit = lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef = 2,number = Inf)
nrDEG11 = na.omit(deg)

#为deg数据框添加几列
#1.加probe_id列，把列名变成一行
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)
#2.加symbol列，火山图要用
deg<-cbind(deg,symbol=exp_sy$symbol)
head(deg)

#3.加change列，标记上下调基因
logFC_t = 1
P.Value_t=0.05
k1=(deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2=(deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,
                "down",
                ifelse(k2,
                       "up",
                       "stable"))
table(change)
deg <- mutate(deg,change)
ex_deg<- deg[k1|k2,]    #这个ex_deg是差异表达矩阵(包含上调和下调的基因)
save(group_list,deg,logFC_t,P.Value_t,ex_deg,file="step4output.Rdata")

#检测差异表达矩阵是否需要归一化
boxplot(exp,las=2)

#单基因表达量
library(ggplot2)
library(ggpubr)
exp_unique = exp_sy
rownames(exp_unique) = exp_unique$symbol
exp_unique = exp_unique[,-1]
gene2=exp_unique["AFTPH",]
geneee = data.frame(gene=t(gene2),group=group)

ggplot(data=geneee, aes(x=group, y=AP3B2, fill=group)) +
  geom_boxplot() +
  geom_point(position = "jitter") +
  stat_compare_means() +
  theme(text = element_text(size = 18),  # 全局字体大小
        axis.text.x = element_text(size = 18),  # 横坐标字体大小
        axis.text.y = element_text(size = 18),  # 纵坐标字体大小
        legend.key.size = unit(1, "cm"),  # 图例框的大小
        legend.text = element_text(size = 18))  # 图例文字的大小

# 创建箱线图对象
boxplot <- ggplot(data=geneee, aes(x=group, y=AFTPH, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  stat_compare_means() +
  theme(text = element_text(size = 18),  # 全局字体大小
        axis.text.x = element_text(size = 18),  # 横坐标字体大小
        axis.text.y = element_text(size = 18),  # 纵坐标字体大小
        legend.key.size = unit(1, "cm"),  # 图例框的大小
        legend.text = element_text(size = 18))  # 图例文字的大小

# 使用ggsave函数保存图像
ggsave(filename = "Boxplot_AFTPH.png", plot = boxplot, dpi = 400, width = 8, height = 8)




library(tidyverse)
#install.packages("tidyverse")
load(file = "step2output.Rdata")
load(file = "step4output.Rdata")

group




df <- exp_sy[match(ex_deg$symbol,exp_sy$symbol),]
df <- df[!duplicated(rownames(df)), ]
rownames(df)<-df$symbol
df <- data.frame(t(df[,-1])) %>%
  cbind(group,.)

#这里还需要构建一个名为nrDEG的矩阵，用来画热图和火山图
nrDEG = ex_deg
rownames(nrDEG) = nrDEG$symbol
head(nrDEG)
#再删除7-9列
nrDEG = nrDEG[,-(7:8)]



#随机森林，在取完交集之后用
library(randomForest)
set.seed(123)
#ind = sample(2,nrow(df2),replace = TRUE,prob = c(0.8,0.2))
#forest = randomForest(group~.,data = df2[ind==1,],ntree=10,nPerm=10,mtry=3,proximity=TRUE,importance=TRUE)
#fit.forest<-randomForest(group~.,data = df2,ntree=300)
#forest
#fit.forest
#plot(fit.forest)
#plot(forest)

#ffff = predict(forest,df2[ind==2,])
#table(observed = df2[ind==2,"group"],predicted=ffff)
#importance(forest,type = 1)

#找出使模型准确率达到最优需要的树的数量
#which.min(fit.forest$err.rate[,1])

#set.seed(1423)
#fit.forest_2 <- randomForest(group~., data = df2, ntree=8)
#fit.forest_2 #(122+168)/(122+168+14+9)



#建立随机森林模型-其中一个代码
ozone.rf <- randomForest(group ~ ., 
                         data=df2, 
                         mtry=3,
                         importance=TRUE, 
                         proximity=TRUE,
                         na.action=na.omit)
Randomforest = ozone.rf
ozone.rf
plot(Randomforest)
importance(ozone.rf)
varImpPlot(ozone.rf)
varImpPlot(Randomforest)
#varImpPlot(fit.forest_2)
varImpPlot(forest)
#result_df <- data.frame(importance(fit.forest_2,type = 2)) %>%
  #arrange(desc(.))
result_df = data.frame(importance(ozone.rf,type=1)) %>% arrange(desc(.))
result_df
#选取几个比较重要的基因放入后续的分析
impo_df = filter(result_df,MeanDecreaseAccuracy>3)
#impo_df <- filter(result_df,MeanDecreaseGini>0.2)
impo_df



# 计算变量重要性并创建一个数据框
importance_df = data.frame(importance(ozone.rf, type=1))

# 按 MeanDecreaseAccuracy 降序排列
importance_df = importance_df %>% arrange(desc(MeanDecreaseAccuracy))

# 选取 MeanDecreaseAccuracy > 3 的基因
impo_df = filter(importance_df, MeanDecreaseAccuracy > 3)

# 如果 impo_df 是一个数据框，它需要被转换成一个随机森林对象的形式
# 以便 varImpPlot 函数可以正确地处理。以下代码仅适用于已经是随机森林对象的 impo_df
# 如果 impo_df 不是随机森林对象，请参考随机森林包文档进行必要的转换。

# 绘制只包含 MeanDecreaseAccuracy > 3 的基因的变量重要性图
if (nrow(impo_df) > 0) {
  # 仅在 impo_df 不为空时绘图
  varImpPlot(Randomforest, sort=TRUE, n.var=nrow(impo_df), main="Variable Importance")
}





# 计算变量重要性并创建一个数据框
importance_df = data.frame(importance(ozone.rf, type=1))

# 按 MeanDecreaseAccuracy 降序排列
importance_df = importance_df %>% arrange(desc(MeanDecreaseAccuracy))

# 选取 MeanDecreaseAccuracy > 3 的基因
impo_df = filter(importance_df, MeanDecreaseAccuracy > 3)
# 设定保存图片的参数，包括文件名、宽度、高度和分辨率
png(filename = "variable_importance.png", width = 6, height = 6, units = "in", res = 400)

# 绘制只包含 MeanDecreaseAccuracy > 3 的基因的变量重要性图
# 并且去掉 MeanDecreaseGini
if (nrow(impo_df) > 0) {
  # 仅在 impo_df 不为空时绘图
  varImpPlot(Randomforest, sort=TRUE, n.var=nrow(impo_df), type=1, main="Variable Importance")
}

# 关闭设备，完成图像保存
dev.off()


# 绘制只包含 MeanDecreaseAccuracy > 3 的基因的变量重要性图
# 并且去掉 MeanDecreaseGini
if (nrow(impo_df) > 0) {
  # 仅在 impo_df 不为空时绘图
  varImpPlot(Randomforest, sort=TRUE, n.var=nrow(impo_df), type=1, main="Variable Importance")
}








#找到这几个重要的基因的信息
impo_exp <- exp_sy[match(rownames(impo_df),exp_sy$symbol),]
rownames(impo_exp) <- impo_exp[,1]

n=impo_exp[,-1]
dim(n)

#这部分先不画
#绘制随机森林筛选出来的前20个hub基因的热图
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n)
library(ggplotify)
heatmap_plot <- as.ggplot(pheatmap(n,show_colnames = F,
                                   show_rownames = T,
                                   cluster_cols = T,
                                   scale = "row",
                                   annotation_col = annotation_col))
rownames(impo_exp)
save(list = ls(),file = "step5output.Rdata")



limma_sigGene = nrDEG11[nrDEG11$change !="normal",1]

library(png)
#这里得造出来一个expp用来存放exp_cy的改变
expp=exp_sy
rownames(expp) = expp$symbol
expp = expp[,-1]
options(stringsAsFactors = F)
#需要画一个前五十个差异表达基因的热图
library(pheatmap)
chayi = ex_deg
# 按照logFC列的绝对值大小排序矩阵
chayi <- chayi %>%
  arrange(desc(abs(logFC)))
rownames(chayi)=chayi$symbol

# 选取chayi矩阵中change列为up的前25行
up_rows <- chayi[chayi$change == "up", ][1:50, ]
down_rows = chayi[chayi$change=="down",][1:50,]
# 合并前25行为up和down的数据
choose <- rbind(up_rows, down_rows)
choose_gene = rownames(choose)
choose_gene = rownames(chayi)

choose_matrix = xxx[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix)


library(pheatmap)
rownames(same1) = same1$symbol
choose_gene1 = head(rownames(same1),85)
xxx = exp_sy
rownames(xxx) = xxx$symbol
xxx = xxx[,-1]
choose_matrix2 = xxx[choose_gene1,]
choose_matrix2 = t(scale(t(choose_matrix2)))
pheatmap(choose_matrix2)
n=pheatmap(choose_matrix2, fontsize_row = 7, fontsize_col = 9)
ggsave("85retu.png", n, width = 8, height = 8, dpi = 400)
# 保存当前热图为PNG格式
#png("heatmap.png")
#print(heatmap)  # 将热图打印到文件
#dev.off()  # 关闭绘图设备

library(pheatmap)
rownames(ex_deg) = ex_deg$symbol
choose_gene1 = head(rownames(ex_deg),1654)
xxx = exp_sy
rownames(xxx) = xxx$symbol
xxx = xxx[,-1]
choose_matrix1 = xxx[choose_gene1,]
choose_matrix1 = t(scale(t(choose_matrix1)))
pheatmap(choose_matrix1)
m=pheatmap(choose_matrix1, fontsize_row = 4, fontsize_col = 8)
ggsave("1000retu.png", m, width = 7, height = 7, dpi = 400)
#差异基因火山图

DEG=deg
logFC_cutoff <- logFC_t

this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='up',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='down',])
)
this_tile
head(DEG)
library(ggplot2)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)
ggsave("1000huoshan.png", g, width = 7, height = 7, dpi = 400)



#单基因表达量
library(ggplot2)
library(ggpubr)
exp_unique = exp_sy
rownames(exp_unique) = exp_unique$symbol
exp_unique = exp_unique[,-1]
gene2=exp_unique["AP2A2",]
geneee = data.frame(gene=t(gene2),group=group)
ggplot(data=geneee, aes(x=group, y=AP2A2, fill=group)) +
  geom_boxplot() +
  geom_point(position = "jitter") +
  stat_compare_means() +
  theme(text = element_text(size = 18),  # 全局字体大小
        axis.text.x = element_text(size = 18),  # 横坐标字体大小
        axis.text.y = element_text(size = 18),  # 纵坐标字体大小
        legend.key.size = unit(1, "cm"),  # 图例框的大小
        legend.text = element_text(size = 18))  # 图例文字的大小


library(ggplot2)
library(ggpubr)

# 假设exp_unique是你的基因表达数据框，group是样本分组
genes_of_interest <- c("AP3B2", "AP2A2", "ARHGEF12", "AFTPH")

# 提取感兴趣的基因
exp_selected <- exp_unique[genes_of_interest,]

# 将数据转换为长格式以适应ggplot
#geneee_long <- reshape2::melt(t(exp_selected))
geneee_long = data.frame(gene=t(exp_selected),group=group)
colnames(geneee_long) <- c("AP3B2", "AP2A2", "ARHGEF12", "AFTPH","group")



# 绘制四个基因的表达量盒须图
par(mfrow=c(2,2))
for (i in 1:length(exp_selected)) {
    ggplot(data=geneee_long, aes(x=group, y=exp_selected[i,], fill=group))

}




p <- ggplot(data=geneee_long, aes(x=group, y=expression, fill=group)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  stat_compare_means() +
  facet_wrap(~ gene, scales = "free") +
  theme(text = element_text(size = 18),  # 全局字体大小
        axis.text.x = element_text(size = 18),  # 横坐标字体大小
        axis.text.y = element_text(size = 18),  # 纵坐标字体大小
        legend.key.size = unit(1, "cm"),  # 图例框的大小
        legend.text = element_text(size = 18))  # 图例文字的大小

# 打印图表
print(p)

