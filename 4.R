
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
#install.packages("GOplot")
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(BiocManager)
#BiocManager::install("topGO")
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
#BiocManager::install("ComplexHeatmap")



info = same1
info = info[,-2]
info = info[,-2]
info = info[,-(4:5)]
info = info[,-5]
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
#gene ID转换
gene <- bitr(info$symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)


GO<-enrichGO(gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)
GO


KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
#gsea
#names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
names(info) <- c('Log2FoldChange','pvalue','padj','SYMBOL')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
GSEA_input
names(GSEA_input) = info_merge$ENTREZID

GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)#GSEA富集分析

barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
yy=dotplot(KEGG,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
           showCategory = 10) 
# 调整图中的字体大小
yy <- yy + theme(
  text = element_text(size = 1),  # 全局字体大小
  axis.title = element_text(size = 15),  # 轴标题
  axis.text = element_text(size = 180),  # 轴文本
  legend.title = element_text(size = 15),  # 图例标题
  legend.text = element_text(size = 16),  # 图例文本
  strip.text = element_text(size = 15),  # facet标签
  axis.text.y = element_text(size = 18)  # 纵坐标文本字体大小
)
print(yy)
ggsave("kegg.png", yy, width = 8, height = 10, dpi = 500)

mm = dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图


nn=dotplot(GO,x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
        showCategory =5,#只显示前5
        split="ONTOLOGY") + #以ONTOLOGY类型分开
  facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分屏绘图
plot <- nn +scale_x_continuous(breaks = seq(0, 0.06, by = 0.01))

print(plot)
ggsave("plot.png", plot, width = 10, height = 12, dpi = 500)

library(ggplot2)

# 假设GO是你的数据框，包含ONTOLOGY列
# 创建dotplot图
p <- dotplot(GO, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free")

# 调整图中的字体大小
nn <- nn + theme(
  text = element_text(size = 1),  # 全局字体大小
  axis.title = element_text(size = 15),  # 轴标题
  axis.text = element_text(size = 180),  # 轴文本
  legend.title = element_text(size = 15),  # 图例标题
  legend.text = element_text(size = 16),  # 图例文本
  strip.text = element_text(size = 15),  # facet标签
  axis.text.y = element_text(size = 16)  # 纵坐标文本字体大小
)

library(dplyr)
library(ggplot2)








+# 显示图形
print(mm)

dotplot(KEGG)
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE
head(info)

write.table(KEGG$ID, file = "C:/Users/瑶瑶/Desktop/KEGG_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
browseKEGG(KEGG,"hsa04144")#选择其中的hsa05166通路进行展示




genedata<-data.frame(ID=info$SYMBOL,logFC=info$Log2FoldChange)
write.table(GO$ONTOLOGY, file = "F:/毕业论文/室管膜瘤/图图们/85GO_ONTOLOGYs.txt", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


GOplotIn_BPP<-GO[1:4,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[5:12,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MFP<-GO[21:28,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BPP$geneID <-str_replace_all(GOplotIn_BPP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MFP$geneID <-str_replace_all(GOplotIn_MFP$geneID,'/',',')
names(GOplotIn_BPP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MFP)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BPP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MFP$Category = "MF"
circ_BPP<-GOplot::circle_dat(GOplotIn_BPP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MFP<-GOplot::circle_dat(GOplotIn_MFP,genedata) 


chord_BP<-chord_dat(data = circ_BPP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MFP,genes = genedata) 
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size =10,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 17) #GO Term字体大小
ggsave("BP.png", BP, width = 20, height = 10, dpi = 500)
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 10,
        lfc.col = c('red','white','blue'), 
        process.label = 17) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 10,
        lfc.col = c('red','white','blue'), 
        process.label = 17)

GOChord(chord_CC,space = 0.02,gene.order = "LogFC",gene.space = 0.25,gene.size = 5,process.label = 10)








pdf("chord.pdf",height = 13,width = 13)#准备好一块画布。
GOChord(chord_BP,space = 0.02,gene.order = 'logFC',gene.space = 0.25,gene.size = 5)#进行绘图。
dev.off()#然后必须关闭并保存这块画布。


library(GOplot)

# 创建弦图
q <- GOChord(data = chord_CC,
             title = 'GO-Cellular Component',
             space = 0.01,
             limit = c(1,1),
             gene.order = 'logFC',
             gene.space = 0.25,
             gene.size = 13,
             lfc.col = c('red','white','blue'),
             process.label = 40)

# 调整标题字体大小
q <- q + theme(plot.title = element_text(size = 40))  # 设置标题字体大小为40

# 将图例放在图片右上角
q <- q + theme(legend.position = c(1, 1), # 图例位置设置为右上角
               legend.justification = c(1, 1), # 图例对齐设置为右上角
               legend.box.just = "right", # 确保整个图例框靠右对齐
               legend.box.margin = margin(0, 0, 0, 0), # 移除图例周围的空白
               legend.margin = margin(0, 0, 0, 0)) # 移除图例内部标签周围的空白

# 打印图形
print(q)
library(ggplot2)
# 设置文件路径到桌面
file_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "go.png")
# 使用ggsave函数保存图形到桌面
ggsave(file_path, plot = p, dpi = 400)

ggsave("GOChord_plot.png", plot = q, width = 10, height = 8, dpi = 400)


install.packages("GOplot")



chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 
