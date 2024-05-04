


exxx = ex_deg #差异表达矩阵

grey_data1 <- read.table("brown.txt", header = TRUE, stringsAsFactors = FALSE)

#rownames(grey_data) = grey_data[, 1]
#grey_data <- grey_data[!duplicated(grey_data[, 1]), ]

same = exxx[match(grey_data1[, 1],exxx$symbol),]
#对于多个探针对应一个基因的情况，取表达量最大的那个探针
same1 = na.omit(same)

write.table(same1, file = "Candidate signature genes.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


df2 = df
# 检查df2的列名是否可以在same1的symbol列找到
matching_cols <- colnames(df2) %in% same1$symbol

# 筛选df2的列名
df2_filtered <- df2[, matching_cols, drop = FALSE]
df2 = df2_filtered

#df2 <- df2[, match(same1$symbol, colnames(df2))]
colnames(df2)
same1$symbol
# 添加名为"group"的列，并将其值设为"group"的值
df2 <- cbind(group = group, df2)
#ex_deg = same1


#随机森林
set.seed(11)
library(caret)
index = createDataPartition(df2$group,
                            p=0.7,
                            list=FALSE)

train = df2[index,]
test = df2[-index,]
set.seed(123)
x = train %>% dplyr::select(-group)
y = factor(train$group)
rf = randomForest(x,y,data=train,importance = TRUE)
rf
dev.off()
dev.new()
plot(margin(rf,train$group),main="观测值被预测正确的概率图")
#使用验证集
test_predict = predict(rf,test)
compare_test = table(test$group,test_predict,dnn=c("Actual","Predicted"))
compare_test
importance_otu = importance[order(importance$MeanDecreaseAccuracy,decreasing = TRUE),]




the_impotance = svm_rfe_ranking[match(same1$symbol,svm_rfe_ranking$var),]





# 获取df2的列名
column_names <- colnames(df2)
file_path <-'C:/Users/瑶瑶/Desktop/6666.txt'
# 将列名写入文本文件
write.table(column_names, file = file_path, row.names = FALSE, col.names = FALSE)




# 加载VennDiagram包
#install.packages("ggvenn")
library(ggvenn)
# 创建Venn图
genev = list(DEGs=ex_deg$symbol,MEbrown=grey_data1[, 1])
#ggvenn(genev,fill_color = c("#0073C2FF","#EFC000FF"),stroke_size = 0.5,set_name_size = 4)
# 假设genelist是一个包含集合的列表
ggvenn(genev, fill_color = c("#FF0000", "#0073C2FF"), stroke_size = 1) +
  theme(
    text = element_text(size = 20),  # 全局字体大小
    plot.title = element_text(size = 15)  # 标题字体大小，如果有的话
  )

# 首先创建维恩图对象
venn_plot <- ggvenn(genev, fill_color = c("#FF0000", "#0073C2FF"), stroke_size = 1) +
  theme(
    text = element_text(size = 20),  # 全局字体大小
    plot.title = element_text(size = 15)  # 标题字体大小，如果有的话
  )

# 使用ggsave函数保存图像
ggsave(filename = "Venn_diagram.png", plot = venn_plot, dpi = 400, width = 5.5, height = 4)



random = data.frame("BEX4","AP3B2","ARHGEF12","GOLIM4","CADM2","ARF5", "AMPH","ANGPT2",               
                    "CTNS","CCSER2","ASPHD2","BRINP1","GNG3","AP2A2","BEX2","AGAP2",                
                    "ACADVL","ARRB1","ARAF","ATP2A2","C14orf93","CHCHD5","ANXA1")
Randomforest = data.frame("CADM2","ARRB1","AFTPH","AP2A2","BEX4","AP3B2","ASPHD2",
                          "BRINP1","CTNS","ACAD10","ARAF","CCSER2","GOLIM4","BEX2",                 
                          "EEF1A2","ANXA1","ACADVL","AGAP1","ARHGEF12","ATP1A3",               
                          "ANGPT2","COL4A5")
MCC = data.frame("ARHGEF25","AP2A2","GNG3","ARHGEF1","AP3B2","AP5M1","ARHGEF12",
                 "AFTPH","ACOX3","CDYL")

svmm = data.frame("AP3B2","ACADVL","ARF5","CCSER2","AFTPH","BEX2","ASPHD2","AMPH",
                  "AP2A2","ARRB1","CTNS","CADM2","ARAF","ANXA1","ANGPT2","BEX4",
                  "GNG3","GOLIM4","ARHGEF12","BRINP1")
genelist = list(随机森林=Randomforest[1,],cytohHubba = MCC[1,])
ggvenn(genelist,fill_color = c("#FF0000","#0073C2FF"),stroke_size = 0.5,set_name_size = 8)
ggvenn(genelist, fill_color = c("#00FF00","#800080"), stroke_size = 0.5, set_name_size = 8)
ggvenn(genelist, fill_color = c("#5F9EA0","#FF7F50"), stroke_size = 0.5, set_name_size = 8)

library(ggvenn)

# 假设genelist是一个包含集合的列表
ggvenn(genelist, fill_color = c("#FF0000", "#0073C2FF"), stroke_size = 1) +
  theme(
    text = element_text(size = 20),  # 全局字体大小
    plot.title = element_text(size = 20)  # 标题字体大小，如果有的话
  )

genelist


write.table(ex_deg, file = "DEGs.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

