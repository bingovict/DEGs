##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############加载数据
load('counts_vst.Rda')  ####此处为标准化后的count数据
load('DEG.Rda')

##############筛选差异基因
##############加载R包
library(dplyr)
deg.genes <- DEG %>% 
             filter(padj < 0.01) %>% 
             filter(abs(log2FoldChange) > 3)
deg.genes <- deg.genes[,1]

test <- counts_vst[deg.genes,]

#############制作样本分组
table(substr(colnames(test),14,16))
group <- data.frame(Type=c(rep('Normal',36),rep('Tumor',410)))
rownames(group) <- colnames(test)

###### 绘制自己的热图
library(pheatmap) 
outFile <- "my_heatmap.pdf"      

#绘制
pdf(file=outFile,width=6,height=5.5)
pheatmap(test,
         annotation=group,
         cluster_cols = T,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = F,# 显示行名
         show_colnames = F,# 显示行名
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

#########使用fpkm数据，绘制热图
load('Symbol.Rda')
test2 <- symbol_tpm[deg.genes,]
save(test2,file='test2.Rda')
#test2 <- cbind(id=rownames(test2),test2)
#write.csv(test2,file = 'deg.genes_fpkm.csv',row.names = F)

#test2 <- test2[,-1]
test2 <- as.matrix(test2)
log_test2 <- log2(test2+1)


outFile <- "my_heatmap2.pdf"      
#绘制
pdf(file=outFile,width=6,height=5.5)
pheatmap(log_test2,
         annotation=group,
         cluster_cols = T,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = F,# 显示行名
         show_colnames = F,# 显示行名
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)
dev.off()
