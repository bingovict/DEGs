##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############加载数据
deg.data <- data.table::fread('DEG.csv',data.table = F)
##############加载所需要R包，上节课已经安装过了
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)

###############获取需要富集分析的基因
#############修改列名
colnames(deg.data)[3] <- 'logFC'
colnames(deg.data)[7] <- 'adj.P.Val'
############新增一列，区分上调和下调的基因
deg.data$Group <- 'Not-significant'
deg.data$Group[which((deg.data$adj.P.Val) < 0.01 & (deg.data$logFC > 2) )] <- 'Up'
deg.data$Group[which((deg.data$adj.P.Val) < 0.01 & (deg.data$logFC < -2) )] <- 'Down'

gene_up <- deg.data$id[deg.data$Group == 'Up']
gene_down <- deg.data$id[deg.data$Group == 'Down']

gene_all <- c(gene_up,gene_down)

###########讲Symbol ID 转换为Entrezid ID  小鼠的为：org.Mm.eg.db
gene_up <- bitr(gene_up, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

gene_down <- bitr(gene_down, 
                  fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Hs.eg.db")

gene_all <- bitr(gene_all, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb="org.Hs.eg.db")


###########开始作KEGG富集分析
kk_up <- enrichKEGG(gene = gene_up$ENTREZID,
                    organism = 'hsa',              ###小鼠的为‘mmu’
                    pvalueCutoff = 0.05)           ###此处的pvaleCutoff可以设置为1
save(kk_up,file = 'kk_up.Rda')

kk_down <- enrichKEGG(gene = gene_down$ENTREZID,
                       organism = 'hsa',           ###小鼠的为‘mmu’
                       pvalueCutoff = 0.05)
save(kk_down,file = 'kk_down.Rda')

kk_all <- enrichKEGG(gene = gene_all$ENTREZID,
                     organism = 'hsa',             ###小鼠的为‘mmu’
                     pvalueCutoff = 0.05) 
save(kk_all,file = 'kk_all.Rda')

#######加载KEGG富集结果
load('kk_up.Rda')
load('kk_down.Rda')
load('kk_all.Rda')

#######让KEGG结果变得可读：https://www.jianshu.com/p/f7e220ccdf4b
kk_up <-  setReadable(kk_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kk_down <- setReadable(kk_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kk_all <- setReadable(kk_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

###### 保存分析的结果
a <- kk_all@result
write.csv(a,'kk_all.csv',row.names = F)



#######查看富集分析显著性
table(kk_up@result$p.adjust < 0.05)     ####实在不行就使用P值
table(kk_down@result$p.adjust < 0.05)
table(kk_all@result$p.adjust < 0.05)

######KEGG结果的可视化
barplot(kk_up)
barplot(kk_down)
barplot(kk_all)


dotplot(kk_up)
dotplot(kk_down)
dotplot(kk_all)

pdf(file = 'kk_all_cneplot.pdf',width = 15,height = 15)
enrichplot::cnetplot(kk_all,circular=TRUE,colorEdge = TRUE)
dev.off()

###########KEGG本地化，下载KEGG包，再也不用担心网络问题
########### KEGG.db_1.0.tar.gz R包已经放在文件夹
#install.packages("KEGG.db_1.0.tar.gz", repos=NULL)### 安装该R包
############加载该R包，妈妈再也不用担心我的网络问题
library(KEGG.db) 
kk_up <- enrichKEGG(gene = gene_up$ENTREZID,
                    organism = 'hsa',              ###小鼠的为‘mmu’
                    pvalueCutoff = 0.05,
                    use_internal_data =T)  ####这句代码的意思是不使用在线网站

