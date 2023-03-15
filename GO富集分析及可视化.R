##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############GO富集分析学习链接：https://www.jianshu.com/p/6b0b5c55058f
##############富集分析输入数据：基因的entrezid ID（富集分析指定使用entrezid ID）
##############富集分析就是看你提供的基因主要集中分布在哪几个通路。
##############Y叔富集分析包，clusterProfiler 让富集分析变得简单易行
##############可以用来做各种富集分析，如KEGG，GO,GSEA富集分析等
#############R包学习地址：https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
##############加载数据
deg.data <- data.table::fread('DEG.csv',data.table = F)
#############安装R包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")

#install.packages('ggthemes')
#install.packages('stringr')

##############加载所需要R包
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
deg.data$Group[which((deg.data$adj.P.Val) < 0.01 & (deg.data$logFC > 3) )] <- 'Up'
deg.data$Group[which((deg.data$adj.P.Val) < 0.01 & (deg.data$logFC < -3) )] <- 'Down'

gene_up <- deg.data$id[deg.data$Group == 'Up']
gene_down <- deg.data$id[deg.data$Group == 'Down']

gene <- c(gene_down)#,gene_down)


###############基因名称ID的转换 clusterprofiler R包中的自带函数一键转换
#    bitr(geneID, fromType,toType, OrgDb)

gene <- bitr(gene,
             fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb = 'org.Hs.eg.db')  ####人类的注释包是org.Hs.eg.db,小鼠是org.Mm.eg.db

keytypes(org.Hs.eg.db) ###查看各种各样的基因ID

########注意在进行ID转换时，会损失部分基因，不是所有基因的ID都能够对应的

###############开始富集分析
ego <- enrichGO(gene = gene$ENTREZID,
                OrgDb = org.Hs.eg.db,  ###此处是富集分析所对应的人类数据库
                ont = 'ALL',           ###Go富集分析的那三个子类 cc bp mf ont参数：One of "BP", "MF", and "CC" 或者‘ALL’
                readable = T)          ###让你的结果可读

#ego_BP <- enrichGO(gene = gene_diff,
#                   OrgDb= org.Hs.eg.db,
#                   ont = "BP",         ###只做BP
#                   readable = T)
#
################查看富集结果
##ONTOLOGY    代表那三个子类 BP CC MF
##Description 富集到三个子类的名称
##GeneRatio   25/411     该通路的差异基因数，对应到数据库中的差异基因总数 
##BgRatio     111/18723  该通路的所有基因数目，所有通路总共的基因数目
##输入560个，只有411个说明部分基因没有被数据库收录，没被收入的基因说明还没有研究清楚他们的功能需要人类继续开发
##count指的时该部分富集的基因数目 P值代表富集的程度

################GO富集结果的保存
#save(ego,file = 'ego.Rda')
a <- ego@result
#write.csv(a,file = 'ego_up.csv',row.names = T)

##############富集分析结果的可视化 
##############可视化链接：https://www.jianshu.com/p/13a8fd09169b
##############            https://www.jianshu.com/p/18a0a696a4cc
##############            https://zhuanlan.zhihu.com/p/133342400
#柱状图
#barplot(ego) 
#气泡图
#dotplot(ego)

#######分层柱状图
barplot(ego,drop = TRUE,showCategory =5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

#####描述条目太长怎么办？
#####必应一下！
#####链接：https://www.jianshu.com/p/f7ff4d2c2976
########## https://www.jianshu.com/p/6486c8e9a4d4
barplot(ego,drop = TRUE,label_format = 100,showCategory =5,split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale='free')+
  theme(legend.text = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 20))
  


#######分层气泡图
dotplot(ego,label_format = 100,showCategory = 5,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')+
  theme(legend.text = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 20))

#######显示差异基因倍数的圈图
colnames(gene)[1] <-'id' 
deg <- deg.data[,c(1,3)]
gene <- merge(deg,gene,by='id')

geneList <-  gene$logFC
names(geneList) <- gene$ENTREZID

pdf(file = 'GO_cneplot.pdf',width = 10,height = 10)
cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
dev.off()

pdf(file = 'GO_cneplot2.pdf',width = 13,height = 10)
cnetplot(ego, showCategory = 3,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
dev.off()

##########不使用差异基因配色
enrichplot::cnetplot(ego,circular=TRUE,colorEdge = TRUE)



