##############安装本次课程所需要的R包。
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("remotes",force = TRUE)

BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data",force = T)

BiocManager::install("BioinformaticsFMRP/TCGAbiolinks",force = T)

BiocManager::install("SummarizedExperiment",force = TRUE)

##############TCGAbiolink包的学习地址
##############https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
######设置工作目录
setwd('D:\\R\\正式科研代码')
######加载R包  如果加载失败，说明没有下载R包成功
library(SummarizedExperiment)
library(TCGAbiolinks)

######使用TCGAbiolink包下载转录组信息，三个步骤，查询，下载，整理简单易行
######第一步，查询TCGA官网需要下载的癌症项目的信息
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",   
                  workflow.type = "STAR - Counts", 
                  legacy = FALSE)

#####project对应的是TCGA官网的癌症类型
#####data.category对应的是需要下载的转录组测序数据
#####data.type对应的是基因表达定量

#####下载，将文件下载到工作目录中，并在目录下创建GDCdata的文件夹
#####下次使用时，不再需要运行该句代码，重复下载
GDCdownload(query = query)

#####整理，将下载好的文件整理成SE格式
SE=GDCprepare(query = query)

names(SE@assays)    ###assays包含的是打包的各种数据类型的基因表达量的信息
names(SE@rowRanges) ###rowRanges包含基因注释信息
names(SE@metadata)

#####获取我们所需要的基因注释信息
mydf <- as.data.frame(rowRanges(SE))
###查看基因类型
table(mydf$gene_type)

#####获取所需的基因表达量的数据类型：tpm_unstrand，fpkm_unstrand，unstranded（counts）
gene_expr_counts <- as.data.frame(assay(SE,i='unstranded'))
gene_expr_fpkm <- as.data.frame(assay(SE,i='fpkm_unstrand'))
gene_expr_tpm <- as.data.frame(assay(SE,i='tpm_unstrand'))

#####Ensemble ID转换
identical(mydf$gene_id,rownames(gene_expr_counts))

gene_expr_counts <- cbind(type=mydf$gene_type,ID=mydf$gene_name,gene_expr_counts)
mydata_counts <- gene_expr_counts[,-1]

gene_expr_fpkm <- cbind(type=mydf$gene_type,ID=mydf$gene_name,gene_expr_fpkm)
mydata_fpkm <- gene_expr_fpkm[,-1]

gene_expr_tpm <- cbind(type=mydf$gene_type,ID=mydf$gene_name,gene_expr_tpm)
mydata_tpm <- gene_expr_tpm[,-1]
#####数据保存
save(mydata_counts,mydata_fpkm,mydata_tpm,file = 'TCGA_STAD.Rda')







