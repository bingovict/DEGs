##############设置工作目录
setwd("D:\\R\\正式科研代码")
##############安装本次课程所需要的R包。
if(!require(DESeq2))BiocManager::install('DESeq2')
##############加载R包
library(DESeq2)
##############加载数据
load('Symbol.Rda')
data <- symbol_counts
#############获取样本分组信息
table(substr(colnames(data),14,16))
group <- c(rep('Normal',36),rep('Tumor',410))

sample <- data.frame(row.names=colnames(data),condition=group)

############数据打包成DESeq格式的数据集
dds <- DESeqDataSetFromMatrix(countData = data,    ###此处放我们的count数据
                              colData = sample,    ###此处放我们的分组信息表格
                              design = ~condition) ###此处放sample表格中的condition列名

############此处的警告是自动将sample表格中的condition这一列的字符转换为因子


############过滤掉低表达的counts值,count函数过滤
dds <- dds[rowSums(counts(dds))>=10,]  #此处的标准自己定，可以设置为1


###########-----------------------分割线------------------------------------------------------
###########数据标准化处理/归一化处理，可用于后续的热图绘制，或者基因表达量的绘图
###########数据标准化的学习连接：https://www.jianshu.com/p/8aa995149744
vsdata <- vst(dds,blind = FALSE) #对差异分析的结果进行归一化

########### PCA分析 主成分分析，用于判断肿瘤样本和正常样本是否有差异
plotPCA(vsdata,intgroup='condition') ####此处的condition是分组信息的列名

###########提取标准化后的数据，可用于后续的基因绘图
counts_vst <- as.data.frame(assay(vsdata))
save(counts_vst,file = 'counts_vst.Rda')
###########----------------------------------------分割线------------------------------------------------------

###########开始差异分析
dds <- DESeq(dds)
#save(dds,file = 'dds.Rda')
#load('dds.Rda')

############构建contrast对象  用于后续差异分析结果的提取，谁是正常样本，谁是肿瘤样本
table(sample$condition)
contrast <- c('condition','Tumor','Normal') ####此处肿瘤组一定要在前面，R包的要求
                                            ####condition为sample表格的列名，tumor和normal是这一列的元素
###########提取差异分析的结果
res <- results(dds,contrast = contrast) 

write.csv(DEG,file='DEG',row.names = F)
DEG <- as.data.frame(res)
DEG <- na.omit(DEG)        ####如果某个基因在很多样品里面的表达量为0的话，差异倍数可能为0，显示NA

####保存差异分析的结果文件
library(tidyr)
library(dplyr)
library(tibble)
DEG <- DEG %>% 
       as.data.frame %>% 
      rownames_to_column('id')
save(DEG,file = 'DEG.Rda')
