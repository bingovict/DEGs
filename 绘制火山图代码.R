##############绘制火山图的代码链接：科研猫：https://cloud.tencent.com/developer/article/1512442
##############设置工作目录
setwd('D:\\R\\正式科研代码')
#############加载数据
load('DEG.Rda')
#write.csv(DEG,'DEG.csv',row.names = F)

#############第一列基因ID，
#############第二列basemean：该基因在所有样本中的均值
#############第三列log2FoldChange 为肿瘤比上正常（实验组比上对照组）的差异，大于0表示上调，小于0表示下调
#############第四列不知道
#############第五列不知道
#############第六列和第七列P值和校正P值
#############加载R包
#install.packages('ggpubr')
#install.packages('ggthemes')

library(ggpubr)
library(ggthemes)

#############修改列名
deg.data <-DEG 
colnames(deg.data)[3] <- 'logFC'
colnames(deg.data)[7] <- 'adj.P.Val'

#############对差异基因校正后的P值进行-log10转换
deg.data$logP <- -log10(deg.data$adj.P.Val)

############绘制最基本的火山图
ggscatter(deg.data,x= 'logFC',y = 'logP')

ggscatter(deg.data,x= 'logFC',y = 'logP') + theme_base() #theme_base() 加方框

############新增一列，区分上调和下调的基因
deg.data$Group <- 'Not-significant'
############logFC大于2（4倍差异) 大于1 （2倍差异）大于0.75 1倍差异

############which函数用法，针对的是向量，返回位置
#x <- c(1,1,1,0,0,1,1)  
#which(x!=1)
#x[which(x!=1)] <- '2'

deg.data$Group[which((deg.data$adj.P.Val) < 0.05 & (deg.data$logFC > 1) )] <- 'Up'
deg.data$Group[which((deg.data$adj.P.Val) < 0.05 & (deg.data$logFC < -1) )] <- 'Down'

table(deg.data$Group)   ####查看上调下调的基因

############绘制新的火山图
ggscatter(deg.data,x = 'logFC',y = 'logP',color = 'Group')

############改变颜色和点的大小
ggscatter(deg.data,x = 'logFC',y = 'logP',color = 'Group',
          palette = c('green','gray','red'), size = 1   )

###########加上线条
ggscatter(deg.data,x = 'logFC',y = 'logP',color = 'Group',
          palette = c('green','gray','red'), size = 1   )+
          geom_hline(yintercept = 1.3)+     #####-log10(0.05)=1.3 添加水平线
          geom_vline(xintercept = c(-1,1))        

###########变成虚线
ggscatter(deg.data,x = 'logFC',y = 'logP',color = 'Group',
          palette = c('green','gray','red'), size = 1   )+
          geom_hline(yintercept = 1.3, linetype = 'dashed')+     #####-log10(0.05)=1.3 添加水平线
          geom_vline(xintercept = c(-1,1), linetype = 'dashed')  


#########将上下调差异倍数前10的基因绘制在图中
####新增一列label
deg.data$Lable <- ''
####将差异表达的P值进行从小到大排序
deg.data <- deg.data[order(deg.data$adj.P.Val),]
####上调的基因中选取P值最小的10个
up.genes <- head(deg.data$id[which(deg.data$Group == 'Up')],10)
####下调的基因中选取P值最小的前10个
down.genes <- head(deg.data$id[which(deg.data$Group == 'Down')],10)

######将上下调的基因合并加入到Label中
deg.top10.genes <- c(as.character(up.genes),as.character(down.genes))

deg.data$Lable[match(deg.top10.genes,deg.data$id)] <- deg.top10.genes

#x <- c("A","B","C")
#y <- c("B","D","E","A","C")
#match(x,y) 
#y[match(x,y)] <- '2'

library(ggplot2)
#####改变坐标轴的颜色和坐标轴的标注，使得图片更美观
ggscatter(deg.data,x = 'logFC',y = 'logP',color = 'Group',
          palette = c('#2f5688','#BBBBBB','#CC0000'), size = 1,
          label = deg.data$Lable,
          font.label = 12,
          repel = T,
          xlab = 'log2FoldChange',
          ylab = '-log10(Adjust P-value)')+
  geom_hline(yintercept = 1.30, linetype = 'dashed')+    
  geom_vline(xintercept = c(-1,1), linetype = 'dashed')+
  theme(legend.text = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4)))


pdf(file = 'my_volcano.pdf')
ggscatter(deg.data,x = 'logFC',y = 'logP',color = 'Group',
          palette = c('#2f5688','#BBBBBB','#CC0000'), size = 1,
          label = deg.data$Lable,
          font.label = 8,
          repel = T,
          xlab = 'log2FoldChange',
          ylab = '-log10(Adjust P-value)')+
  geom_hline(yintercept = 1.30, linetype = 'dashed')+    
  geom_vline(xintercept = c(-2,2), linetype = 'dashed')+
  theme(legend.text = element_text(size = 12))
dev.off()
