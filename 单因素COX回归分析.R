#############--------------------------项目实战课-单因素COX回归-----------------------------########################################
##############学习网站：https://www.jianshu.com/p/617db057df37
##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############加载数据
load('metadata.Rda')
##############加载R包
library(survival)      
##############数据调整
rt <- as.data.frame(t(expr_LIHC))
rt <- log2(rt+1)
identical(rownames(rt),rownames(metadata))   

rt <- cbind(metadata[,c(3,2)],rt) 
colnames(rt)[1:2] <- c('futime','fustat')

##############单因素Cox回归分析
#https://blog.csdn.net/weixin_43569478/article/details/108079548
#https://cloud.tencent.com/developer/article/1805497
#https://www.jianshu.com/p/8b257d1ff818
#单因素Cox回归分析筛选变量依赖于survival包的coxph函数，PH检验依赖于cox.zph函数
#通过R计算各个基因的风险比，并检验变量是否符合cox等比例风险模型的前提PH假设
#得到符合z值和P值小于0.05和符合PH假定检验（cox.zph检验P值>0.05）的基因，最终输出
#基因名，单因素COX的B值，风险比HR
#coef 回归系数
#exp(coef) HR值
#se 标准误
#Concordance C指数


#单因素cox分析
pFilter=0.01

#gene <- colnames(rt)[3]

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt) ###模型公式对因素拟合
  coxSummary = summary(cox)                              ###查看详细的数据
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP) )
  }
}

#保存单因素分析结果文件
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
surSigExp=rt[,sigGenes]
save(surSigExp,file = 'surSigExp.Rda')
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

