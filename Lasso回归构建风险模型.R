##############--------------------------项目实战课-LASSO回归-----------------------------########################################
##############学习网站：https://www.jianshu.com/p/33425451f1f2
##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############加载数据
load('surSigExp.Rda')
rt <- surSigExp
##############加载R包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

#############对数据进行分组#############
set.seed(1234567)
inTrain<-createDataPartition(rt$fustat, p=0.5, list=F) ####切割样本，随机划分训练集和验证集
train<-rt[inTrain,]
test<-rt[-inTrain,]
trainOut <- cbind(id=row.names(train), train)
testOut <- cbind(id=row.names(test), test)

#lasso回归
x <- as.matrix(train[,c(3:ncol(train))])                ###x对应的是需要筛选的变量
y <- data.matrix(Surv(train$futime,train$fustat))       ###y对应的是，生存时间和结局
#模型构建
fit <- glmnet(x, y, family = "cox", maxit = 1000)       ####在LASSO的基础上直接做COX，可以计算风险评分
#cross-validation交叉验证，挑选lambda值 调优参数
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)    ####maxit代表模拟1000次，默认是10倍交叉验证

#使用两个λ值筛选核心基因
#两条虚线代表两个λ值，一个是lambda.min,一个是lambda.1se
#使用lambda.min的准确率更高一点，基因数量更多一点。
coef <- coef(fit, s = cvfit$lambda.min)  ###获取交叉验证误差最小的值
index <- which(coef != 0)                ###获取该基因的位置
actCoef <- coef[index]                   ###获取基因的回归系数
lassoGene <- row.names(coef)[index]      ###获取模型的基因
geneCoef <- cbind(Gene=lassoGene, Coef=actCoef) ###获取基因名及其回归系数


#ֵ输出train组的风险值
trainFinalGeneExp <- train[,lassoGene]   ####获取训练集的模型基因的表达量
#myFun=function(x){crossprod(as.numeric(x), actCoef)}
#trainScore=apply(trainFinalGeneExp, 1, myFun)
trainScore <- predict(cvfit, newx=as.matrix(train[,c(3:ncol(train))]), s="lambda.min", type="response") ###获取训练集的风险评分
outCol <- c("futime", "fustat", lassoGene)
risk <- as.vector(ifelse(trainScore>median(trainScore), "high", "low"))
train <- cbind(train[,outCol], riskScore=as.vector(trainScore), risk)
trainRiskOut <- cbind(id=rownames(train), train)

#输出test组的风险值
testFinalGeneExp=test[,lassoGene]
#testScore=apply(testFinalGeneExp, 1, myFun)
testScore <- predict(cvfit, newx=as.matrix(test[,c(3:ncol(test))]), s="lambda.min", type="response")
outCol <- c("futime", "fustat", lassoGene)
risk <- as.vector(ifelse(testScore>median(trainScore), "high", "low"))
test <- cbind(test[,outCol], riskScore=as.vector(testScore), risk)
testRiskOut <- cbind(id=rownames(test), test)

#做一个生存分析，查看训练集中中高低风险组之间是否有生存差异
diff <- survdiff(Surv(futime, fustat) ~risk, data=train)
pValue <- 1-pchisq(diff$chisq, df=1)  ###提取P值
#做一个生存分析，查看验证集中高低风险组之间是否有生存差异
diffTest <- survdiff(Surv(futime, fustat) ~risk, data=test)
pValueTest <- 1-pchisq(diffTest$chisq, df=1)

#ROC曲线下面积
predictTime <- 36    #预测时间
roc <- timeROC(T=train$futime, delta=train$fustat, ###T代表生存时间
            marker=trainScore, cause=1,            ###cause=1 不存在竞争风险，使用数字1
            weighting='aalen',
            times=c(predictTime), ROC=TRUE)        ###time 为想计算的ROC曲线的时间节点

#weighting：计算方法，默认是weighting="marginal"，KM模型
#weighting="cox" 和weighting="aalen"分别为COX模型和additive Aalen 模型。

rocTest <- timeROC(T=test$futime, delta=test$fustat, 
                marker=testScore, cause=1,
                weighting='aalen',
                times=c(predictTime), ROC=TRUE)	

###满足ROC曲线和生存分析的结果，就输出分组的结果
if((pValue<0.01) & (roc$AUC[2]>0.7) & (pValueTest<0.05) & (rocTest$AUC[2]>0.65)){
  #输出分组结果
  write.table(trainOut,file="train.data.txt",sep="\t",quote=F,row.names=F)
  write.table(testOut,file="test.data.txt",sep="\t",quote=F,row.names=F)
  #lasso图形
  pdf("lambda.pdf")
  plot(fit, xvar = "lambda", label = TRUE)
  dev.off()
  pdf("cvfit.pdf")
  plot(cvfit)
  abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
  dev.off()
  #输出模型公式和风险文件
  write.table(geneCoef, file="geneCoef.txt", sep="\t", quote=F, row.names=F)
  write.table(trainRiskOut,file="trainRisk.txt",sep="\t",quote=F,row.names=F)
  write.table(testRiskOut,file="testRisk.txt",sep="\t",quote=F,row.names=F)
  #所有样品的风险值
  allRiskOut=rbind(trainRiskOut, testRiskOut)
  write.table(allRiskOut,file="allRisk.txt",sep="\t",quote=F,row.names=F)
}




