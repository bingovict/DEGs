##############--------------------------项目实战课-LASSO回归-----------------------------########################################
##############学习网站：https://www.jianshu.com/p/33425451f1f2
##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############加载数据
load('metadata.Rda')
load('surSigExp.Rda')
##############数据调整
gene <- colnames(surSigExp)[3:ncol(surSigExp)]
expr_LIHC <- expr_LIHC[gene,]

##############数据随机划分训练集和验证集
library(caret)
set.seed(12345679)
sam <- createDataPartition(metadata$event, p = .5,list = FALSE)

train <- expr_LIHC[,sam]
test <- expr_LIHC[,-sam]

train_meta <- metadata[sam,]
test_meta <- metadata[-sam,]

#查看切割的数据，分期是否均衡
table(train_meta$stage)
prop.table(table(train_meta$stage))

##############构建LASSO回归模型
##############训练集建模
library(glmnet)
x <-  t(train)
y <-  train_meta$event
set.seed(12345)
cv_fit <- cv.glmnet(x=x, y=y)
plot(cv_fit)

#系数图
fit <- glmnet(x=x, y=y)
plot(fit,xvar = "lambda")

#使用两个λ值筛选核心基因
#两条虚线代表两个λ值，一个是lambda.min,一个是lambda.1se
#使用lambda.min的准确率更高一点，基因数量更多一点。
#基因与系数存放于模型的子集beta中
model_lasso_min <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min)
#model_lasso_1se <- glmnet(x=x, y=y,lambda=cv_fit$lambda.1se)

head(model_lasso_min$beta,20)
choose_gene_min <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]

save(choose_gene_min,file = 'lasso_gene.Rda')

######用测试集验证模型的准确性
lasso.prob <- predict(cv_fit, newx=t(test), s=cv_fit$lambda.min)
re=cbind(event = test_meta$event ,lasso.prob)
re=as.data.frame(re)
colnames(re)=c('event','prob_min')
re$event=as.factor(re$event)

######绘制ROC曲线
library(pROC)
library(ggplot2)
m <- roc(test_meta$event, re$prob_min)
g <- ggroc(m,legacy.axes = T,size = 1,color = "#2fa1dd")
auc(m)

g + theme_minimal() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               colour = "grey", linetype = "dashed")+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of min = ",format(round(as.numeric(auc(m)),2),nsmall = 2)),color = "#2fa1dd")















