##############学习网站：https://www.jianshu.com/p/2214a9e51837
##############          https://www.jianshu.com/p/9b53cf2409b2
##############设置工作目录
setwd('D:\\R享\\正式科研代码')
##############加载数据
load('metadata.Rda')
##############          任意指标的生存分析，km-plot
library(survival)       ####没安装该包自己安装
library(survminer)
library(ggplot2)
#KM-plot是一种单因素的生存分析，可用于研究一个因素对于生存时间的影响
#同时，Kaplan-Meier方法只能针对分类变量（治疗A vs 治疗B，男 vs 女）
#不能分析连续变量对生存造成的影响。

############# 性别作为分类变量做生存分析
fit <- survfit(Surv(time, event)~gender, data=metadata)  
ggsurvplot(fit,pval=TRUE)

##############这句代码是直接网上检索来的
ggsurvplot(
  fit,
  data = metadata,
  size = 1,                 # 线条大小
  palette =
    c("#E7B800", "#2E9FDF"),# 分组颜色
  conf.int = TRUE,          # 是否添加置信区间
  pval = TRUE,              # 是否添加p值
  risk.table = TRUE,        # 是否添加风险表
  risk.table.col = "strata",# 分线表颜色
  legend.labs =
    c("Male", "Female"),    # 图例标签
  risk.table.height = 0.25, # 生存曲线图下所有生存表的高度，数值0-1之间
  ggtheme = theme_bw()      # 主题
)


##############连续变量（指标）做生存分析
########ifelse函数用法
#Usage
#ifelse(test, yes, no)
#x <- c(1,4,8,9,5,4,2,1,4,5)
#ifelse(x<5,'T','F')

#median(metadata$age,na.rm = T)
#group <-  ifelse(metadata$age>60,">60","<=60")
#table(group)
#metadata <- cbind(metadata,group)

#fit2 <- survfit(Surv(time, event)~group, data=metadata)
#ggsurvplot(fit2,pval =TRUE)
#ggsurvplot(fit2,pval =TRUE, data = metadata, risk.table = TRUE)


##############任意基因做生存分析
a <- 'TMOD1'
metadata$Gene <- ifelse(as.numeric(expr_LIHC[a,]) > median(as.numeric(expr_LIHC[a,])) ,'High','Low')
fit3 <- survfit(Surv(time, event)~Gene, data=metadata)
ggsurvplot(fit3,pval =TRUE, 
           data = metadata, risk.table = TRUE,
           legend.title = a, 
           legend.labs = c("High", "Low"))
  


#############非常强大的R包，tinyarray,生信技能树的老师：小洁老师
#############简书：小洁忘了怎么分身：https://www.jianshu.com/p/f31a5ff62474
#############如果使用该R包发文章了，请在文章中对生信技能树进行致谢
#############生信技能树团队在网上写了大量的教学帖子和教程，心中只有感谢！
#if(!require(devtools))install.packages("devtools")
#if(!require(tinyarray))devtools::install_github("xjsun1221/tinyarray",upgrade = F)

library(tinyarray)
library(patchwork)
UP <- c("ALB","WT1","SPP1","KRT7","CDKN2A","PRAME")
DOWN <- c("APOA1","SPINK5","ALDOB","FLG","GCG","TMOD1")
myplot <-  exp_surv(exprSet_hub = expr_LIHC[c("ATP12A","APOB"),],
              meta = metadata)
wrap_plots(myplot)
