##############家人们投一下币，点一下关注吧！！
##############QQ群：1121444037  为了防止有人进来打广告，请不要使用小号加群！！！管理员可能不会给你通过！
##############由于QQ群开放讨论，请注意防骗。
##############2022-11-24
##############所有课程都是紧密联系的，如果没看过之前的课程请不要观看这节课！！！！
##############以后的课用到的R语言技能都是基础篇和中级篇学过的！没有简单的学过R语言只想套代码的后，可以放弃了
##############出课程的目的一是给大家快速入门R语言，二用来纪录我自己的学习内容。
##############-------------------------分割线------------------------------------------############
##############-------------------------项目实战课————数据清洗########################################
##############设置工作目录
setwd('D:\\R\\正式科研代码')
######加载数据
load('TCGA_STAD.Rda')
#####加载R包
library(dplyr)
#####基因去重复，保留基因表达量最大的重复基因
STAD_counts <- mydata_counts %>% 
               mutate(mean=rowMeans(.[,-1])) %>% 
               arrange(desc(mean)) %>% 
               distinct(ID,.keep_all = T) %>% 
               select(-mean)

STAD_fpkm <- mydata_fpkm %>% 
             mutate(mean=rowMeans(.[,-1])) %>% 
             arrange(desc(mean)) %>% 
             distinct(ID,.keep_all = T) %>% 
             select(-mean)

STAD_tpm <- mydata_tpm %>% 
  mutate(mean=rowMeans(.[,-1])) %>% 
  arrange(desc(mean)) %>% 
  distinct(ID,.keep_all = T) %>% 
  select(-mean)
######样本排序
rownames(STAD_counts) <- STAD_counts$ID
table(substr(colnames(STAD_counts),14,16))

Tumor <- grep('01A',colnames(STAD_counts))
Tumor_expr <- STAD_counts[,Tumor]

Normal <- grepl('11A',colnames(STAD_counts))
Normal_expr <- STAD_counts[,Normal]

identical(rownames(Normal_expr),rownames(Tumor_expr))
symbol_counts <- cbind(Normal_expr,Tumor_expr)

######同样的处理fpkm数据
rownames(STAD_fpkm) <- STAD_fpkm$ID
table(substr(colnames(STAD_fpkm),14,16))

tumor <- grep('01A',colnames(STAD_fpkm))
tumor_expr <- STAD_fpkm[,tumor]

normal <- grepl('11A',colnames(STAD_fpkm))
normal_expr <- STAD_fpkm[,normal]

identical(rownames(normal_expr),rownames(tumor_expr))
symbol_fpkm <- cbind(normal_expr,tumor_expr)

######处理tpm数据
rownames(STAD_tpm) <- STAD_tpm$ID
table(substr(colnames(STAD_tpm),14,16))

tumor <- grep('01A',colnames(STAD_tpm))
tumor_expr <- STAD_tpm[,tumor]

normal <- grepl('11A',colnames(STAD_tpm))
normal_expr <- STAD_tpm[,normal]

identical(rownames(normal_expr),rownames(tumor_expr))
symbol_tpm <- cbind(normal_expr,tumor_expr)

save(symbol_counts,symbol_fpkm,symbol_tpm,file = 'Symbol.Rda')
