##############--------------------------项目实战课-生存资料的下载和处理-----------------------------########################################
##############设置工作目录
setwd('D:\\R\\正式科研代码')
##############读取文件
clinical <- data.table::fread('TCGA-STAD.GDC_phenotype.tsv',data.table = F)
survival <- data.table::fread('TCGA-STAD.survival.tsv',data.table = F)

##############两个数据框合并,去除生存状态缺失，生存时间缺失和小于30天的样本
colnames(clinical)[1] <- 'sample'
library(dplyr)


metadata <- survival %>% 
                          inner_join(clinical,by='sample') %>% 
                          filter(substr(sample,14,16)== '01A') %>% 
                          filter(OS.time >= 30) %>% 
                          filter(OS.time != '') %>% 
                          filter(OS != '') 

###############此处根据自己的需要获取相应的临床信息
metadata = metadata[,c(
  'sample',
  'OS',
  'OS.time',
  'age_at_initial_pathologic_diagnosis',
  'gender.demographic',
  'tumor_stage.diagnoses'
)]

###TNM分期
# 'pathologic_M',
# 'pathologic_N',
# 'pathologic_T'


colnames(metadata)=c('ID','event','time','age','gender','stage')


################生存时间改成月份
metadata$time <- metadata$time/30
range(metadata$time)

################对stage进行处理
################经典R包 stringr：https://cran.r-project.org/web/packages/stringr/index.html
library(stringr)
table(metadata$stage,useNA = 'always')

################去除分期不明的样本
metadata <-  metadata    %>% 
                               filter(stage != 'not reported') %>% 
                               filter(stage!='')

################将stage删掉，并将字母大写
metadata$stage <-  metadata$stage %>% 
                                  str_remove("stage ") %>% 
                                  str_to_upper() #小写字母转换为大写
                                  
metadata$stage <- str_remove(metadata$stage,'A|B|C')
rownames(metadata) <- metadata$ID

######################根据临床样本，过滤基因表达矩阵的样本
load('test2.Rda')


rownames(metadata)[1]
nchar(rownames(metadata)[1])

colnames(test2)[1]
nchar(colnames(test2)[1])

##########修改test2的列名
index <- substr(colnames(test2),1,16)
colnames(test2) <- index

a <- intersect(colnames(test2),rownames(metadata))
metadata <- metadata[a,]
expr_LIHC <- test2[,a]

identical(rownames(metadata),colnames(expr_LIHC))
save(metadata,expr_LIHC,file = 'metadata.Rda')
