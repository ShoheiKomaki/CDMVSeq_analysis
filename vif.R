print("load packages")
library(ggplot2)
library(reshape2)
library(data.table)
library(corrplot)
library(car)

# vif
cov=read.table("path/to/covariate/data",head=T,sep="\t")
cov$smkStatus2=as.numeric(as.factor(cov$smkStatus))
cov$s_sex2=as.numeric(as.factor(cov$s_sex))
tmp1=NULL
formula=paste0("smkStatus2~s_age+s_sex2+",paste("PC",1:20,sep="",collapse="+"))
cov$smkStatus2=as.numeric(as.factor(cov$smkStatus))
tmpA=vif(lm(formula,data=cov))
tmpB=data.frame(vif=tmpA)
tmpB$response=x
tmpB$predictor=rownames(tmpB)
tmp1=rbind(tmp1,tmpB)
tmp2=as.data.frame(tmp1[,c("response","predictor","vif")])
write.table(tmp2,"vif.tsv",row.names=F,quote=F,sep="\t")
