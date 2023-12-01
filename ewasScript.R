
print("load packages")
library(ggplot2)
library(reshape2)
library(data.table)
library(corrplot)

print("load datasets")
print("loading DNAm data")
dnam=as.data.frame(fread("path/to/dna/methylation/data",head=T))

print("loading covariates data")
cov=read.table("path/to/cohort/data",head=T)

print("sort datasets") # exclude FS
ids=intersect(cov[cov$smkStatus=="NS"|cov$smkStatus=="CS","immid"],colnames(dnam)[c(-3:-1)])
dnam=dnam[,c("chrom","start","end",ids)]
cov=cov[ids,]

print("PCA")
tmp1=na.omit(dnam[,c(-3:-1)])
tmp2=tmp1[apply(tmp1,1,sd)>0,]
tmp3=prcomp(t(tmp2),scale=T)
tmp4=tmp3$x[,c(paste0("PC",1:20))]
cov=cbind(cov,tmp4)

# EWAS
se=function(x){sd(x)/sqrt(length(x))}
ewas=function(x){
 chr=as.character(x[1])
 pos=as.character(x[3])
 m=as.numeric(x[4:length(x)])
 data=cov
 data$m=m
 data=data[data$smkStatus=="NS" | data$smkStatus=="CS",]
 data$smkStatus=factor(data$smkStatus,levels=c("NS","CS"))
 cr=length(m[!is.na(m)])/length(m)
 medianNS=median(data[data$smkStatus=="NS","m"],na.rm=T)
 seNS=se(na.omit(data[data$smkStatus=="NS","m"]))
 medianCS=median(data[data$smkStatus=="CS","m"],na.rm=T)
 seCS=se(na.omit(data[data$smkStatus=="CS","m"]))
 res=c(chr,pos,cr,medianNS,seNS,medianCS,seCS)

 model=as.formula(formula)
 tryCatch({lm.res=lm(model,data=data,na.action=na.omit)},error=function(e){lm.res <<- NULL})
 if(is.null(lm.res)){
  res1=rep(NA,5)
 }else if (!is.na(lm.res$coefficients["smkStatusCS"])){
  # get stats
  coef=summary(lm.res)$coefficients[2,1]
  coefSE=summary(lm.res)$coefficients[2,2]
  coefLow=confint(lm.res)["smkStatusCS",1]
  coefHigh=confint(lm.res)["smkStatusCS",2]
  pval=summary(lm.res)$coefficients[2,4]
  res1=c(coef,coefSE,coefLow,coefHigh,pval)
 }else{
  res1=rep(NA,5)
 }
 # get adjusted DNAm
 model=as.formula(sub("smkStatus+","",formula))
 tryCatch({lm.res=lm(model,data=data,na.action=na.omit)},error=function(e){lm.res <<- NULL})
 if(is.null(lm.res)){
  res2=rep(NA,4)
 }else if (!is.na(lm.res$coefficients["s_age"])){
  data=data[!is.na(data$m),]
  data$adjm=as.numeric(lm.res$residuals)
  medianNS=median(data[data$smkStatus=="NS","adjm"])
  seNS=se(na.omit(data[data$smkStatus=="NS","adjm"]))
  medianCS=median(data[data$smkStatus=="CS","adjm"])
  seCS=se(na.omit(data[data$smkStatus=="CS","adjm"]))
  res2=c(medianNS,seNS,medianCS,seCS)
 }else{
  res2=rep(NA,4)
 }
 return(c(res,res1,res2))
}

formula=paste0("m~smkStatus+s_age+s_sex+",paste("PC",1:20,sep="",collapse="+"))
print(paste0("formula: ",formula))
tmpA=as.data.frame(t(apply(dnam,1,ewas)))
colnames(tmpA)=c("chr","pos","cr","medianNS","seNS","medianCS","seCS",
 "coef","coefSE","coefLow","coefHigh","pval",
 "medianAdjNS","seAdjNS","medianAdjCS","seAdjCS")
write.table(tmpA,paste0("ewasRes_",sub("m~","",formula),".tsv"),row.names=F,quote=F,sep="\t")
