library(reshape2)
library(data.table)

# set functions
eqtm.f=function(x){
 chr.pos=paste0(as.character(x[1]),":",as.character(x[2]))
 m=as.numeric(x[3:length(x)])
 model="e ~ m"
 model=as.formula(model)
 res.lm=NA
 rm(res.lm)
 tryCatch({res.lm=lm(model, na.action=na.omit)},
  error = function(e) {res.lm <<- NA})
 if (head(!is.na(res.lm))[1]) {
  if (is.na(as.character(res.lm[1]$coefficients[2]))){
   res=c(chr.pos,"NA")
  } else {
   res1=summary(res.lm)$coefficient
   res3=confint(res.lm)
   pval    =res1[2,4]
   res=c(chr.pos,pval)
  }
  } else {
   #browser()
   res=c(chr.pos,"NA")
  }
  return(res)
}

print("Run")
targets=data.frame(topCpG=c("3:98251047","5:373315","7:110733941","14:93552119","19:17000552"),gene=c("GPR15","AHRR","LRRN3","ITPK1","F2RL3"))
for(n in 1:nrow(targets)){
 top.cpg=targets[n,"topCpG"]
 gene=targets[n,"gene"]
 chr=sub(":.*$","",top.cpg)
 toppos=sub("^.*:","",top.cpg)
 print(paste0("target CpG: ",top.cpg))

 for(cell in c("cd4t","mono")){
  print(paste0("target cell: ",cell))
  # cell="cd4t"
  print("loading beta")
  fn.beta="path/to/DNAm/data"
  beta=as.data.frame(fread(fn.beta,header=T))
  colnames(beta)=sub("-",".",colnames(beta))
  beta=beta[beta$callrate >= 0.5,]
  beta=beta[!is.na(beta$pos) &beta$pos > as.numeric(toppos)-500000 & beta$pos < as.numeric(toppos)+500000,]

  # load FPKM
  print("loading fpkm")
  fpkm.fn=paste0("path/to/expression/data")
  fpkm=read.table(fpkm.fn,header=T,stringsAsFactors=F)
  gene.trackid=fpkm$tracking_id
  geneids=read.table("path/to/GeneSymbol-EnsemblID/table",header=F)
  gene.name=as.character(subset(geneids,V2=track.id)$V1)
  fpkm$tracking_id=gene.name
  colnames(fpkm)[1]="gene.name"
  fpkm=fpkm[fpkm$gene.name == gene,]
  fpkm.log=log10(fpkm[2:length(fpkm)] + 0.1)

  # load smoking status
  print("loading covariates")
  smk=read.table("path/to/covariates/data",header=T,stringsAsFactors=F)
  smk$icid=sub("-",".", smk$icid)
  rownames(smk)=smk$icid
  smk=subset(smk,smkStatus=="NS" | smkStatus=="CS")

  # id-sort for EWAS
  # beta & cov
  print("reshaping datasets")
  smk.id=rownames(smk)
  beta.id=colnames(beta)[7:ncol(beta)]
  immid.ewas=intersect(smk.id,beta.id)
  smk.ewas=smk[immid.ewas,]
  beta.ewas<-beta[,c("chr","pos",immid.ewas)]
  sex<-smk.ewas$s_sex
  sex[sex=="F"]<-0
  sex[sex=="M"]<-1
  age<-as.numeric(smk.ewas$s_age)
  smk<-smk.ewas$smkStatus
  smk[smk=="NS"]<-0
  smk[smk=="CS"]<-1
  var<-smk

  # id-sort for eqtm
  fpkm.id=colnames(fpkm.log)
  immid.eqtm=intersect(beta.id,fpkm.id)
  beta.eqtm<-beta[,c("chr","pos",immid.eqtm)]
  fpkm.log<-fpkm.log[,immid.eqtm]

  # remove gene/cpg without variation
  beta.eqtm<-beta.eqtm[apply(beta.eqtm[,3:ncol(beta.eqtm)],1,var,na.rm=TRUE) != 0,]
  fpkm.log<-fpkm.log[apply(fpkm.log,1,var,na.rm=TRUE)!= 0,]

  # load hm450 list
  hm450=read.table("path/to/HM450/CpG/list",header=T,stringsAsFactors=F)

  print("performing eqtm")
  e<-as.numeric(fpkm.log)
  res2=apply(beta.eqtm,1,eqtm.f)
  res2=as.data.frame(t(res2))
  colnames(res2)=c("chr.pos","eQTM.pval")

  print("checking LD")
  top.cpg.m=beta.eqtm[beta.eqtm$pos==toppos,c(-2:-1)]
  LD=NULL
  for(j in 1:nrow(beta.eqtm)){
  	ld=cor(as.numeric(top.cpg.m),as.numeric(beta.eqtm[j,c(-2:-1)]),use="pairwise.complete.obs")^2
	  LD=rbind(LD,ld)
  }
  LD=as.data.frame(as.numeric(LD))
  colnames(LD)="LD"
  cutpoint=c(0,0.2,0.4,0.6,0.8,1,2)
  bins=cut(LD$LD,breaks=cutpoint,right=F)
  levels(bins)=c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0","1.0")
  LD$LD.bin=bins
  LD$"chr.pos"=paste0(beta.eqtm$chr,":",beta.eqtm$pos)

  df=merge(res2,LD,by="chr.pos")

  hm450$"chr.pos"=paste0(hm450$chr,":",hm450$pos)
  df2=merge(df,hm450[,c("chr.pos","HM450")],all.x=T,all.y=F)
  df2[is.na(df2$HM450),"HM450"]="off"

  write.table(df2,paste0(cell,"_chr",chr,"_",toppos,"_",gene,"summary.tsv"),row.names=F,quote=F,sep="\t")
 }
}
