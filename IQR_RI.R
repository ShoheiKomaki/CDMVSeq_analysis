print("load packages")
library(data.table)

print("load dataset")
tmpMono=as.data.frame(fread("path/to/randomsampled/monocyte/DNAmdata",head=T))
tmpcd4t=as.data.frame(fread("path/to/randomsampled/CD4T/DNAmdata",head=T))

sampleAndCalcRI=function(x){
    tmpA=sample(x,n,replace=T)
    rihigh=as.numeric(as.character(quantile(tmpA,0.95)))
    rilow=as.numeric(as.character(quantile(tmpA,0.05)))
    return(rihigh-rilow)
}

# mono
tmp2=as.data.frame(t(apply(tmpMono[,7:ncol(tmpMono)],1,function(x){
    n<<-10;tmp10=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-20;tmp20=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-30;tmp30=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-40;tmp40=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-50;tmp50=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-60;tmp60=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-70;tmp70=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-80;tmp80=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-90;tmp90=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-100;tmp100=IQR(replicate(100,sampleAndCalcRI(x)))
    return(as.character(round(c(tmp10,tmp20,tmp30,tmp40,tmp50,tmp60,tmp70,tmp80,tmp90,tmp100),digits=4)))
})))
colnames(tmp2)=paste0("n",seq(from=10,to=100,by=10))
write.table(tmp2,"mono_riCalc.tsv",row.names=F,quote=F,sep="\t")

# cd4t
tmp2=as.data.frame(t(apply(tmpcd4t[,7:ncol(tmpcd4t)],1,function(x){
    n<<-10;tmp10=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-20;tmp20=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-30;tmp30=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-40;tmp40=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-50;tmp50=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-60;tmp60=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-70;tmp70=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-80;tmp80=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-90;tmp90=IQR(replicate(100,sampleAndCalcRI(x)))
    n<<-100;tmp100=IQR(replicate(100,sampleAndCalcRI(x)))
    return(as.character(round(c(tmp10,tmp20,tmp30,tmp40,tmp50,tmp60,tmp70,tmp80,tmp90,tmp100),digits=4)))
})))
colnames(tmp2)=paste0("n",seq(from=10,to=100,by=10))
write.table(tmp2,"cd4t_riCalc.tsv",row.names=F,quote=F,sep="\t")

