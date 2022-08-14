#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# function description
# Calculate flexible IF for variant-gene pairs
#===============

library(foreach)
library(doParallel)
library(stringr)

cores=detectCores()
cl <- makeCluster(cores[1]-30)    #not to overload your computer
registerDoParallel(cl)

args = commandArgs(trailingOnly=TRUE)
pair<-read.table(args[1],sep="\t")
hic_dir=args[2]
ext=args[3]
out_dir=args[4]

print("======================")
print(paste("pairs:",args[1]))
print(paste("hic_dir:",hic_dir))
print(paste("ext:",ext))
print(paste("out_dir:",out_dir))
print("======================")

straw_path="/straw"
system(paste("rm -r ",out_dir,sep=" "))
system(paste("mkdir",out_dir,sep=" "))


Calculate_bin <- function(x) {
  x<-ifelse(x/1000<5,1000,ifelse(x/1000>=5 & x/1000<10,5000,ifelse(x/1000>=10 & x/1000<50,50000,ifelse(x/1000>=50 & x/1000<100,100000,ifelse(x/1000>=100 & x/1000<500,500000,ifelse(x/1000>=500,1000000,x))))))
  x
}

finalMatrix <- foreach (count=1:nrow(pair), .combine=rbind) %dopar% {
  library(stringr)
  file=hic_dir
  
    chrom1=str_split(pair[count,1],"chr")[[1]][2]
    start1=pair[count,2]
    end1=pair[count,3]
    bin1=Calculate_bin(end1-start1)
    chrom2=str_split(pair[count,4],"chr")[[1]][2]
    start2=pair[count,5]
    end2=pair[count,6]
    bin2=Calculate_bin(end2-start2)
    bin_m<-ifelse(bin1>=bin2,as.integer(bin1),as.integer(bin2))
    print(bin_m)
    start1_ext=format(start1-as.numeric(bin1)*as.numeric(ext),scientific=FALSE)
    end1_ext=format(end1+as.numeric(bin1)*as.numeric(ext),scientific=FALSE)
    start2_ext=format(start2-as.numeric(bin2)*as.numeric(ext),scientific=FALSE)
    end2_ext=format(end2+as.numeric(bin2)*as.numeric(ext),scientific=FALSE)
    system(paste(straw_path,"KR",hic_dir,paste(chrom1,start1,end1,sep=":"),paste(chrom2,start2,end2,sep=":"),"BP",bin_m,">",paste(out_dir,paste("IF.target.tmp",count,sep="_"),sep="/"),sep=" "))
    system(paste(straw_path,"KR",hic_dir,paste(chrom1,start1_ext,end1_ext,sep=":"),paste(chrom2,start2_ext,end2_ext,sep=":"),"BP",bin_m,">",paste(out_dir,paste("IF.all.tmp",count,sep="_"),sep="/"),sep=" "))
    if(!file.size(paste(out_dir,paste("IF.target.tmp",count,sep="_"),sep="/")) == 0){
      target<-read.table(paste(out_dir,paste("IF.target.tmp",count,sep="_"),sep="/"))
      target_IF<-sum(target[,3])
      row_count_target<-nrow(target)
      bg<-read.table(paste(out_dir,paste("IF.all.tmp",count,sep="_"),sep="/"))
      bg_IF<-sum(bg[,3])
      row_count_gb<-nrow(bg)
      rate_IF<-(target_IF/row_count_target)/(bg_IF/row_count_gb)
    }
    else{
      rate_IF=0
    }
    system(paste("rm",paste(out_dir,paste("IF.target.tmp",count,sep="_"),sep="/"),sep=" "))
    system(paste("rm",paste(out_dir,paste("IF.all.tmp",count,sep="_"),sep="/"),sep=" "))
    rate_IF<-as.data.frame(rate_IF)
    rate_IF
}

stopCluster(cl)
#system(paste("echo",file,">>",paste(out_dir,paste(file,"rate_IF.txt",sep="_"),sep="/"),sep=" "))
write.table(finalMatrix,paste(out_dir,"rate_IF.txt",sep="/"),sep="\t",quote = F,row.names = F,col.names = F)
