#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# function description
# Calculate flexible IF for ICT fusion gene pairs
#===============

library("stringr")

args = commandArgs(trailingOnly=TRUE)
fusion<-read.table(args[1],sep="\t")
hic_dir=args[2]
ext=args[3]
out_dir=args[4]

print("======================")
print(paste("fusion:",args[1]))
print(paste("hic_dir:",hic_dir))
print(paste("ext:",ext))
print(paste("out_dir:",out_dir))
print("======================")

straw_path="/homeb/tangli/data/hic_file/straw"
system(paste("rm -r ",out_dir,sep=" "))
system(paste("mkdir",out_dir,sep=" "))


Calculate_bin <- function(x) {
  x<-ifelse(x/1000<5,1000,ifelse(x/1000>=5 & x/1000<10,5000,ifelse(x/1000>=10 & x/1000<50,50000,ifelse(x/1000>=50 & x/1000<100,100000,ifelse(x/1000>=100 & x/1000<500,500000,ifelse(x/1000>=500,1000000,x))))))
  x
}

for(file in list.files(hic_dir)){
  print(file)
  system(paste("echo",file,">>",paste(out_dir,paste(file,"rate_IF.txt",sep="_"),sep="/"),sep=" "))
  
for(i in 1:nrow(fusion)){
  chrom1=str_split(fusion[i,1],"chr")[[1]][2]
  start1=fusion[i,2]
  end1=fusion[i,3]
  bin1=Calculate_bin(end1-start1)
  chrom2=str_split(fusion[i,4],"chr")[[1]][2]
  start2=fusion[i,5]
  end2=fusion[i,6]
  bin2=Calculate_bin(end2-start2)
  bin_m<-ifelse(bin1>=bin2,as.integer(bin2),as.integer(bin1))
  print(bin_m)
  start1_ext=format(start1-as.numeric(bin1)*as.numeric(ext),scientific=FALSE)
  end1_ext=format(end1+as.numeric(bin1)*as.numeric(ext),scientific=FALSE)
  start2_ext=format(start2-as.numeric(bin2)*as.numeric(ext),scientific=FALSE)
  end2_ext=format(end2+as.numeric(bin2)*as.numeric(ext),scientific=FALSE)
  system(paste(straw_path,"KR",paste(hic_dir,file,sep="/"),paste(chrom1,start1,end1,sep=":"),paste(chrom2,start2,end2,sep=":"),"BP",bin_m,">",paste(out_dir,"IF.target.tmp",sep="/"),sep=" "))
  system(paste(straw_path,"KR",paste(hic_dir,file,sep="/"),paste(chrom1,start1_ext,end1_ext,sep=":"),paste(chrom2,start2_ext,end2_ext,sep=":"),"BP",bin_m,">",paste(out_dir,"IF.all.tmp",sep="/"),sep=" "))
  if(!file.size(paste(out_dir,"IF.target.tmp",sep="/")) == 0){
    target<-read.table(paste(out_dir,"IF.target.tmp",sep="/"))
    target_IF<-sum(target[,3])
    row_count_target<-nrow(target)
    bg<-read.table(paste(out_dir,"IF.all.tmp",sep="/"))
    bg_IF<-sum(bg[,3])
    row_count_gb<-nrow(bg)
    rate_IF<-(target_IF/row_count_target)/(bg_IF/row_count_gb)
  }
  else{
    rate_IF=0
  }
  system(paste("echo",rate_IF,">>",paste(out_dir,paste(file,"rate_IF.txt",sep="_"),sep="/"),sep=" "))
  }
 
}

