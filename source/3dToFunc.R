#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# function description
# 3dToFunc for variant-gene pairs 
#===============

# function for installing needed packages
installpkg <- function(x){
    if(x %in% rownames(installed.packages())==FALSE) {
        if(x %in% rownames(available.packages())==FALSE) {
            paste(x,"is not a valid package - please check again...")
        } else {
            install.packages(x)           
        }

    } else {
        paste(x,"package already installed...")
    }
}

# install necessary packages
required_packages  <- c("optparse","stringr")
lapply(required_packages,installpkg)

getExtension <- function(file){ 
    ex <- strsplit(basename(file), split="\\.")[[1]]
    return(ex[-1])
}

Gtex_tissue<-c("Blood","Brain","Breast","Colon","Esophagus","Liver","Lung","Muscle","Ovary","Pancreas","Prostate","Skin","Stomach","Thyroid","Uterus")
Tissue_match<-function(tissue){
  for(t in Gtex_tissue){
    if(grepl(t,tissue,ignore.case=TRUE)){
      return(t)
      break
    }
  }
} 

Calculate_bin <- function(x) {
  x<-ifelse(x/1000<5,1000,ifelse(x/1000>=5 & x/1000<10,5000,ifelse(x/1000>=10 & x/1000<50,50000,ifelse(x/1000>=50 & x/1000<100,100000,ifelse(x/1000>=100 & x/1000<500,500000,ifelse(x/1000>=500,1000000,x))))))
  x
}

Range<- function(x){(x-min(x))/(max(x)-min(x))}

## Parse inputs
## ============================================================
option_list = list(
  make_option(c("--tissue"), type="character", default=NULL, 
              help="Calculate 3dToFunc scores for a specific tissue. If no tissue provided, all the tissues will be calculated by default.", metavar="character"),
  make_option(c("--expr_mutant"), type="character", default="icgc_db", 
              help="The path of cancer/mutant RNA-seq expression, in which gene and counts were seperated by tab. If no file provided, the cancer expression data from icgc will be used.", metavar="character"),
  make_option(c("--expr_control"), type="character", default="gtex_db", 
              help="The path of normal/control RNA-seq expression, in which gene and counts were seperated by tab. If no file provided, the normal expression data from gtex will be used.", metavar="character"),
  make_option(c("--hic_mutant"), type="character", default=NULL, 
              help="The path of cancer/mutant hic matrix in .hic format. This file is required.", metavar="character"),
  make_option(c("--hic_control"), type="character", default=NULL, 
              help="The path of normal/control hic matrix in .hic format. This file is required.", metavar="character"),
  make_option(c("--pair"), type="character", default="pair_db", 
              help="The path of variant-gene pair in .txt format, seperated by tab. If no file provided, the predicted v-g pair from PCAWG will be used.", metavar="character"),             
  make_option(c("--variants"), type="character", default="all", 
              help="The type of variants to calculate, please input SV, SNP or ICT. If no input provided, SNP will be calculated by default.", metavar="character"),           
  make_option(c("--extend"), type="integer", default="5", 
              help="The fold of length for extending variant-gene pair to calculate the flexible IF. The fold of 5 will be set by default.", metavar="integer"),
  make_option(c("--out_dir"), type="character", default="../", 
              help="The fold of length for extending variant-gene pair to calculate the flexible IF. The fold of 5 will be set by default.", metavar="character")             

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$hic_mutant) | is.null(opt$hic_control)){
  print_help(opt_parser)
  stop("At least two argumens must be supplied (.hic file).n", call.=FALSE)
};

if(is.null(opt$tissue)){
  print_help(opt_parser)
  stop("Please input the name of tissue for calculation.", call.=FALSE)
}else if(is.null(Tissue_match(opt$tissue)) & opt$pair=="pair_db"){
  stop("Please input a correct tissue name (Blood,Brain,Breast,Colon,Esophagus,Liver,Lung,Muscle,Ovary,Pancreas,Prostate,Skin,Stomach,Thyroid,Uterus).", call.=FALSE)
}else if(!is.null(Tissue_match(opt$tissue)) & opt$pair=="pair_db"){
    tissue<-Tissue_match(opt$tissue)
}else{
    tissue<-opt$tissue
};

if(getExtension(opt$hic_mutant)!="hic" | getExtension(opt$hic_control)!="hic"){
  stop("Please check the format of hic file.", call.=FALSE)
};

if(opt$expr_mutant=="icgc_db"){
    print("cancer/mutant RNA-seq expression: icgc cancer data")
    print("Reading icgc db, please wait...")
    RNA_fpkm<-read.table(paste("../data/icgc_db/icgc_db_",tissue,".txt",sep=""),sep="\t",header = T)
}else if(opt$expr_mutant!="icgc_db" & file.exists(opt$expr_mutant)){
    print(paste("cancer/mutant RNA-seq expression:",opt$expr_mutant))
    RNA_fpkm<-read.table(opt$expr_mutant)
}else{
    stop(paste(opt$expr_mutant,"dosen't exist."), call.=FALSE)
};

if(opt$expr_control=="gtex_db"){
    print("normal/control RNA-seq expression: gtex data")
    print("Reading gtex db, please wait...")
    GTEx_fpkm<-read.table(paste("../data/gtex_db/gtex_db_",tissue,".txt",sep=""),sep="\t",header = T)
}else if(opt$expr_control!="gtex_db" & file.exists(opt$expr_control)){
    print(paste("normal/control RNA-seq expression:",opt$expr_control))
    GTEx_fpkm<-read.table(opt$expr_control)
}else{
    stop(paste(opt$expr_control,"dosen't exist."), call.=FALSE)
};

if(file.exists(opt$hic_mutant)){
    print(paste("cancer/mutant hic matrix:",opt$hic_mutant))
}else{
    stop(paste(opt$hic_mutant,"dosen't exist."), call.=FALSE)
};

if(file.exists(opt$hic_control)){
    print(paste("normal/control hic matrix.",opt$hic_control))
}else{
    stop(paste(opt$hic_control,"dosen't exist."), call.=FALSE)
};

if(opt$variants!="SNP" & opt$variants!="SV"  & opt$variants!="ICT"){
    stop(paste("Please input the correct type of variants (SNP, SV, or ICT)."), call.=FALSE)
}else{
    variants<-paste(opt$variants,"G",sep="_")
    print(paste("The type of variants to calculate:",variants))
}

if(opt$pair=="pair_db"){
    print("V-G pair: predicted v-g pairs from PCAWG.")
    for(i in list.files(paste("../data/",variants,sep="/"))){
        if(grepl(tissue,i,ignore.case=TRUE)){
            pair<-read.table(paste("../data",variants,i,sep="/"))
            break;
        }
    }
}else if(opt$pair!="pair_db" & file.exists(opt$pair)){
    print(paste("V-G pair:",opt$pair))
    pair<-read.table(opt$pair)
}else{
    stop(paste(opt$pair,"dosen't exist."), call.=FALSE)
};

if(is.numeric(opt$extend)){
    print(paste("Extension fold:",opt$extend))
}else{
    stop(paste("The extension fold should be an integer from 1-10."), call.=FALSE)
}

if(opt$out_dir=="../"){
    print("Output directory: ../")
}else{
    system(paste("mkdir",opt$out_dir,sep=" "))
}

## Calculate expression change
## ============================================================
if(opt$expr_mutant!="icgc_db"){
        sd<-sd(RNA_fpkm[,2])
        RNA_fpkm[,3]<-RNA_fpkm[,2]+sd/nrow(RNA_fpkm)
        RNA_fpkm[,4]<-RNA_fpkm[,2]-sd/nrow(RNA_fpkm)
        row.names(RNA_fpkm)<-RNA_fpkm[,1]
        RNA_fpkm<-RNA_fpkm[,-1]
    }
if(opt$expr_control!="gtex_db"){
        sd2<-sd(GTEx_fpkm[,2])
        GTEx_fpkm[,3]<-GTEx_fpkm[,2]+sd2/nrow(GTEx_fpkm)
        GTEx_fpkm[,4]<-GTEx_fpkm[,2]-sd2/nrow(GTEx_fpkm)
        row.names(GTEx_fpkm)<-GTEx_fpkm[,1]
        GTEx_fpkm<-GTEx_fpkm[,-1]
    }
if(opt$variants=="SNP" | opt$variants=="SV"){
    for(i in 1:nrow(pair)){
        cancer_expr<-t(na.omit(RNA_fpkm[pair[i,2],]))
        normal_expr<-t(na.omit(GTEx_fpkm[pair[i,2],]))
        if(length(cancer_expr)>1 & length(normal_expr)>1){
            p<-t.test(cancer_expr,normal_expr)
            pair[i,"expr_change"]<-p$p.value
        }else{
            pair[i,"expr_change"]<-"NaN"
        }
    }
}else{
    for(i in 1:nrow(pair)){
        cancer_expr_1<-t(RNA_fpkm[pair[i,1],])
        normal_expr_1<-t(GTEx_fpkm[pair[i,1],])
        p_1<-t.test(cancer_expr_1,normal_expr_1)
        cancer_expr_2<-t(RNA_fpkm[pair[i,2],])
        normal_expr_2<-t(GTEx_fpkm[pair[i,2],])
        p_2<-t.test(cancer_expr_2,normal_expr_2)
        pair[i,"expr_change"]<-(p_1$p.value+p_1$p_2.value)/2
    }
}

## Calculate flexible IF
## ============================================================
ext<-opt$extend
system(paste("mkdir","./tmp",sep=" "))
tmp_dir<-"./tmp"

for(i in 1:nrow(pair)){
  #chrom1=str_split(pair[i,3],"chr")[[1]][2]
  chrom1=pair[i,3]
  start1=ifelse((pair[i,5]-pair[i,4])>5000, pair[i,4],ifelse(pair[i,4]>2500, pair[i,4]-2500,0))
  end1=ifelse((pair[i,5]-pair[i,4])>5000, pair[i,5],pair[i,5]+2500)
  bin1=Calculate_bin(end1-start1)
  #chrom2=str_split(pair[i,6],"chr")[[1]][2]
  chrom2=pair[i,6]
  start2=ifelse((pair[i,8]-pair[i,7])>5000, pair[i,7],ifelse(pair[i,7]>2500, pair[i,7]-2500,0))
  end2=ifelse((pair[i,8]-pair[i,7])>5000, pair[i,8],pair[i,8]+2500)
  bin2=Calculate_bin(end2-start2)
  bin_m<-ifelse(bin1>=bin2,as.integer(bin2),as.integer(bin1))
  print(bin_m)
  start1_ext=format(start1-as.numeric(bin1)*as.numeric(ext),scientific=FALSE)
  end1_ext=format(end1+as.numeric(bin1)*as.numeric(ext),scientific=FALSE)
  start2_ext=format(start2-as.numeric(bin2)*as.numeric(ext),scientific=FALSE)
  end2_ext=format(end2+as.numeric(bin2)*as.numeric(ext),scientific=FALSE)
  start1_ext<-ifelse(start1_ext<0,0,start1_ext)
  start2_ext<-ifelse(start2_ext<0,0,start2_ext)

  # flexible IF for cancer Hi-C
  system(paste("./C++/straw","KR",opt$hic_mutant,paste(chrom1,start1,end1,sep=":"),paste(chrom2,start2,end2,sep=":"),"BP",bin_m,">",paste(tmp_dir,"IF.target.tmp",sep="/"),sep=" "))
  system(paste("./C++/straw","KR",opt$hic_mutant,paste(chrom1,start1_ext,end1_ext,sep=":"),paste(chrom2,start2_ext,end2_ext,sep=":"),"BP",bin_m,">",paste(tmp_dir,"IF.all.tmp",sep="/"),sep=" "))
  
  if(!file.size(paste(tmp_dir,"IF.target.tmp",sep="/")) == 0){
    target<-read.table(paste(tmp_dir,"IF.target.tmp",sep="/"))
    sd1<-ifelse(nrow(target)==1,0.05,sd(target[,3]))
    target_IF<-sum(target[,3])
    row_count_target<-nrow(target)
    bg<-read.table(paste(tmp_dir,"IF.all.tmp",sep="/"))
    bg_IF<-sum(bg[,3])
    sd2<-ifelse(nrow(bg)==1,0.01,sd(bg[,3]))
    row_count_gb<-nrow(bg)

    rate_IF_cancer_1<-(target_IF/row_count_target)/(bg_IF/row_count_gb)
    rate_IF_cancer_2<-((target_IF-sd1)/row_count_target)/((bg_IF-sd2)/row_count_gb)
    rate_IF_cancer_3<-((target_IF+sd1)/row_count_target)/((bg_IF+sd2)/row_count_gb)
    rate_IF_cancer<-c(rate_IF_cancer_1,rate_IF_cancer_2,rate_IF_cancer_3)
  }
  else{
    rate_IF_cancer=c(0,0,0)
  }

  # flexible IF for normal Hi-C
  system(paste("./C++/straw","KR",opt$hic_control,paste(chrom1,start1,end1,sep=":"),paste(chrom2,start2,end2,sep=":"),"BP",bin_m,">",paste(tmp_dir,"IF.target.tmp",sep="/"),sep=" "))
  system(paste("./C++/straw","KR",opt$hic_control,paste(chrom1,start1_ext,end1_ext,sep=":"),paste(chrom2,start2_ext,end2_ext,sep=":"),"BP",bin_m,">",paste(tmp_dir,"IF.all.tmp",sep="/"),sep=" "))
  
  if(!file.size(paste(tmp_dir,"IF.target.tmp",sep="/")) == 0){
    target<-read.table(paste(tmp_dir,"IF.target.tmp",sep="/"))
    sd1<-ifelse(nrow(target)==1,0.05,sd(target[,3]))
    target_IF<-sum(target[,3])
    row_count_target<-nrow(target)
    bg<-read.table(paste(tmp_dir,"IF.all.tmp",sep="/"))
    bg_IF<-sum(bg[,3])
    row_count_gb<-nrow(bg)
    sd2<-ifelse(nrow(bg)==1,0.01,sd(bg[,3]))

    rate_IF_normal_1<-(target_IF/row_count_target)/(bg_IF/row_count_gb)
    rate_IF_normal_2<-((target_IF-sd1)/row_count_target)/((bg_IF-sd2)/row_count_gb)
    rate_IF_normal_3<-((target_IF+sd1)/row_count_target)/((bg_IF+sd2)/row_count_gb)
    rate_IF_normal<-c(rate_IF_normal_1,rate_IF_normal_2,rate_IF_normal_3)
  }else{
    rate_IF_normal=c(0,0,0)
  }
    p<-t.test(rate_IF_cancer,rate_IF_normal)
    pair[i,"IF_change"]<-p$p.value

    print(c(i,p$p.value))

  }

## Nonlinear Least Square Curve Fitting
## ============================================================
eDecay <- function(t, ampl, tau) (ampl*exp(-t/tau))

pair<-pair[pair$expr_change!="NaN",]
pair<-pair[pair$IF_change!="NaN",]
pair<-na.omit(pair)
pair$expr_change<-as.numeric(pair$expr_change)
pair$IF_change<-as.numeric(pair$IF_change)
pair$pred1<-eDecay(pair$IF_change, 151.17,38.43)
pair$pred2<-eDecay(pair$IF_change, 205.74,110.93)
pair$diff<-abs(log2(pair$expr_change)+log2(pair$IF_change))
pair$diff2<-Range(pair$diff)

for(i in 1:nrow(pair)){
  p<-t.test(pair$diff,mu=pair[i,"diff"],alternative="less")
  pair[i,"p.value"]=p$p.value
}

write.table(pair,"~/Downloads/test_3.txt",sep="\t")