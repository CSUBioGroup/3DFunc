# 3dToFunc

## **Contents**
- [Prerequisites](#Prerequisites)
- [Prepare your inputs](#Prepare-your-inputs)
- [Parameters](#Parameters)
- [Example usage](#Example-usage)

## **Prerequisites**

[R](https://www.r-project.org/)(>=4.0.2) with library optparse, stringr, and [C++](https://gcc.gnu.org/projects/cxx-status.html)

Please run the commandline to check your dependencies
```bash
# Run the bash script with commandline
chmod a+x ./Check_dependencies.sh
./Check_dependencies.sh
```
## **Prepare your inputs**

3dToFunc requires two Hi-C data (.hic format) and two RNA-seq expression data (plain .txt format seperated by tab) as input. <br>
* For Hi-C data, we recommend using the .hic files from two comparing sampels, such as, one samples is cancer or mutant, the other sample is normal or control. There are several public data portal for downloading the .hic data: 4DN (https://data.4dnucleome.org/), ENCODE (https://www.encodeproject.org/).
* For RNA-seq expression data, we recommend using the .txt from two comparing sampels, such as, one samples is cancer or mutant, the other sample is normal or control. Here we show the format should be like,
```bash
# The first column should be Ensemble ID, the second column should be counts
ENSG00000000003	40.6903128161
ENSG00000000005	0.403147391846
ENSG00000000419	23.8598458463
ENSG00000000457	2.55279776628
ENSG00000000460	1.09035063407
```


We have integrated the cancer data from ICGC, and the normal data from GTEx, including tissues: Blood, Brain,Breast, Colon, Esophagus, Liver, Lung, Muscle, Ovary, Pancreas,Prostate, Skin, Stomach, Thyroid, Uterus. Please uncompress the integrated data before using them. The program will use the expression data of corresponding tissue when there is no specific expression data provided.
```bash
# Run the bash script with commandline
gunzip ./data/icgc_db/*
gunzip ./data/gtex_db/*
```
* For variant-gene pair, the .txt format file with 8 columns should be prepared. <br>
column #1: name of variant<br>
column #2: name of gene<br>
column #3: chromosome of variant<br>
column #4: start position of variant<br>
column #5: end position of variant<br>
column #6: chromosome of gene<br>
column #7: start position of gene<br>
column #8: end position of gene<br>
```bash
# Here is an example
chr1_100003083_G_T	ENSG00000156876	chr1	100003083	100003083	100549119	100598511	
chr1_100020068_A_C	ENSG00000156876	chr1	100020068	100020068	100549119	100598511	
chr1_100023572_T_G	ENSG00000156876	chr1	100023572	100023572	100549119	100598511	
chr1_100032823_G_T	ENSG00000156876	chr1	100032823	100032823	100549119	100598511	
chr1_100033503_T_C	ENSG00000156876	chr1	100033503	100033503	100549119	100598511
```
We have integrated the predicted v-g pairs from PCAWG, if you don't input the pairs, the program will use the integrated data. Please uncompress the integrated data before using them. 
```bash
# Run the bash script with commandline
gunzip ./data/SNP_G/*
gunzip ./data/SV_G/*
gunzip ./data/ICT_G/*
```

## **Parameters**
**--tissue** *(required), type="character", default=NULL*. Calculate 3dToFunc scores for a specific tissue. If no tissue provided, all the tissues will be calculated by default.<br><br>
**--hic_mutant** *(required), type="character", default=NULL*. The path of cancer/mutant hic matrix in .hic format. This file is required.<br><br>
**--hic_control** *(required), type="character", default=NULL*. The path of normal/control hic matrix in .hic format. This file is required.<br><br>
**--variants** *(required), type="character", default="SNP"*. The type of variants to calculate, please input SV, SNP or ICT. If no input provided, SNP will be calculated by default.<br><br>
**--expr_mutant** *(optional), type="character", default="icgc_db"*. The path of cancer/mutant RNA-seq expression, in which gene and counts were seperated by tab. If no file provided, the cancer expression data from icgc will be used. <br><br>
**--expr_control** *(optional), type="character", default="gtex_db"*. The path of normal/control RNA-seq expression, in which gene and counts were seperated by tab. If no file provided, the normal expression data from gtex will be used.<br><br>
**--pair** *(optional), type="character", default="pair_db"*. The path of variant-gene pair in .txt format, seperated by tab. If no file provided, the predicted v-g pair from PCAWG will be used.<br><br>
**--extend** *(optional), type="integer", default="5"*. The fold of length for extending variant-gene pair to calculate the flexible IF. The fold of 5 will be set by default.<br><br>
**--out_dir** *(optional), type="character", default="../"*. The fold of length for extending variant-gene pair to calculate the flexible IF. The fold of 5 will be set by default. <br> <br>


## **Example usage**

### Calculate the 3dToFunc scores for cardiovescular disease
1. Download the example mutant and control .hic file
```bash
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262966&format=file&file=GSM3262966%5FD80%5FHiC%5FRep1%2Ehic -O cardiovescular_mut_hic.hic
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3262956&format=file&file=GSM3262956%5FD00%5FHiC%5FRep1%2Ehic -O cardiovescular_ctl_hic.hic
```
2. Download the example mutant and control expression file
```bash
wget https://data.cyverse.org/dav-anon/iplant/home/litang/3dToFunc/RNA_id_D80_2_rpkm.txt -O cardiovescular_mut_exp.txt
wget https://data.cyverse.org/dav-anon/iplant/home/litang/3dToFunc/RNA_id_D0_2_rpkm.txt -O cardiovescular_ctl_exp.txt
```
3. Download the SNP-gene pair file
```bash
wget https://data.cyverse.org/dav-anon/iplant/home/litang/3dToFunc/GTEx_ge_heart_merge_2.txt -O cardiovescular_snp_gene.txt
```
4. Run the 3dToFunc.R script
```bash
cd /directory_of_3dToFunc/source/
Rscript --vanilla 3dToFunc.R --tissue heart --hic_mutant cardiovescular_mut_hic.hic  --hic_control cardiovescular_ctl_hic.hic --expr_mutant cardiovescular_mut_exp.txt --expr_control cardiovescular_ctl_exp.txt --variants SNP --pair cardiovescular_snp_gene.txt --out_dir ../
```
The output file "3dToFunc_output.txt" will be saved in the directory of 3dToFunc.
