# 3dToFunc

## Contents
- [Prerequisites](#Prerequisites)
- [Prepare your inputs](#Prepare-your-inputs)
- [Parameters](#Parameters)
- [Example usage](#Example-usage)
- [Trouble shooting](#Trouble-shooting)

## Prerequisites

[R](https://www.r-project.org/)(>=4.0.2)
[C++](https://gcc.gnu.org/projects/cxx-status.html)

Please run the commandline to check your dependencies
```bash
# Run the commandline
chmod a+x ./Check_dependencies.sh
./Check_dependencies.sh

```
## Prepare your inputs



## Parameters
--tissue, type="character", default=NULL, 
    Calculate 3dToFunc scores for a specific tissue. If no tissue provided, all the tissues will be calculated by default.
--expr_mutant, type="character", default="icgc_db", 
    The path of cancer/mutant RNA-seq expression, in which gene and counts were seperated by tab. If no file provided, the cancer expression data from icgc will be used.
--expr_control, type="character", default="gtex_db", 
    The path of normal/control RNA-seq expression, in which gene and counts were seperated by tab. If no file provided, the normal expression data from gtex will be used.
--hic_mutant, type="character", default=NULL, 
    The path of cancer/mutant hic matrix in .hic format. This file is required.
--hic_control, type="character", default=NULL, 
    The path of normal/control hic matrix in .hic format. This file is required.
--pair, type="character", default="pair_db", 
    The path of variant-gene pair in .txt format, seperated by tab. If no file provided, the predicted v-g pair from PCAWG will be used.            
--variants, type="character", default="all", 
    The type of variants to calculate, please input SV, SNP or ICT. If no input provided, SNP will be calculated by default.          
--extend, type="integer", default="5", 
    The fold of length for extending variant-gene pair to calculate the flexible IF. The fold of 5 will be set by default.
--out_dir, type="character", default="../", 
    The fold of length for extending variant-gene pair to calculate the flexible IF. The fold of 5 will be set by default.            


## Example usage

### 1. If you don't prepare any RNA-seq/Hi-C datasets,   
Please prepare unfiltered alignment file of all PETs as input, the UV Rate will be printed on the screen.
```bash
# Run the commandline
UVRate_calculation.sh align_file.sam /tmpDir

```

### 2. PC Calculation
Please prepare at least one loop files (.bedpe format) and place them in a folder, high-quality ChIP-seq/CUT&Run peak file(.bed format).
```bash
# Run the commandline
Peak-occupancy.sh /directory_of_loop_folder peak.bed prefix /outDir

```


## Trouble shooting
* Stdin error<br>
Please unzip all the annotation/loop files before running the script.

* Bedtools not found<br>
Please install bedtools and add the path into system.