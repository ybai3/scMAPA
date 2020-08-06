
# scMAPA

Alternative polyadenylation (APA) causes shortening or lengthening of the 3ʹ-untranslated region (3ʹ-UTR), widespread in complex tissues. To detect APA and identify cell-type-specific APA in multi-cluster setting, we developed a model-based method, scMAPA. First part of scMAPA is coded as shell scripts, which can 1) divide the aligned read data by the cell cluster; 2) convert BAM to Bedgraph file; 3) estimate the abundance of 3ʹ-UTR long and short isoform of genes in each cluster-bulk data using linear regression and quadratic programming implemented in DaPars. Second part of scMAPA is coded as a R package, which can 4) fit a logistic regression model for each gene and estimate the significance of APA; 5) Identify cluster-specific 3'UTR shortening and lengthening; 6) Do visualization to show the APA dynamics. 

The data used in this example is an intermediate analysis result on a subset of mouse cortex and midbrain dorsal data. The source of the data could be found here: [Zeisel et al. DOI:https://doi.org/10.1016/j.cell.2018.06.021] In this tutorial, we only use three cell clusters: Neurons, Immunes, and Oligos. 

## Part 1 BAM file processing and estimation of long/short isoforms

### Download BAM processing shell scripts
<!-- -->

    git clone https://github.com/ybai3/scMAPA.git
   
There are 3 shell scripts named step1.sh, step2.sh, step3.sh, corresponding to step 1), 2), and 3), respectively.

### Split BAM files

In this step, we use step1.sh to split BAM files by cell cluster information provided by user. To run step1.sh, please first load Python 3.7 with dependent modules. Step1.sh takes two positional arguments. The first argument indicates the directory to the file storing cluster information and the second argument indicates the directory to the BAM file.

Cluster information should be stored in a csv file with two columns. First column contains the cell barcodes that match the barcodes in the BAM file. Second columns contains the cluster ID. For example, the first several lines of a cluster csv file should look like: 
```
  Barcode,Cluster
  ACGTGATGAGCGTT-1,Oligos
  ACGACCCTAGACTC-1,Oligos
  ATATAGTGTCACCC-1,Oligos
  AGTCAGACGAGAGC-1,Oligos
  GATCTTTGGTCATG-1,Oligos
  AATCCTACGCATAC-1,Oligos
```
Here is an exmaple of running step1.sh:
```
module load gcc/8.2.0
module load python/bioconda-3.7-2019.03
chmod +x step1.sh
./step1.sh /Path/to/cluster.csv /Path/to/merged.bam
```
The BAM file for each cluster will be output to the working directory. 

The Python code used to split BAM file was kindly shared by Dr. Ming Tang at http://doi.org/10.5281/zenodo.3946832. It was originally designed for scATAC Seq and we modified it so that it can be used to split RNA Seq BAM files. 

### Convert BAM files to Bedgraph format

Then, we use step2.sh to convert all BAM files to bedgraph format. Please load samtools and bedtools before running step2.sh. Step2.sh takes one argument, which indicates the directory to the folder containing the BAM files. Please make sure that the input folder does not contain other BAM files, otherwise all BAM files would be converted to Bedgraph. The bedgraph files will be output to working directory. 

Here is an example of running step2.sh:
```
module load gcc/8.2.0
module load samtools/1.9
module load bedtools/2.29.0
chmod +x step2.sh
./step2.sh /Path/to/folder/containing/bam/files
```
We found that the most common publicly available BAM files from cellranger have chromosome names in the following notation: e.g. 1, 22, X, Y. To match the chromosome names in the reference bed file, we wrote command in step2.sh to automatically add "chr" in front of all chromosome names. Please note that if your BAM file's chromosome notation already include prefix 'chr', the command 'awk '{print "chr" $0}' $BEDGRAPH > $NEWBEDGRAPH' in step2.sh should be removed and the whole script should be modified accordingly. To be on safer side, you can check the output bedgraph to see if chromosome names are in right format. 

### Estimate long/short isoforms

Step3.sh will call DaPars developed by Dr. Zheng Xia in the following paper: https://doi.org/10.1038/ncomms6274. The original script of DaPars was written in Python 2.7, the version we used was slightly modified so that it is compatible with Pyhton 3.X. To run step3.sh, please load bedtools and Python with dependent modules. 

step3.sh takes one argument, which indicates the directory to the folder containning the configuration files for each cell cluster. The number of configuration files should equal to the number of cell clusters. A configuration file is a text file that clarifies the location of input files and parameters to DaPars. Here is an example of configuration file for Oligo cells:
```
Annotated_3UTR=/Path/to/mm10_refseq_extracted_3UTR.bed
Group1_Tophat_aligned_Wig=/Path/to/Cluster_Oligos.bedgraph
Group2_Tophat_aligned_Wig=/Path/to/Cluster_Oligos.bedgraph
Output_directory=PDUI_Oligos/
Output_result_file=PDUI_Oligos

#Parameters
Num_least_in_group1=1
Num_least_in_group2=1
Coverage_cutoff=30
FDR_cutoff=0.05
PDUI_cutoff=0.5
Fold_change_cutoff=0.59
```
In first line user should clarify the path to the reference bed file. In scMAPA, we provide two reference bed files, mm10_refseq_extracted_3UTR.bed for mouse and hg38_refseq_extracted_3UTR.bed for human. To generate reference bed for other species or assemblies, please refer to http://lilab.research.bcm.edu/dldcc-web/lilab/zheng/DaPars_Documentation/html/DaPars.html.

The input to the second line and third line should be same: the path to one bedgraph file generated by step2.sh. In fourth and fifth line, user can set the output directory and the name of output file.

In parameter's section, we use the default parameters from DaPars. For the meaning of each parameter, please refer to the original paper. 

Here is an example of running step3.sh:

```
module load gcc/8.2.0
module load python/anaconda3.7-5.3.1_genomics
module load bedtools/2.27.1
module load r/3.5.1

chmod +x step3.sh
./step3.sh /ihome/hpark/hyp15/mshuda_hyp15/idViralSeq_mShuda/mergeBAM/test_scMAPA/scMAPA/scMAPA_BAMprocess/configure
```
In the configure folder, we have 3 text files: Config_Immune.txt, Config_Neurons_all.txt, Config_Oligos.txt, corresponding to 3 cell clusters.

## Part 2 Model fitting, APA detection, and identification of cluster-specific 3'UTR dynamics

### Load dependencies

Please install the following R packages before using scMAPA: multcomp, lmtest, wec, Matrix.utils, nnet, stringr, gplots, ggplot2, grDevices. 
```{r, echo=TRUE, message=FALSE}
if(!require(devtools)) install.packages("devtools")
require(devtools)
if(!require(scMAPA)) devtools::install_github("ybai3/scMAPA") 
require(scMAPA)
library(multcomp)
library(lmtest)
library(wec)
library(Matrix.utils)
library(nnet)
library(stringr)
library(gplots)
library(ggplot2)
library(scMAPA)
```

### Read in PDUI matrices output from shell scripts

First, we use readin() function to input the PDUI matrices from shell scripts. Path= should be set to the folder containing (only) output files of all clusters generated by previous step. NAcutoff refers to the number of NAs among clusters should be tolerated for each gene. For example, if we have 6 clusters, genes with more than NAcutoff=3 of NAs among clusters should be filtered out. Genes with number of NAs <= NAcutoff will remain. 
```{r Load object, echo=TRUE, message=FALSE}
setwd(system.file("extdata", package = "scMAPA", mustWork = T))
a <- readin(path = getwd(), NAcutoff = 1, CPMcutoff = 10)
a[1:5,]
```

                                         Immune_long_exp Immune_short_exp Neurons_long_exp Neurons_short_exp Oligos_long_exp Oligos_short_exp
    ENSMUST00000000327.12|Clec10a|chr11|+             339                0              204                 0              27                0
    ENSMUST00000000449.8|Mkrn2|chr6|+                  21                6              154               173              31               21
    ENSMUST00000000619.7|Clcn4|chr7|-                   1               77               81              1678               9              388
    ENSMUST00000000696.6|Cd52|chr4|-                  288                0              186                 0              29                0
    ENSMUST00000000755.14|Sult5a1|chr8|-               29                0              366                 1              54                1


As shown above, readin() returns filtered isoform specific count matrix. Each cluster has two columns, corresponding to long/short isoforms. All further statistical analysis would be conducted on this matrix.

### Fit model for each gene and estimate gene-level significance of APA

In this example, to save time, we first randomly select 100 genes from 2222 genes. 
```{r}
a <- a[sample(1:nrow(a), 100),]
dim(a)
```

    [1] 2222    6


Then, we use estimateSig() function to estimate gene-level significance of APA. This function estimates significance of APA event at gene level using either model-based approach or test-based approach (argument mode = "model", "test"). Model-based approach will fit a logistic regression model for each gene and estimate significance using Likelihood ratio test. Test-based approach will condunct Fisher's exact test on all possible pairs of clusters (i.e. 1 vs. 2; 1 vs. 3; 2 vs. 3 if 3 clusters) for each gene and use FDR to adjust for multiple-testing issue. Test-based approach defines significant APA genes if they pass significance test in any pair of the clusters. We recommend to use model-based approach since all downstream analyses are based on model-based result.
```{r}
b <- estimateSig(ISOMatrix = a, mode = "model", FDR_P_cutoff = 0.05)
names(b)
```

    [1] "siglist_FDRp"  "ECoeff_Mat"    "ECoeffSig_Mat"

Model-based mode will return a list object containing three elements: 
1. siglist_FDRp contains significant APA event IDs and their gene-level FDR controlled P values.
```{r}
b$siglist_FDRp[1:5]
```

    ENSMUST00000000619.7|Clcn4|chr7|- ENSMUST00000000769.13|Serpinf1|chr11|-   ENSMUST00000000925.9|Smarcb1|chr10|-     ENSMUST00000000985.6|Oxa1l|chr14|+ 
                          4.017108e-02                           8.424570e-51                           3.912039e-31                           6.661349e-03 
       ENSMUST00000001156.7|Cfp|chrX|- 
                          8.211515e-04 
    
2. ECoeff_Mat contains estimated coefficients, P values, and SE from logsitic regression model for all genes. and 3. ECoeffSig_Mat, same matrix but only contains significant APA events. This will be used for identification of cluster-specific APA and visualization. 
```{r}
b$ECoeffSig_Mat[1:3,]
```

                                                                        Genes      Immune.coef   Neurons.coef Oligos.coef   Immune.pval   Neurons.pval  Oligos.pval Immune.se  Neurons.se   Oligos.se
    ENSMUST00000000619.7|Clcn4|chr7|-           ENSMUST00000000619.7|Clcn4|chr7|-   -1.136818     0.176077    -0.556795     2.443141e-01  1.664931e-02  0.05779045   0.9764202  0.07353817   0.2934695
    ENSMUST00000000769.13|Serpinf1|chr11|- ENSMUST00000000769.13|Serpinf1|chr11|-   -3.091410     2.155614     0.650714     4.905856e-26  4.394695e-18  0.23257649   0.2929338  0.24868416   0.5451016
    ENSMUST00000000925.9|Smarcb1|chr10|-     ENSMUST00000000925.9|Smarcb1|chr10|-  -12.564490     5.622854    -14.585091    9.866261e-01  9.851048e-01  0.98447556 749.5591201  301.17936633 749.5589398

### Identify cluster-specific 3' UTR shortening and lengthening

After identifying genes with significant APA dynamics among all clusters, we would like to go one step further to identify the clusters with significantly more or less long isoforms than the across-cluster average for the gene. We call this step as gene-cluster-level identification, it identifies clusters whose estimated regression model coefficient is significantly (FDR < 0.05 using Wald test) and strongly (absolute value of coefficient > log(2), corresponds to 2 fold change in odds ratio) deviated from 0. Both these two parameters could be changed by users.
```{r}
c <- IdentifyClusterAPA(ECoeffSig_Mat = b$ECoeffSig_Mat, FDR_P_cutoff = 0.05, CoeffCutoff = log(2))
names(c)
```

    [1] "Immune"  "Neurons" "Oligos"


As shown above, IdentifyClusterAPA() returns a list consists of tables for every cluster. Tables contain the gene IDs, APA event IDs of 3'UTR shortening or lengthening that passed filters. For example:
```{r}
c$Immune[1:5,]
```

          GeneSymbol         UTR          TranscriptID
    1       Ykt6        lengthening   ENSMUST00000002818.8
    2      Letm1        lengthening   ENSMUST00000005431.5
    3       Neu1        lengthening   ENSMUST00000007253.5
    4      Uchl5        lengthening   ENSMUST00000018333.12
    5        Hn1        lengthening   ENSMUST00000021083.6

### Visualization

To give a overview of APA dynamics, scMAPA can draw heatmap of all gene-cluster-level significant APA events by using clusterAPAheatmap function.
```{r}
clusterAPAheatmap(ECoeffSig_Mat = b$ECoeffSig_Mat, FDR_P_cutoff = 0.05, CoeffCutoff = log(2))
```
![](scMAPA_RPackage/vignettes/image/heatmap.png)

The color indicates the degree and direction of 3' UTR processing. In this example, we can observe more 3' UTR lengthening in Neurons comparing to other two cell types. 

Usually, there is a certain set of genes that may play an important role or act as marker for specific cluster. scMAPA could visualize the 3' UTR dynamic of user-defined significant APA. APAgenes parameter in APAdotplot function is not case sensitive. The input to APAgenes parameter should be a character vector consistent with second component of APA event ID (Genes column of ECoeffSig_Mat or siglist_FDRp).
```{r}
APAdotplot(ECoeffSig_Mat = b$ECoeffSig_Mat, FDR_P_cutoff = 0.05, CoeffCutoff = log(2), APAgenes = c("Rab6b","Ensa","Btf3","Hbegf","Yif1a","Fuca1","Man2b1","Paip2",
                                                                                                    "Alg5","Ehbp1l1","Tcf25","Atat1","Nfkb2","Nmral1","Mthfs","Lgals1",
                                                                                                    "Dynlt1c","Uap1l1","Snrpb","Eed","Cwc25","Pigp","St3gal6",
                                                                                                    "Gtf2a2","Eloc","Spcs2"))
```
![](scMAPA_RPackage/vignettes/image/dotplot.png)
