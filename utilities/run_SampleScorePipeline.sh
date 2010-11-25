#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q R32.q,R8.q,R4.q,R4hi.q
#$ -r yes
#$ -t 1-26:1
#$ -N sampleScorePipeline

IDLIST=$1
BFILE=$2
CHR=$SGE_TASK_ID
#CHR=$4
#TMP_DIR=/tmp
TMP_DIR=$3
FINAL_DIR=$3
FORGE_PATH=/storage/adata/FORGE/Development/FORGE/forge.pl
ANNOTATION_FOLDER=/storage/adata/FORGE/SNPtoGeneAnnotation/Ensemble_gene_SNP_v59

echo Making PLINK Binary files for the TRAINNING set
plink  --bfile $BFILE \
    --keep $IDLIST \
    --all \
    --out $TMP_DIR/$1.trainning.chr$CHR \
    --make-bed \
    --chr $CHR \
    --silent 
echo Files written to $TMP_DIR/$1.trainning.chr$CHR.*
    
echo Making PLINK Binary files for the TARGET set        
plink  --bfile $BFILE \
    --remove $IDLIST \
    --all \
    --out $TMP_DIR/$1.target.chr$CHR \
    --make-bed \
    --chr $CHR \
    --silent     
echo Files written to $TMP_DIR/$1.target.chr$CHR.*
    
echo Running association analysis for the TRAINNING set
plink --bfile $TMP_DIR/$1.trainning.chr$CHR \
    --logistic \
    --out $TMP_DIR/$1.trainning.chr$CHR \
    --chr $CHR \
    --silent 
echo Files written to $TMP_DIR/$1.trainning.chr$CHR.*

# define output tag
OUT_TAG=$1.target.chr$CHR

perl $FORGE_PATH -bfile $TMP_DIR/$1.target.chr$CHR \
	-a $TMP_DIR/$1.trainning.chr$CHR.assoc.logistic \
	-chr $CHR \
        -m $ANNOTATION_FOLDER/ensemblv59_SNP_2_GENE.chr$CHR.txt \
	-sample_score \
	-o $TMP_DIR/$OUT_TAG 

R --no-save --no-readline --slave \
	fam=$TMP_DIR/$OUT_TAG.sample_score.dat \
	scores=$TMP_DIR/$OUT_TAG.sample_score.mat \
	out=$TMP_DIR/$OUT_TAG.sample_score.genep  < /storage/adata/FORGE/Development/FORGE/Sample_Score.R
echo Output written to $TMP_DIR/$OUT_TAG.sample_score.genep

R --no-save --no-readline --slave \
	fam=$TMP_DIR/$OUT_TAG.sample_score.dat \
	scores=$TMP_DIR/$OUT_TAG.sample_score.mat \
	out=$TMP_DIR/$OUT_TAG.sample_score \
	test_family=binomial  < /storage/adata/FORGE/Development/FORGE/run_glmnet_SampleScore.R
echo Output written to $TMP_DIR/$OUT_TAG.sample_score with exentension *.variance and *.betas

#mv $TMP_DIR/$1.trainning.chr$CHR.* $FINAL_DIR
#mv $TMP_DIR/$1.target.chr$CHR.* $FINAL_DIR
#mv $TMP_DIR/$OUT_TAG.* $FINAL_DIR
