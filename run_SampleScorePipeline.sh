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
#CHR=$SGE_TASK_ID
CHR=$4
#TMP_DIR=/tmp
TMP_DIR=$3
FINAL_DIR=$3
FORGE_PATH=/Users/inti/DATA/My_Soft/FORGE_rep/forge.pl
ANNOTATION_FOLDER=~/DATA/Functional_Antotation/Ensemble_gene_SNP_v59

plink  --bfile $BFILE \
    --keep $IDLIST \
    --all \
    --out $TMP_DIR/$1.trainning.chr$CHR \
    --make-bed \
    --chr $CHR
        
        
plink  --bfile $BFILE \
    --remove $IDLIST \
    --all \
    --out $TMP_DIR/$1.target.chr$CHR \
    --make-bed \
    --chr $CHR
        
        
plink --bfile $TMP_DIR/$1.trainning.chr$CHR \
    --logistic \
    --out $TMP_DIR/$1.trainning.chr$CHR \
    --chr $CHR

# define output tag
OUT_TAG=$1.target.chr$CHR

perl $FORGE_PATH -bfile $TMP_DIR/$1.target.chr$CHR \
	-a $TMP_DIR/$1.trainning.chr$CHR.assoc.logistic \
	-chr $CHR \
        -m $ANNOTATION_FOLDER/ensemblv59_SNP_2_GENE.chr$CHR.txt \
	-sample_score \
	-o $TMP_DIR/$OUT_TAG 

R --no-save --no-readline \
	fam=$TMP_DIR/$OUT_TAG.sample_score.dat \
	scores=$TMP_DIR/$OUT_TAG.sample_score.mat \
	out=$TMP_DIR/$OUT_TAG.sample_score.genep  < /storage/adata/FORGE/Development/FORGE/Sample_Score.R

R --no-save --no-readline \
	fam=$TMP_DIR/$OUT_TAG.sample_score.dat \
	scores=$TMP_DIR/$OUT_TAG.sample_score.mat \
	out=$TMP_DIR/$OUT_TAG.sample_score test_family=binomial  < /storage/adata/FORGE/Development/FORGE/run_glmnet_SampleScore.R

mv $TMP_DIR/$1.trainning.chr$CHR.* $FINAL_DIR
mv $TMP_DIR/$1.target.chr$CHR.* $FINAL_DIR
mv $TMP_DIR/$OUT_TAG.* $FINAL_DIR
