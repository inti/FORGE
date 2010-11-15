#############
# Script to run glmnet model with gene scores results
# example run R --no-save --no-readline --slave fam=filein.sample_score.dat  scores=filein.sample_score.mat out=output.RData  < run_glmnet_SampleScore.R


test_family="binomial"
qt="F"
pheno_name="TRAIT"
cov="F"
dataout="F"
# this block parses command line option in to the program.
# an option like : option=value in transformed in a variable option with value = value
commandArgs()
for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
  } else {
    assign(ta[[1]][1],TRUE)
    cat("assigned ",ta[[1]][1]," the value of TRUE\n")
  }
}

##########################

pseudo_r_square<-function(phe,p){
# see ref 1 and 2 for Efron's method and 
# ref 2 to for Mckelvey & Zavoina's method
# References:
# 1- Efron, B. Regression and ANOVA with zero-one data: measures of residual variation. J Am Stat Assoc 1978; 73: 113-121.
# 2- Veall, M. R. & Zimmermann, K. F. Evaluating Pseudo-R 2?~@~Ys for binary probit models. Quality and Quantity 1994; 28: 151-164.
        cat("\tPseudo-R-squares are tricky. Please check these refs for more details\n")
        cat("\t\tEfron, B. Regression and ANOVA with zero-one data: measures of residual variation. J Am Stat Assoc 1978; 73: 113-121.\n")
        cat("\t\tVeall, M. R. & Zimmermann, K. F. Evaluating Pseudo-R 2?~@~Ys for binary probit models. Quality and Quantity 1994; 28: 151-164")
        efron <-function(phe,p){
                1-sum((phe - p)^2)/sum((phe - mean(phe))^2)
        }

        mckelvey_zavoina<-function ( p ){

                EV =  sum((p - mean( p ))^2)
                N = length( p )
                EV/(EV + N)
        }

        return(list (efron = efron(phe,p), mckelvey_zavoina = mckelvey_zavoina( p ) ))
}

cat("loading libraries\n")
library(glmnet)
cat("loading data\n")
cat("   '-> Reading phenotypes from [",fam,"]\n",sep=" ")
dat<-read.table(fam,header=T)
cat("   '-> Sene Scores from [",scores,"]\n",sep=" ")
sc<-read.table(scores,row.names=1)
m<-t(sc[,2:ncol(sc)])
N_predictors<-ncol(m)

cat("Removing non-needed data\n")
rm(sc)
# make phenotypes 0 or 1
cat("running cross validation and fitting lasso\n")
if (qt == "F"){
  cat("   '-> Using a binary trait\n",sep=" ")
  phenotype<-as.factor(dat$TRAIT == 1)
}
if (qt == "T"){
  cat("   '-> Using a continous trait\n",sep=" ")
  test_family<-"gaussian"
  phenotype<-dat[,pheno_name]
}
if (cov != "F"){
        cat("Using covariates in the regression analysis\n",sep=" ")
        cat("   '-> Reading covariantes from [",cov,"]\n",sep=" ")
        cov2<-read.table(cov)
        cov2<-as.matrix(cov2[,3:ncol(cov2)])
        cov2[cov2==-9]<-NA
        m<-rbind(m,cov2)
}
fit1.cv<-cv.glmnet(m,factor(phenotype),family=test_family)
#recode the phenotype
new_phe<-dat$TRAIT
new_phe[new_phe == 1]<--1
new_phe[new_phe == 2]<-1
# calculate the predicted phenotypes
cat("Predicting phenotype\n")
p<-predict(fit1.cv$glmnet.fit,m[,1:N_predictors],s=fit1.cv$lambda.min)
cat("Calculting pseudo R-square\n")
pseudoR2<-pseudo_r_square(new_phe,p)
pseudoR2<-as.data.frame(pseudoR2)
rownames(pseudoR2)<-scores
betas<-fit1.cv$glmnet.fit$beta[, which.min(fit1.cv$cvm)]
betas<-as.data.frame(betas)
out<-as.character(out)
cat("Writting gene beta values to [",out,"] Saving data image\n",sep="\n")
write.table(betas,file=paste(out,".betas",sep=""),quote=F,col.names=F,row.names=T,sep="\t")
write.table(pseudoR2,file=paste(out,".variance",sep=""),quote=F,col.names=T,row.names=T,sep="\t")
if (dataout != "F"){
        cat("Saving data image\n")
        save.image(file=paste(out,".RData",sep=""))
}
cat("Done with analysis\n",sep=" ")