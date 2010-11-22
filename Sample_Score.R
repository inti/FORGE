header<-"N"
test_family="binomial"
qt="F"
pheno_name="TRAIT"
cov="F"
perm="F"
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
cat("Reading samples and phenotypes from [",fam,"]\n",sep=" ")
dat<-read.table(fam,header=T)
cat("Reading gene scores from [",scores,"]\n",sep=" ")
sc<-read.table(scores,row.names=1) # genes in the rows
cat("   '-> [",ncol(sc),"] genes read\n",sep=" ")
m<-as.matrix(sc[,2:ncol(sc)]) # genes in the row
m_var<-apply(m,1, function(x) (max(x) - min(x)))
m2<-m[which( m_var != 0),]
n_cols_removed<-abs(nrow(m) - nrow(m2))
sc<-m2
cat("Removing [",n_cols_removed,"] genes because their score have variance 0 or NA\n",sep=" ")	

cat("Making have mean 0 and sd 1\n",sep=" ")
sc<-apply(sc,1,function(row) (row - mean(row))/sd(row) ) # gene back in columns
sc<-t(sc) # genes in rows
N<-nrow(sc)
cat(N,"Genes to be analysed\n",sep=" ")
phenotype<-0
dat[dat==-9]<-NA
if (qt == "F"){
  cat("   '-> Using a binary trait\n",sep=" ")
  phenotype<-as.factor(dat$TRAIT == 1)
}
if (qt == "T"){
  cat("   '-> Using a continous trait\n",sep=" ")
  test_family<-"gaussian"
  phenotype<-dat[,pheno_name]
}
expected_coeff<-1
perm_p<-"F"
cat("Running GLM model for genes\n",sep=" ")
if (cov == "F") {
  assoc<-apply(sc,1, function(row) { list( coefficients= summary(glm(phenotype ~ row,family=test_family))$coeff)}  )
  
  if (perm != "F"){
    cat("Going to run [",perm,"] permutations\n",sep=" ")
    perm<-as.numeric(perm)
    rand_phe<-sample(phenotype)
    perm_p<-sapply(1:length(assoc), function(x)
                   {
                    if(round(x / 25) - (x / 25) == 0){
                      cat("   '-> Done with [",x,"] Genes\n",sep=" ")
                    }
                    coeff<-as.data.frame(assoc[[x]])
                    if (nrow(coeff) != expected_coeff + 1){
                     return(NA)
                    }else{
                     observed <-coeff[2,4];
                     counter<-0;
                     seen<-0;
                     for(counter in 1:perm){
                      perm_coeff<-summary(glm(rand_phe ~ as.numeric(sc[x,]),family=test_family))$coeff
                      if (perm_coeff[2,4] <= observed){
                        counter<-counter+1
                        seen<-seen+1
                        if (seen == 10){
                          seen<-(seen+1)/(counter+1)
                          return(seen)
                        }
                      }
                      rand_phe<-sample(phenotype)
                     }
                     seen<-(seen+1)/(perm+1)
                     return(seen);
                    }
                   }
                  )
  }
} else {
  cat("Using covariates in the regression analysis\n",sep=" ")
  cat("   '-> Reading covariantes from [",cov,"]\n",sep=" ")
  cov2<-read.table(cov)
  cov2<-as.matrix(cov2[,3:ncol(cov2)])
  cov2[cov2==-9]<-NA
  assoc<-apply(sc,1, function(row) { if (var(row) == 0 || is.na(var(row)) == "TRUE" ){ list (coefficients = matrix(nrow=1,ncol=2))  } else {list( coefficients= summary(glm(phenotype ~ row + cov2,family=test_family))$coeff)} } )
  expected_coeff<-ncol(cov2) + 1
    if (perm != "F"){
    cat("Going to run [",perm,"] permutations\n",sep=" ")
    perm<-as.numeric(perm)
    rand_phe<-sample(phenotype)
    perm_p<-sapply(1:length(assoc), function(x)
                   {
                    if(round(x / 5) - (x / 5) == 0){
                      cat("   '-> Done with [",x,"] Genes\n",sep=" ")
                    }
                    coeff<-as.data.frame(assoc[[x]])
                    if (nrow(coeff) != expected_coeff + 1){
                     return(NA)
                    }else{
                     observed <-coeff[2,4];
                     counter<-0;
                     seen<-0;
                     for(counter in 1:perm){
                      perm_coeff<-summary(glm(rand_phe ~ as.numeric(sc[x,]) + cov2,family=test_family))$coeff
                      if (perm_coeff[2,4] <= observed){
                        counter<-counter+1
                        seen<-seen+1
                        if (seen == 10){
                          seen<-(seen+1)/(counter+1)
                          return(seen)
                        }
                      }
                      rand_phe<-sample(phenotype)
                     }
                     seen<-(seen+1)/(perm+1)
                     return(seen);
                    }
                   }
                  )
  }
}
cat("Formatting output\n",sep=" ")
assoc<-sapply(assoc,function(g) { x<-g[[1]]; if (nrow(x) != expected_coeff + 1){ return(c(NA,NA,NA))}else{ x[2,c(1,2,4)]} } )
assoc<-t(assoc)
colnames(assoc)<-c("estimate","std_error","pvalue")
assoc<-as.data.frame(assoc)
assoc$empi_p<-perm_p
#assoc$gene_symbol<-as.vector(sc[,"V2"])
cat("Writting to file [",out,"]\n",sep=" ")
write.table(assoc,file=out,col.names=T,row.names=T,sep="\t",quote=F)
cat("Done with analysis\n",sep=" ")
