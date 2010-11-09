header<-"N"
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
dat<-read.table(fam,header=T)
sc<-read.table(scores,row.names=1)

N<-nrow(sc)
cat(N,"Genes to be analysed\n",sep=" ")

assoc<-apply(sc[,2:ncol(sc)],1, function(row) { list( coefficients= summary(glm(dat$TRAIT == 1 ~ row))$coeff) } )
assoc<-sapply(assoc,function(g) { x<-g[[1]]; if (nrow(x) == 1){ return(c(NA,NA,NA))}else{ x[2,c(1,2,4)]} } )
assoc<-t(assoc)
colnames(assoc)<-c("estimate","std_error","pvalue")
assoc<-as.data.frame(assoc)
assoc$gene_symbol<-as.vector(sc[,"V2"])
write.table(assoc,file=out,col.names=T,row.names=T,sep="\t",quote=F)
