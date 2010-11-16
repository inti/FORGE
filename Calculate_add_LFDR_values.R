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
#######


library(locfdr)
if (header == "T"){
	data<-read.table(filein,sep = "\t",skip = 1) 
} else {
	data<-read.table(filein,sep = "\t") 
}
# recalculate the FORGE p-values. This because usually does not work on double precission
data[,9]<-pchisq(data[,10],df= data[,11],lower.tail=F); 
# recalculate the Galwey p-values
data[,12]<-pchisq(data[,13],df=2* data[,14],lower.tail=F); 
# re-order that data
data<-data[order(data[,9]),]
data$new_sidak<-apply(data,1, function(x) { if (as.numeric(x[8]) == "0") { return(x[7])} else { return(x[8]) }; if (x[8]==1){return(0.999999999999)} else {return(x[8])} })

if (length(which(as.numeric(data$new_sidak) == 1)) > 0) { data[which(as.numeric(data$new_sidak) == 1),]$new_sidak<-0.999999999999 }
if (length(which(as.numeric(data[,9]) == 1)) > 0) { data[which(as.numeric(data[,9]) == 1),9]<-0.999999999999 }

data[which(data[,9]==1),9]<-0.999999999999
data_fisher<-data[which(data[,9] != "NA"),9]
data_fisher_fdr<-locfdr(qnorm(data_fisher,lower.tail=F),plot=0)
data[which(data[,9] != "NA"),"FORGE_locfdr"]<-data_fisher_fdr$fdr

data_sidak<-data[which(data[,8] != "NA"),"new_sidak"]
data_sidak_fdr<-locfdr(qnorm(as.numeric(data_sidak),lower.tail=F),plot=0)
data$sidak_locfdr<-data_sidak_fdr$fdr

data[which(data[,12]==1),12]<-0.999999999999
data_galwey<-data[which(data[,12] != "NA"),12]
data_galwey_fdr<-locfdr(qnorm(as.numeric(data_galwey),lower.tail=F),plot=0)
data[which(data[,12] != "NA"),"galwey_locfdr"]<-data_galwey_fdr$fdr

write.table(data,file=fileout, quote=F, row.names=F,col.names=T,sep = "\t")
