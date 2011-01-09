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
	data<-read.table(filein,sep = "\t",header=T) 
} else {
	data<-read.table(filein,sep = "\t") 
}
# recalculate the FORGE p-values. This because usually does not work on double precission
data[which(data[,24] != 1),9]<-pchisq(data[which(data[,24] != 1),10],df= data[which(data[,24] != 1),11],lower.tail=F); 
if (length(which(as.numeric(data[,24]) == 1)) > 0) { 
	data[which(data[,24] == 1),9]<- data[which(data[,24] == 1),7] # add fisher p-values  
	data[which(data[,24] == 1),8]<- data[which(data[,24] == 1),7] # add sidak p-values
	
	data[which(data[,24] == 1),14]<- data[which(data[,24] == 1),7] # add Z_P_fix 
	
	data[which(data[,24] == 1),17]<- data[which(data[,24] == 1),7] # add Z_P_random 
}

# re-order that data
data<-data[order(data[,14]),]

#SIDAK
data$new_sidak<-data[,8]
data[which(data$new_sidak==1),"new_sidak"]<-0.999999999999
data[which(data$new_sidak==0),"new_sidak"]<-data[which(data$new_sidak==0),7]*data[which(data$new_sidak==0),24]

data_sidak<-data[which(is.na(data[,"new_sidak"]) == FALSE),"new_sidak"]
data_sidak_fdr<-locfdr(qnorm(as.numeric(data_sidak),lower.tail=F),plot=0)
data[which(is.na(data[,"new_sidak"]) == FALSE),"SIDAK_locfdr"]<-data_sidak_fdr$fdr

#FISHER
data$new_fisher<-data[,9]
data[which(data$new_fisher==1),"new_fisher"]<-0.999999999999
data[which(data$new_fisher==0),"new_fisher"]<-data[which(data$new_fisher==0),7]*data[which(data$new_fisher==0),24]

data_fisher<-data[which(is.na(data[,"new_fisher"]) == FALSE),"new_fisher"]
data_fisher_fdr<-locfdr(qnorm(as.numeric(data_fisher),lower.tail=F),plot=0)
data[which(is.na(data[,"new_fisher"]) == FALSE),"FISHER_locfdr"] <-data_fisher_fdr$fdr

#Z-fix
data$new_Z_FIX_P<-data[,14]
data[which(data$new_Z_FIX_P==1),"new_Z_FIX_P"]<-0.999999999999
data[which(data$new_Z_FIX_P==0),"new_Z_FIX_P"]<-data[which(data$new_Z_FIX_P==0),7]*data[which(data$new_Z_FIX_P==0),24]

data_Z_FIX_P<-data[which(is.na(data[,"new_Z_FIX_P"]) == FALSE),"new_Z_FIX_P"]
data_Z_FIX_P_fdr<-locfdr(qnorm(as.numeric(data_Z_FIX_P),lower.tail=F),plot=0)
data[which(is.na(data[,"new_Z_FIX_P"]) == FALSE),"Z_FIX_P_locfdr"] <-data_Z_FIX_P_fdr$fdr

#Z-RANDOM
data$new_Z_RANDOM_P<-data[,17]
data[which(data$new_Z_RANDOM_P==1),"new_Z_RANDOM_P"]<-0.999999999999
data[which(data$new_Z_RANDOM_P==0),"new_Z_RANDOM_P"]<-data[which(data$new_Z_RANDOM_P==0),7]*data[which(data$new_Z_RANDOM_P==0),24]

data_Z_RANDOM_P<-data[which(is.na(data[,"new_Z_RANDOM_P"]) == FALSE),"new_Z_RANDOM_P"]
data_Z_RANDOM_P_fdr<-locfdr(qnorm(as.numeric(data_Z_RANDOM_P),lower.tail=F),plot=0)
data[which(is.na(data[,"new_Z_RANDOM_P"]) == FALSE),"Z_RANDOM_P_locfdr"]<-data_Z_RANDOM_P_fdr$fdr

write.table(data,file=fileout, quote=F, row.names=F,col.names=T,sep = "\t")
