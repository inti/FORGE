outfile<-"forge.out.txt"
report=100
flush=10
max_step<-10000
step<-10
MAX<-1000000
distance<-20
save_annot<-"FALSE"
affy2rsid<-"FALSE"
header<-"FALSE"
snp_annot<-NULL
gene_annot<-NULL

bfile="hm3_ceu.chr22"
assoc_file="~/DATA/WTCCC/WEB/basic/WTCCC_BD_imputed.txt_plusQC.filtered.pvals"
report=10
#snp_annot="29April2011.snp_annot.txt"
#gene_annot="29April2011.gene_annot.txt"


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
report<-as.numeric(report)
flush<-as.numeric(flush)
max_step<-as.numeric(max_step)
step<-as.numeric(step)
MAX<-as.numeric(MAX)
distance<-as.numeric(distance)*1000
header<-as.logical(header)
######
cat("Loading necessary libraries\n")
cat("snpMatrix\n")
library(snpMatrix)
cat("corpcor\n")
library(corpcor)
cat("biomaRt\n")
library(biomaRt)
cat("IRanges\n")
library(IRanges)
cat("local functions [ functions.R ]\n")
source("functions.R")
cat("Done with libraries\n")

cat("Reading Genotypes from [ ",bfile," ]\n",sep="")

bed<-paste(bfile,".bed",sep="")
bim<-paste(bfile,".bim",sep="")
fam<-paste(bfile,".fam",sep="")

genotypes<-read.plink(bed,bim,fam)

if (is.null(gene_annot) ==FALSE){
    cat("Reading gene annotation from [ ", gene_annot," ]\n",sep="")
    ENSEMBL_GENES<-read.table(gene_annot,header=T,sep="\t")
} else {
    cat("Getting gene annotation from ensembl\n")
    genemart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    ENSEMBL_GENES<-getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position","end_position", "band","gene_biotype"), mart = genemart)
}
cat("Reading SNP p-values from [ ",assoc_file," ]\n",sep="")
assoc<-read.table(assoc_file,header=header)
cat("   '-> [ ",nrow(assoc)," ] rows read\n",sep="")

if(ncol(assoc) < 2){
  cat("Only found one columns in file")
  exit()
}
if (ncol(assoc) == 2){
  colnames(assoc)<-c("SNP","P")
}

if (affy2rsid != "FALSE"){
    cat("Matching Affy ids with file [ ",affy2rsid," ]\n")
    affy2rsid<-read.table(affy2rsid)
    colnames(affy2rsid)<-c("affy","SNP")
    cat("   '-> Matching ids\n")
    with_affy_ids<-merge(affy2rsid,assoc,by.x="SNP",by.y="SNP")
    cat("   '-> Adding ids to data\n")
    assoc<-rbind(assoc,with_affy_ids[,c("SNP","P")])
}

# get set of SNPs with stats and genotypes
genotyped_snps<-as.data.frame(colnames(genotypes))
colnames(genotyped_snps)<-c("SNP")
working_snps<-merge(genotyped_snps,assoc,by.x="SNP",by.y="SNP")
cat("[ ",nrow(working_snps)," ] SNPs with genotypes and statitics\n",sep="")

cat("Deleting genotypes of unused SNPs\n",sep="")
#genotypes<-genotypes[,working_snps$SNP]
snp<-NULL
if (is.null(snp_annot) ==FALSE){
    cat("Reading SNP annotation from [ ", snp_annot," ]\n",sep="")
    snp<-read.table(snp_annot,header=T,sep="\t")
} else {
    cat("Annotating SNPs\n")
    snpmart = useMart("snp", dataset = "hsapiens_snp")
    snp<-getBM(c("refsnp_id","chr_name","chrom_start", "chrom_strand","mapweight","allele"), filters="refsnp",values=working_snps$SNP,mart = snpmart)
}

if (save_annot!="FALSE"){
    cat("Saving Gene and SNP annotation to [ ",save_annot,".gene_annot.txt ] and [ ",save_annot,".snp_annot.txt ]\n",sep="")
    snp_out<-paste(save_annot,".snp_annot.txt",sep="")
    gene_out<-paste(save_annot,".gene_annot.txt",sep="")
    write.table(snp,file=snp_out,sep="\t",col.names=T,row.names=F,quote=F)
    write.table(ENSEMBL_GENES,file=gene_out,sep="\t",col.names=T,row.names=F,quote=F)
}
# reduce SNP set to the working set
snp<-merge(snp,working_snps,by.x="refsnp_id",by.y="SNP")

DATA_OUT<-NULL
header<-c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","band","gene_biotype","SIDAK","FISHER","Z_FIX","Z_RANDOM","SEEN_SIDAK","SEEN_FISHER","SEEN_Z_FIX","SEEN_Z_RANDOM","N","I2","Q","tau_squared","N_SNPs")
write.table(t(header), file=outfile,append=F,col.names=F,row.names=F,quote=F,sep="\t")
# get list of all chromosomes
all_chrs<-unique(ENSEMBL_GENES$chromosome_name)
gene_counter=0
# loop over all chromosomes
for(chr in all_chrs){ 
    # select SNPs in the chromosome
    snp_chr<-snp[which(snp$chr_name == chr),]
    # check if there are SNPs mapped to the chromsome
    if ( nrow(snp_chr) == 0){ 
        next
    }
    cat("Woking on chromsome [ ",chr," ]\n",sep="")
    # make SNP range
    snp_range<-IRanges(start=snp_chr$chrom_star,end=snp_chr$chrom_star+1)

    # select genes in the chromosome
    gene_chr<-ENSEMBL_GENES[which(ENSEMBL_GENES$chromosome_name==chr),]
    # Map SNPs to genes
    cat("Map SNPs to Genes\n")
    gene_range<-IRanges(start=gene_chr$start_position - distance ,end=gene_chr$end_position + distance)
    map<-findOverlaps(query=gene_range,subject=snp_range,type="any",select="all")
    map<-as.matrix(map)
    cat("[ ",length(unique(map[,1]))," ] genes mapped to SNPs with genotypes and association results\n");
    cat("Starting to Analyse Genes\n")
    # loop over each gene and perform the analyses
    for ( my_gene_index in unique(map[,1])){
        report_advance(gene_counter,report,"Genes")
        gene_counter<-gene_counter+1
        gene_data<-gene_chr[my_gene_index,]
        # get snps in the gene
        gene_snps<-snp_chr[map[which(map[,1]==my_gene_index),2],]
        if ( nrow(gene_snps) == 0){ 
            #cat("No SNPs mapped to gene [ ",my_gene," ]\n",sep="")
            next
        } else {
            cat(nrow(gene_snps)," SNPs mapped to gene [ ",gene_data$ensembl_gene_id," ",gene_data$hgnc_symbol," ]\n",sep="")  
        }
        if (nrow(gene_snps) == 1){
            MND_P<-rep(gene_snps$P,4)
            SEEN<-rep(-1,4)
            TOTAL<-0
            ready<-unlist(c(gene_data,MND_P,SEEN,TOTAL,-1,-1,-1,1))
        } else {
            # extract genotypes
            gene_genotypes<-as(genotypes[,gene_snps$refsnp_id],'numeric')
            gene_genotypes[is.na(gene_genotypes)]<-4
            # calculate pearson correlation
            #gene_genotypes_cor<-cor(gene_genotypes,use="pairwise.complete.obs")
            gene_genotypes_cor<-cor.shrink(gene_genotypes,verbose=FALSE)
            gene_genotypes_cor<-unclass(gene_genotypes_cor)
            # set NA due to missing data to 0
            gene_genotypes_cor[is.na(gene_genotypes_cor)]<-0
            # make the correlation matrix positive definite
            gene_genotypes_cor<-my_make_positive_definite(gene_genotypes_cor)
            # calculate gene p-values
            # apply simulatin based method

            w<-rep(1/length(gene_snps$P),length(gene_snps$P))
            min_p<-min(gene_snps$P)
            sidak<-1-(1 - min_p)^nrow(gene_genotypes_cor)
            fisher<-modified_fisher(p=gene_snps$P,w=w,cor=gene_genotypes_cor)
            w<-rep(1/nrow(gene_genotypes_cor),nrow(gene_genotypes_cor))
            z<-qnorm(gene_snps$P,lower.tail=F)

            Z_methods<-z_fix_and_random_effects(z=z,w=w,cov=gene_genotypes_cor)
            SEEN<-rep(0,4)
            TOTAL<-0    
            running_step<-step
            while (min(SEEN) < 10){
                my_sim_z<-rmvnorm(running_step,sigma=gene_genotypes_cor)
                index=1
                for (index in 1:nrow(my_sim_z)){ 
                    z_sim<-my_sim_z[index,]
                    sim_p<-pchisq(z_sim^2,df=1,lower.tail=F)
                    z_sim<-qnorm(sim_p,lower.tail=F)
                    sim_Z_methods<-z_fix_and_random_effects(z=z_sim,w=w,cov=gene_genotypes_cor); 
                    sim_min_p<-min(sim_p)
                    sim_fisher<-modified_fisher(p=sim_p,w=w,cor=gene_genotypes_cor)
                    if (sim_min_p <= min_p){
                        SEEN[1]<-SEEN[1]+1
                    }
                    if (sim_fisher$p <= fisher$p){
                        SEEN[2]<-SEEN[2]+1          
                    }
                    if (sim_Z_methods$P_FIX <= Z_methods$P_FIX){
                        SEEN[3]<-SEEN[3]+1
                    }
                    if (sim_Z_methods$P_RANDOM <= Z_methods$P_RANDOM){
                        SEEN[4]<-SEEN[4]+1
                    }
                } # for
                TOTAL<-TOTAL+running_step
                running_step<-running_step*10
                if (running_step > max_step){
                    running_step<-max_step  
                }
                #cat(TOTAL,SEEN,"\n")
                if (TOTAL > MAX){ break;}
            }# while 
            MND_P<-(SEEN+1)/(TOTAL+1)
            ready<-unlist(c(gene_data,MND_P,SEEN,TOTAL,Z_methods$I2,Z_methods$Q,Z_methods$tau_squared,nrow(gene_genotypes_cor)))
        } # else more than 1 SNP 
        if (is.null(DATA_OUT) == TRUE) { 
            DATA_OUT<-as.data.frame(t(ready))  
        } else {
            DATA_OUT<-rbind(DATA_OUT,t(ready))
        }
        if(nrow(DATA_OUT) > flush){
            write.table(DATA_OUT, file=outfile,append=T,col.names=F,row.names=F,quote=F,sep="\t")
            DATA_OUT<-NULL
        }
    }
    write.table(DATA_OUT, file=outfile,append=T,col.names=F,row.names=F,quote=F,sep="\t")
} # for loop chromosomes

