#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false", dest="verbose", help="Print little output"),
    make_option(c("--bfile"), type="character", default=NULL, help="Files in PLINK binary format, corresping file with extension *.bed, *.bim and *.fam.", metavar=NULL),
    make_option(c("--assoc_file"), type="character", default=NULL, help="SNP association results. If header is provided please use also --header", metavar=NULL),
    make_option(c("-gc", "--gc_correction"), action="store_true", default=FALSE, help="Apply genomic control correction"),
    make_option(c("--lambda"), type="integer", default=1, help="Lambda value for genomic control correction"),
    make_option(c("-o","--outfile"), type="character", default="forge.out.txt", help="Output file Name", metavar=NULL),
    make_option(c("--report"), type="integer", default=100, help="How often to report advance", metavar=NULL),
    make_option(c("--flush"), type="integer", default=100, help="How often print out results", metavar=NULL),
    make_option(c("--max_step"), type="integer", default=10000, help="In MND simulation, max number of simulation to perform on each step", metavar=NULL),
    make_option(c("--max_mnd"), type="integer", default=10000000, help="Max number of MND simulations possible", metavar=NULL),
    make_option(c("--target"), type="integer", default=10, help="How many times the observed statistics must be seen", metavar=NULL),
    make_option(c("-d","--distance"), type="integer", default=20, help="Max SNP-to-gene distance allowed (in kb)", metavar=NULL),
    make_option(c("--save_annot"), type="character", default='FALSE', help="Save gene and SNP annotation data used", metavar=NULL),
    make_option(c("--affy2rsid"), type="character", default='FALSE', help="Affy id to rsid mapping file", metavar=NULL),
    make_option(c("--header"), action="store_true", default=FALSE, help="SNP association file has header", metavar=NULL),
    make_option(c("--snp_annot"), type="character", default=NULL, help="File with SNP annotations", metavar=NULL),
    make_option(c("--gene_annot"), type="character", default=NULL, help="File with gene annotation", metavar=NULL),
    make_option(c("--annot_folder"), type="character", default="annotation/ensemblv62/", help="Folder with RData SNP annotations. default: annotation/ensemblv62/ ", metavar=NULL)
)

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
######
cat("Loading necessary libraries\n")
cat("snpMatrix\n")
suppressPackageStartupMessages(library(snpMatrix))
cat("corpcor\n")
suppressPackageStartupMessages(library(corpcor))
if ( is.null(opt$annot_folder) == TRUE ){
    cat("biomaRt\n")
    suppressPackageStartupMessages(library(biomaRt))
}
cat("IRanges\n")
suppressPackageStartupMessages(library(IRanges))
cat("local functions [ functions.R ]\n")
source("functions.R")
cat("Done with libraries\n")

cat("Reading Genotypes from [ ",opt$bfile," ]\n",sep="")

bed<-paste(opt$bfile,".bed",sep="")
bim<-paste(opt$bfile,".bim",sep="")
fam<-paste(opt$bfile,".fam",sep="")

genotypes<-read.plink(bed,bim,fam)


cat("Reading SNP p-values from [ ",opt$assoc_file," ]\n",sep="")
assoc<-read.table(opt$assoc_file,header=opt$header)
cat("   '-> [ ",nrow(assoc)," ] rows read\n",sep="")

if(ncol(assoc) < 2){
  cat("Only found one columns in file")
  exit()
}
if (ncol(assoc) == 2){
  colnames(assoc)<-c("SNP","P")
}

if (opt$gc_correction == TRUE) {
  cat("Evaluating genomic inflation factor\n")
  l<-lambda(assoc$P)
  cat("   - lambda = [ ",l," ]\n")
  if (l>1){
    assoc$P<-apply_gc(p=assoc$P,lambda = l)      
    l<-lambda(assoc$P)
    cat("   - After GC lambda = [ ",l," ]\n")
  }
}


if (opt$affy2rsid != "FALSE"){
    cat("Matching Affy ids with file [ ",opt$affy2rsid," ]\n")
    affy2rsid<-read.table(affy2rsid)
    colnames(affy2rsid)<-c("affy","SNP")
    cat("   '-> Matching ids\n")
    with_affy_ids<-merge(affy2rsid,assoc,by.x="SNP",by.y="SNP")
    cat("   '-> Adding ids to data\n")
    assoc<-rbind(assoc,with_affy_ids[,c("SNP","P")])
}

ENSEMBL_GENES<-NULL
if (is.null(opt$gene_annot) ==FALSE){
    cat("Reading gene annotation from [ ", opt$gene_annot," ]\n",sep="")
    ENSEMBL_GENES<-read.table(opt$gene_annot,header=T,sep="\t",colClasses=c('character','character','character','numeric','numeric','character','character'))
} else if (is.null(opt$annot_folder) == TRUE) {
    cat("Getting gene annotation from ensembl\n")
    genemart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    ENSEMBL_GENES<-getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position","end_position", "band","gene_biotype"), mart = genemart)
} 

snp<-NULL
if (is.null(opt$snp_annot) ==FALSE){
    cat("Reading SNP annotation from [ ", opt$snp_annot," ]\n",sep="")
    snp<-read.table(opt$snp_annot,header=T,sep="\t",colClasses=c('character','character','numeric','numeric','numeric','character'))
} else if (is.null(opt$annot_folder) == TRUE ){
    all_chrs<-sort(unique(ENSEMBL_GENES$chromosome_name))
    # get set of SNPs with stats and genotypes
    genotyped_snps<-as.data.frame(colnames(genotypes))
    colnames(genotyped_snps)<-c("SNP")
    working_snps<-merge(genotyped_snps,assoc,by.x="SNP",by.y="SNP")
    
    cat("[ ",nrow(working_snps)," ] SNPs with genotypes and statitics\n",sep="")
    snpmart = useMart("snp", dataset = "hsapiens_snp")
    snp<-getBM(c("refsnp_id","chr_name","chrom_start", "chrom_strand","mapweight","allele"), filters="refsnp",values=working_snps$SNP,mart = snpmart)
} 


if (opt$save_annot!="FALSE"){
    cat("Saving Gene and SNP annotation to [ ",opt$save_annot,".gene_annot.txt ] and [ ",opt$save_annot,".snp_annot.txt ]\n",sep="")
    snp_out<-paste(opt$save_annot,".snp_annot.txt",sep="")
    gene_out<-paste(opt$save_annot,".snp_annot.txt",sep="")
    write.table(snp,file=snp_out,sep="\t",col.names=T,row.names=F,quote=F)
    write.table(ENSEMBL_GENES,file=gene_out,sep="\t",col.names=T,row.names=F,quote=F)
}

DATA_OUT<-NULL
header<-c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","band","gene_biotype","SIDAK","FISHER","Z_FIX","Z_RANDOM","SEEN_SIDAK","SEEN_FISHER","SEEN_Z_FIX","SEEN_Z_RANDOM","N","I2","Q","tau_squared","N_SNPs")
write.table(t(header), file=opt$outfile,append=F,col.names=F,row.names=F,quote=F,sep="\t")

assoc_genotyped<-NULL
if ( is.null(opt$annot_folder) == FALSE ) {
    cat("Reading SNP annotation from folder [ ",opt$annot_folder," ]\n")
    # get set of SNPs with stats and genotypes
    genotyped_snps<-as.data.frame(colnames(genotypes))
    colnames(genotyped_snps)<-c("SNP")
    assoc_genotyped<-merge(assoc,genotyped_snps,assoc,by.x="SNP",by.y="SNP")
}
gene_counter=0
# get list of all chromosomes
all_chrs<-sort(unique(ENSEMBL_GENES$chromosome_name))
# loop over all chromosomes 
all_chrs=22
for(chr in all_chrs){ 
    if ( is.null(opt$annot_folder) == FALSE ) {
        # loop over all chromosomes
        #rm(snp)
        cat("Loading SNP data for chromsome [ ",chr," ]\n",sep="")
        file<-paste(opt$annot_folder,"chromosome.",chr,".annot.RData",sep="")
        load(file)
        # reduce SNP set to the working set
        working_snps<-merge(snp,assoc_genotyped,by.x="refsnp_id",by.y="SNP")
        snp<-unique(working_snps)
        cat("    - [ ",nrow(snp)," ] SNPs with genotypes and statitics\n",sep="")
    }
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
    gene_range<-IRanges(start=gene_chr$start_position - opt$distance*1000 ,end=gene_chr$end_position + opt$distance*1000)
    map<-findOverlaps(query=gene_range,subject=snp_range,type="any",select="all")
    map<-as.matrix(map)
    cat("[ ",length(unique(map[,1]))," ] genes mapped to SNPs with genotypes and association results\n");
    cat("Starting to Analyse Genes\n")
    # loop over each gene and perform the analyses
    for ( my_gene_index in unique(map[,1])){
        report_advance(gene_counter,opt$report,"Genes")
        gene_counter<-gene_counter+1
        gene_data<-gene_chr[my_gene_index,]
        # get snps in the gene
        gene_snps<-snp_chr[map[which(map[,1]==my_gene_index),2],]
        if ( nrow(gene_snps) == 0){ 
            #cat("No SNPs mapped to gene [ ",my_gene," ]\n",sep="")
            next
        } else {
            #cat(nrow(gene_snps)," SNPs mapped to gene [ ",gene_data$ensembl_gene_id," ",gene_data$hgnc_symbol," ]\n",sep="")  
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
  
            gene_snps[which(gene_snps$P > (1 - .Machine$double.eps)),"P"]<-1 - .Machine$double.eps
            w<-rep(1/length(gene_snps$P),length(gene_snps$P))
            min_p<-min(gene_snps$P)
            sidak<-1-(1 - min_p)^nrow(gene_genotypes_cor)
            fisher<-modified_fisher(p=gene_snps$P,w=w,cor=gene_genotypes_cor)
            w<-rep(1/nrow(gene_genotypes_cor),nrow(gene_genotypes_cor))
            z<-qnorm(gene_snps$P,lower.tail=F)

            Z_methods<-z_fix_and_random_effects(z=z,w=w,cov=gene_genotypes_cor)
            SEEN<-rep(0,4)
            TOTAL<-0    
            running_step<-10
            while (min(SEEN) < opt$target){
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
                if (running_step > opt$max_step){
                    running_step<-opt$max_step  
                }
                #cat(TOTAL,SEEN,"\n")
                if (TOTAL > opt$max_mnd){ break;}
            }# while 
            MND_P<-(SEEN+1)/(TOTAL+1)
            ready<-unlist(c(gene_data,MND_P,SEEN,TOTAL,Z_methods$I2,Z_methods$Q,Z_methods$tau_squared,nrow(gene_genotypes_cor)))
        } # else more than 1 SNP 
        if (is.null(DATA_OUT) == TRUE) { 
            DATA_OUT<-as.data.frame(t(ready))  
        } else {
            DATA_OUT<-rbind(DATA_OUT,t(ready))
        }
        if(nrow(DATA_OUT) > opt$flush){
            write.table(DATA_OUT, file=opt$outfile,append=T,col.names=F,row.names=F,quote=F,sep="\t")
            DATA_OUT<-NULL
        }
    }
    write.table(DATA_OUT, file=opt$outfile,append=T,col.names=F,row.names=F,quote=F,sep="\t")
} # for loop chromosomes