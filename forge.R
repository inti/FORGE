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
    make_option(c("--annot_folder"), type="character", default="annotation/ensemblv62/", help="Folder with RData SNP annotations. default: annotation/ensemblv62/ ", metavar=NULL),
    make_option(c("--chr"), type="character", default=NULL, help="Analyze this chromosome", metavar=NULL),
    make_option(c("--asymp"), action="store_true", default=TRUE,  help="Perform analysis using asymptotic approximations", metavar=NULL),
    make_option(c("--mnd"), action="store_true", default=FALSE,  help="Perform analysis using MND simulation approximation", metavar=NULL)
)

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
######
source("functions.R")
print_OUT("local functions [ functions.R ]")
print_OUT("Loading necessary libraries")
print_OUT("snpMatrix")
suppressPackageStartupMessages(library(snpMatrix))
print_OUT("corpcor")
suppressPackageStartupMessages(library(corpcor))
if ( is.null(opt$annot_folder) == TRUE ){
    print_OUT("biomaRt")
    suppressPackageStartupMessages(library(biomaRt))
}
print_OUT("IRanges")
suppressPackageStartupMessages(library(IRanges))
print_OUT("Done with libraries")

print_OUT(paste("Reading Genotypes from [ ",opt$bfile," ]",sep=""))

bed<-paste(opt$bfile,".bed",sep="")
bim<-paste(opt$bfile,".bim",sep="")
fam<-paste(opt$bfile,".fam",sep="")

genotypes<-read.plink(bed,bim,fam)


print_OUT(paste("Reading SNP p-values from [ ",opt$assoc_file," ]",sep=""))
assoc<-read.table(opt$assoc_file,header=opt$header)
print_OUT(paste("   '-> [ ",nrow(assoc)," ] rows read",sep=""))

if(ncol(assoc) < 2){
  print_OUT("Only found one columns in file")
  exit()
}
if (ncol(assoc) == 2){
  colnames(assoc)<-c("SNP","P")
}
assoc<-assoc[which(is.na(assoc$P) == FALSE),]
if (opt$gc_correction == TRUE) {
  print_OUT("Evaluating genomic inflation factor")
  assoc[which(assoc$P > (1 - .Machine$double.eps)),"P"]<-1 - .Machine$double.eps
  l<-lambda(assoc$P)
  print_OUT(paste("   - lambda = [ ",l," ]"))
  if (l>1){
    assoc$P<-apply_gc(p=assoc$P,lambda = l)      
    l<-lambda(assoc$P)
    print_OUT(paste("   - After GC lambda = [ ",l," ]"))
  }
}


if (opt$affy2rsid != "FALSE"){
    print_OUT(paste("Matching Affy ids with file [ ",opt$affy2rsid," ]"))
    affy_map_rsid<-read.table(opt$affy2rsid)
    colnames(affy_map_rsid)<-c("affy","SNP")
    print_OUT("   '-> Matching ids")
    with_affy_ids<-merge(affy_map_rsid,assoc,by.x="SNP",by.y="SNP")
    print_OUT("   '-> Adding ids to data")
    assoc<-rbind(assoc,with_affy_ids[,c("SNP","P")])
}

ENSEMBL_GENES<-NULL
if (is.null(opt$gene_annot) ==FALSE){
    print_OUT(paste("Reading gene annotation from [ ", opt$gene_annot," ]",sep=""))
    ENSEMBL_GENES<-read.table(opt$gene_annot,header=T,sep="\t",colClasses=c('character','character','character','numeric','numeric','character','character'))
} else if (is.null(opt$annot_folder) == TRUE) {
    print_OUT("Getting gene annotation from ensembl")
    genemart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    ENSEMBL_GENES<-getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position","end_position", "band","gene_biotype"), mart = genemart)
} 

snp<-NULL
if (is.null(opt$snp_annot) ==FALSE){
    print_OUT(paste("Reading SNP annotation from [ ", opt$snp_annot," ]",sep=""))
    snp<-read.table(opt$snp_annot,header=T,sep="\t",colClasses=c('character','character','numeric','numeric','numeric','character'))
} else if (is.null(opt$annot_folder) == TRUE ){
    all_chrs<-sort(unique(ENSEMBL_GENES$chromosome_name))
    # get set of SNPs with stats and genotypes
    genotyped_snps<-as.data.frame(colnames(genotypes))
    colnames(genotyped_snps)<-c("SNP")
    working_snps<-merge(genotyped_snps,assoc,by.x="SNP",by.y="SNP")
    
    print_OUT(paste("[ ",nrow(working_snps)," ] SNPs with genotypes and statitics",sep=""))
    snpmart = useMart("snp", dataset = "hsapiens_snp")
    snp<-getBM(c("refsnp_id","chr_name","chrom_start", "chrom_strand","mapweight","allele"), filters="refsnp",values=working_snps$SNP,mart = snpmart)
} 


if (opt$save_annot!="FALSE"){
    print_OUT(paste("Saving Gene and SNP annotation to [ ",opt$save_annot,".gene_annot.txt ] and [ ",opt$save_annot,".snp_annot.txt ]",sep=""))
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
    print_OUT(paste("Reading SNP annotation from folder [ ",opt$annot_folder," ]"))
    # get set of SNPs with stats and genotypes
    genotyped_snps<-as.data.frame(colnames(genotypes))
    colnames(genotyped_snps)<-c("SNP")
    assoc_genotyped<-merge(assoc,genotyped_snps,assoc,by.x="SNP",by.y="SNP")
}
gene_counter=0
# get list of all chromosomes
all_chrs<-sort(unique(ENSEMBL_GENES$chromosome_name))
# loop over all chromosomes 
if (is.null(opt$chr) == FALSE){
    all_chrs<-opt$chr  
    print_OUT(paste("Will analyze chromosome [",all_chrs," ]"))
}
for(chr in all_chrs){ 
    if ( is.null(opt$annot_folder) == FALSE ) {
        print_OUT(paste("Loading SNP data for chromsome [ ",chr," ]",sep=""))
        file<-paste(opt$annot_folder,"chromosome.",chr,".annot.RData",sep="")
        load(file)
        # reduce SNP set to the working set
        working_snps<-merge(snp,assoc_genotyped,by.x="refsnp_id",by.y="SNP")
        snp<-unique(working_snps)
        print_OUT(paste("    - [ ",nrow(snp)," ] SNPs with genotypes and statitics",sep=""))
    }
    # select SNPs in the chromosome
    snp_chr<-snp[which(snp$chr_name == chr),]
    # check if there are SNPs mapped to the chromsome
    if ( nrow(snp_chr) == 0){ 
        next
    }
    print_OUT(paste("Woking on chromsome [ ",chr," ]",sep=""))
    # make SNP range
    snp_range<-IRanges(start=snp_chr$chrom_star,end=snp_chr$chrom_star+1)

    # select genes in the chromosome
    gene_chr<-ENSEMBL_GENES[which(ENSEMBL_GENES$chromosome_name==chr),]
    # Map SNPs to genes
    print_OUT("Map SNPs to Genes")
    gene_range<-IRanges(start=gene_chr$start_position - opt$distance*1000 ,end=gene_chr$end_position + opt$distance*1000)
    map<-findOverlaps(query=gene_range,subject=snp_range,type="any",select="all")
    map<-as.matrix(map)
    print_OUT(paste("[ ",length(unique(map[,1]))," ] genes mapped to SNPs with genotypes and association results"));
    print_OUT("Starting to Analyse Genes")
    # loop over each gene and perform the analyses
    for ( my_gene_index in unique(map[,1])){
        report_advance(gene_counter,opt$report,"Genes")
        gene_counter<-gene_counter+1
        gene_data<-gene_chr[my_gene_index,]
        # get snps in the gene
        gene_snps<-snp_chr[map[which(map[,1]==my_gene_index),2],]
        if ( nrow(gene_snps) == 0){ 
            next
        } else {
            #cat(nrow(gene_snps)," SNPs mapped to gene [ ",gene_data$ensembl_gene_id," ",gene_data$hgnc_symbol," ]\n",sep="")  
        }

        if (opt$mnd == TRUE){
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
                min_p<-min(gene_snps$P,na.rm=T)
                sidak<-1-(1 - min_p)^nrow(gene_genotypes_cor)
                fisher<-modified_fisher(p=gene_snps$P,w=w,cor=gene_genotypes_cor)
                w<-rep(1/nrow(gene_genotypes_cor),nrow(gene_genotypes_cor))
                z<-qnorm(gene_snps$P,lower.tail=F)
                gene_genotypes_cor<-gates_transform_corr(gene_genotypes_cor)
                Z_methods<-z_fix_and_random_effects(z=z,w=w,cov=gene_genotypes_cor)
                SEEN<-rep(0,4)
                TOTAL<-0    
                running_step<-10
                while (min(SEEN) < opt$target){
                    my_sim_z<-rmvnorm(running_step,sigma=gene_genotypes_cor)
                    index=1
                    sim_counter<-0
                    for (index in 1:nrow(my_sim_z)){ 
                        z_sim<-my_sim_z[index,]
                        sim_p<-pchisq(z_sim^2,df=1,lower.tail=F)
                        z_sim<-qnorm(sim_p,lower.tail=F)
                        sim_Z_methods<-z_fix_and_random_effects(z=z_sim,w=w,cov=gene_genotypes_cor); 
                        sim_min_p<-min(sim_p,na.rm=T)
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
                        sim_counter<-sim_counter+1
                        if (min(SEEN) > 9){
                          running_step<-sim_counter
                          break;  
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
        } 
        if (opt$asymp == TRUE) {
            if (nrow(gene_snps) == 1){
                MND_P<-rep(gene_snps$P,4)
                ready<-unlist(c(gene_data,rep(gene_snps$P,5),-1,-1,-1,1))
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
                min_p<-min(gene_snps$P,na.rm=T)
                sidak<-1-(1 - min_p)^nrow(gene_genotypes_cor)
                fisher<-modified_fisher(p=gene_snps$P,w=w,cor=gene_genotypes_cor)
                w<-rep(1/nrow(gene_genotypes_cor),nrow(gene_genotypes_cor))
                z<-qnorm(gene_snps$P,lower.tail=F)
                gene_genotypes_cor<-gates_transform_corr(gene_genotypes_cor)
                Z_methods<-z_fix_and_random_effects(z=z,w=w,cov=gene_genotypes_cor)

                ready<-unlist(c(gene_data,min_p,sidak,fisher,Z_methods$P_FIX,Z_methods$P_RANDOM,Z_methods$I2,Z_methods$Q,Z_methods$tau_squared,nrow(gene_genotypes_cor)))                    
            }
        }
        # merge data into output and print out

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

print_OUT("Well Done!!")