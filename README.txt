This is the README for the script to calculate the gene wide p-value using the Fisher method.
Yout should first read the "Doc" file present on this folder. 

PREREQUISITES
You will need to have Perl install in you system. Additionally you will need the following Perl libraries
PDL
PDL::GSL::CDF
PDL::Stats::Basic
PDL::GSL::RNG
Statistics::RankCorrelation

In addition the software make use of the GSL library for some calculation, see attached documentation (doc.txt) for details on that.

FILES ON THE FOLDER
In the example folder you will find a serie of files which you can use to test.

example.assoc: example file for the -a option. It has some association results for the IL23R gene.
example.bed, example.bim and example.fam: genotype data in PLINK binary format. It has some genotypes for the IL23R gene.
example.out: output file example
example.out.correlation: out file with correlation values.
example.map and example.ped: genotypes in PAD and MAP format 
example.snpmap: example of snp-to-gene annotation



In addition I provide here utility (hopefully useful) files
- The Affy id to rsid mapping (affy_to_rsID.tab), based on the SNP Arrays tables from UCSC web browser. The tables I used had information from dbSNP130.

RUNNING THE SCRIPT.
To run the script with the example files and utility files run the command below on the terminal. In this case I only have data for the IL23R gene so I will restrict the analysis to it.
perl forgev0.7.pl -bfile example/example --a example/example.assoc --m example/example.snpmap --out test -genes IL23R


The standard output should look like.
Thu Jun 24 11:40:44 2010	Max SNP-to-gene distance allowed [ 20 ] kb
Thu Jun 24 11:40:44 2010	Read Gene List command line [ IL23R  ]
Thu Jun 24 11:40:44 2010	Reading association file: [ example/example.assoc ]
Thu Jun 24 11:40:44 2010	[ 19 ] SNPs with association data
Thu Jun 24 11:40:44 2010	Reading SNPs info from [ example/example.bim ]
Thu Jun 24 11:40:44 2010	[ 20 ] SNPs on BED file
Thu Jun 24 11:40:44 2010	Reading samples info from [ example/example.fam ]
Thu Jun 24 11:40:44 2010	[ 4806 ] samples read
Thu Jun 24 11:40:44 2010	Reading gene to SNP mapping file from [ example/example.snpmap ]
Thu Jun 24 11:40:44 2010	  '->[ 1 ] Genes will be analyzed 
Thu Jun 24 11:40:44 2010	  '->[ 17 ] SNPs mapped to Genes and with association results will be analyzed 
Thu Jun 24 11:40:44 2010	Output file will be written to [ test ]
Thu Jun 24 11:40:44 2010	Reading genotypes from [ example/example.bed ]
Thu Jun 24 11:40:44 2010	Binary file is on SNP-major format
Thu Jun 24 11:40:52 2010	Starting to Calculate gene p-values
close() on unopened filehandle GEN0 at /usr/lib64/perl5/5.8.5/x86_64-linux-thread-multi/IO/Handle.pm line 366.
Thu Jun 24 11:40:52 2010	Well Done!!


the output file should look like this line
Ensembl_ID	Hugo_id	gene_type	chromosome	start(ensemblv51+20kb)	end(ensemblv51+20kb)	min_p	min_p_sidak_corrected_by_number_of_tests	fisher_combined_gene-pvalue	chi-square	degrees_of_freedom	Galwey_pvalue	Galwey_chi_square	Galwey_degrees_freedom	number_of_snps	number_effective_tests(Gaoetal;PDMID:19434714)
ENSG00000162594	IL23R	protein_coding	1	67632083	67725662	3.491e-02	4.132e-01	2.128e-01	8.42305	6.04941	1.119e-01	28.17326	10.116917	 15



If you used the -print_cor flag the beggining of the correlation file like 
rs6660226 rs11209018 0.090
rs6660226 rs10489628 0.311
rs6660226 rs6664119 0.044
rs6660226 rs2201841 0.447
rs6660226 rs41339545 0.430
rs6660226 rs10489629 0.405
rs6660226 rs10889671 0.288
rs6660226 rs790633 0.097
rs6660226 rs790632 0.072
rs6660226 rs11465803 0.278


You may want to use genotype data from the HapMap (you can get it from here http://ftp.hapmap.org/phase_3/?N=D) for the -bfile option. In general this is quicker than using most current GWAS sample's genotypes but the price comes on possible diff in allele frequency between your sample and the HapMap. The software can handle large studies, it will just take longer. On my personal experience about a day for a 5000 samples study with ~350000 SNPs. To speed up things you can use the --chr flag ( remember to use the chromosome name as it appears in the SNP-to-gene mapping file) and run the analysis on each chromosome separately or used the --gene_list option to give a list of genes identifiers.

Known problems
- In some Perl installation the GetOpt::Long module behaves differently and it does not accept multiple argument per option. You will not it because of the error message

Error in option spec: "genes=s{1,1000}"

If this happens you need to change the options specification

'genes=s{1,1000}' => \@genes,

have to look like

   'genes=s' => \@genes,


This is the end of the README.
Any questions or feedback please e-mail Inti Pedroso <intipedroso@gmail.com> or Gerome Breen <gerome.breen@kcl.ac.uk>.


Inti Pedroso
24th June 2010
