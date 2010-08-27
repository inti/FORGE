#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use PDL;
use PDL::GSL::CDF;
use PDL::Primitive;
use PDL::NiceSlice;
use PDL::Stats::Basic; 
# use PDL::GSL::RNG;
use IO::File;
use IO::Seekable;
use Fcntl;
use Data::Dumper;
use Statistics::RankCorrelation;
use Pod::Usage;

my $VERSION = "0.92";


our ( $help, $man, $out, $snpmap, $bfile, $assoc, $gene_list, @genes,
      $all_genes, $analysis_chr, $report, $spearman, $affy_to_rsid,
      $v,  $lambda,$print_cor, $pearson_genotypes,$distance, $sample_score, $ped, $map);

GetOptions(
   'help|h' => \$help,
   'man' => \$man,
   'ped=s' => \$ped,
   'map=s' => \$map,
   'bfile=s'    => \$bfile, 
   'out|o=s'   => \$out, #name of the output file
   'assoc|a=s' => \$assoc,
   'gene_list|g=s'   => \$gene_list,
   'genes=s' => \@genes,
   'all_genes'       => \$all_genes,
   'chr=s'	=> \$analysis_chr,  
   'snpmap|m=s'      => \$snpmap,
   'report=i'  => \$report,
   'correlation|cor=s' => \$spearman,
   'affy_to_rsid=s' => \$affy_to_rsid,
   'verbose|v' => \$v,
   'lambda=f' => \$lambda,
   'print_cor' => \$print_cor,
   'pearson_genotypes' => \$pearson_genotypes,
   'distance|d=i' => \$distance, 
   'sample_score' => \$sample_score,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);
pod2usage(0) if (not defined $assoc);
open (LOG,">$out.log") or print_OUT("I can not open [ $out.log ] to write to") and exit(1);
print_OUT("FORGE version [ $VERSION ]");
print_OUT("LOG file will be written to [ $out.log ]");

# define distance threshold,
defined $distance or $distance = 20;
print_OUT("Max SNP-to-gene distance allowed [ $distance ] kb");

# $report is use to specify how often report the advance when reading input files
defined $report or $report = 50_000;
# tell if analysis is restricted to a specific chromosome
defined $analysis_chr and print_OUT("Restricting analysis to chromosome [ $analysis_chr ]");
#  defined lambda value for Genomic Control correction
defined $lambda or $lambda = 1;
# tell if user wants to print the *.correlation file 
defined $print_cor and print_OUT("Defined -print_cor: I will print the *.correlation file (it is bulky)");
#set output file if not set already.
defined $out or $out = "gene_based_fisher_v041.OUT";
# tell if user wants to correct p-values by genomic control 
if ($lambda != 1){ print_OUT("SNP p-value will be corrected with lambda = [ $lambda ]");}

open (OUT,">$out") or print_OUT("I can not open [ $out ] to write to") and exit(1);
print OUT  "Ensembl_ID\tHugo_id\tgene_type\tchromosome\tstart\tend\tmin_p\tmin_p_SIDAK\tFORGE\tFORGE_chi-square\tFORGE_df\tn_snps\tn_effective_tests\n";


# i will read the gene_list and i will load data for just this genes to speed up.
if ( not defined $all_genes and not defined @genes and not defined $gene_list){
  $all_genes = 1;
  print_OUT("WARNING: You did not provide an option for the set of genes to be analyzed. I will analyze all genes covered by the SNP association file. Check documentation for options -genes and -gene_list otherwise");
}

if ( not defined $all_genes ) { # in case user want to analyze all genes
  if ( not defined @genes ) { # in case user gave a list of genes in the command line
    print_OUT("Reading Gene List from [ $gene_list ]");
    # read file with gene list and store gene names.
    open( GL, $gene_list ) or print_OUT("I can not open [ $gene_list ]") and exit(1);
    @genes = <GL>;
    chomp(@genes);
    close(GL);
  } else { 
	print_OUT("Read Gene List command line [ @genes  ]"); 
  }
} else {
	print_OUT("Going to analyze all genes on [ $snpmap ] file."); 
}

# Now lets going to read the affy id to rsid mapping. This is used to keep all ids in the
# same nomenclature
my %affy_id = ();
if ( defined $affy_to_rsid ) { # if conversion file is defined
   print_OUT("Reading AFFY to rsID mapping from [ $affy_to_rsid ]");
   open( AFFY, $affy_to_rsid ) or print_OUT("I can not open [ $affy_to_rsid ]") and exit(1);
   while (my $affy = <AFFY>){
      chomp($affy);
      my @b = split(/\t+/,$affy);
      $affy_id{$b[0]} = $b[1];
   }
   close(AFFY);
}

# Read file with genetic association results.
print_OUT("Reading association file: [ $assoc ]");
# create hash to store SNP information
my %assoc_data = ();
open( ASSOC, $assoc ) or print_OUT("I can not open [ $assoc ]") and exit(1);
my $line = 0;
my %header = ();
while (  my $a = <ASSOC> ) {
  $a =~ s/^\s+//;
  $a =~ s/^\t+//;
  $a =~ s/\s+/\t/g;
  # here I get the header line of the file. With the name of the columns I can use the cols
  # SNP and P to extract the information.
  my @data = split( /[\s+\t+]/, $a );
   if ( $line == 0 ){
      %header = %{get_header(\@data)};
      $line++;
      next;
   }
   # In case there is not cols with SNP and P names.
   exists $header{"SNP"} or print_OUT("Not p-value columns available, these are the headers i found [ " . (keys %header) . " ]") and exit(1);
   exists $header{"P"} or print_OUT("Not p-value columns available, these are the headers i found [ " . (keys %header) . " ]") and exit(1);
   # if there is a cols specifying the association test done. Only use result from ADD tests, this is only for compatibility with PLINK 
  if ( exists $header{"TEST"}){
      next if ( $data[$header{"TEST"}] ne "ADD");
   }
  #if (defined $v ) { print $data[$header{"SNP"}]," ", $data[$header{"P"}],"\n"; }
	  # if there is an affy id convert it to rsid.
   if ( defined $affy_to_rsid ) {
      if ($data[$header{"SNP"}] !~ m/^rs/){
         if (exists $affy_id{$data[$header{"SNP"}]}){ $data[$header{"SNP"}] = $affy_id{$data[$header{"SNP"}]};}
      }
   }
   # skip if no P-value or p-value equal NA
   next if ( $data[$header{"P"}] eq "NA");
   next if ( $data[$header{"P"}] eq "");
   
   #generate a pseudo-hash for each snp with the association info
   # correct for genomic control if a lambda > 1 was specified.
  my $A2 = "NA";
  my $A1 = "NA";
  my $OR = "NA";
  if (exists $header{"A2"}){ $A2 = $data[$header{"A2"}];}
  if (exists $header{"A1"}){ $A1 = $data[$header{"A1"}];}
  if (exists $header{"OR"}){ $OR = $data[$header{"OR"}];}
   $assoc_data{ $data[$header{"SNP"}] } = {
                				'pvalue' => 1 ,
						'id' => $data[$header{"SNP"}],
						'or' => $OR,
						'a1' => $A1,
						'a2' => $A2,
        				};
   
   if ($lambda == 1) {
	$assoc_data{ $data[$header{"SNP"}] }->{'pvalue'} = $data[$header{"P"}];
   } elsif ($lambda > 1) {# transform the p-value on a chi-square, correct it by the inflation factor and transform it again on a p-value 
   	$assoc_data{ $data[$header{"SNP"}] }->{ 'pvalue' }= 1-  gsl_cdf_chisq_P( gsl_cdf_chisq_Pinv( $data[$header{"P"}], 1 )/$lambda, 1 );
   } else { 
	print_OUT("\nPlease check the lambda value is correct\n");
	exit(1);
   }
   
   
}
close(ASSOC);

if (scalar keys %assoc_data == 0){ 
	print_OUT("\nNo SNPs with genetic association to used in the analysis\n");
	exit(1);
}
print_OUT("[ " . scalar (keys %assoc_data) . " ] SNPs with association data");

#read snp-to-gene mapping and store in a hash with key equal gene name and value
# an array with the snps in the gene.
my @bim = ();
my @fam = ();
my $ped_map_genotypes;
if (defined $bfile) {
  	# read the bim file with snp information and fam file with sample information
  	@bim = @{ read_bim("$bfile.bim") };
  	
  	@fam = @{ read_fam("$bfile.fam") };
} elsif (defined $ped and defined $map){
	my ($fam_ref,$bim_ref);
	($fam_ref,$ped_map_genotypes,$bim_ref) = read_map_and_ped($ped,$map);
	@fam = @$fam_ref;
	@bim = @$bim_ref;
}
my %bim_ids = ();
my $index = 0;
map {
  $bim_ids{$_->{snp_id}} = $index;
  $index++;
} @bim;
print_OUT("Reading gene to SNP mapping file from [ $snpmap ]");

open( MAP, $snpmap ) or print_OUT("Can not open [ $snpmap ] file") and exit(1);
my %gene = ();
my %snp_to_gene = ();
while ( my $read = <MAP> ) {
   chomp($read);
   # the line is separate in gene info and snps. the section are separated by a tab.
   my ($chr,$start,$end,$ensembl,$hugo,$gene_status,$gene_type,$description,@m) = split(/\t+/,$read);
   my @first_snp_n_fields =  split(/\:/,$m[0]);
  if (4 !=  scalar @first_snp_n_fields){ $description .= splice(@m,0,1); }

  # get all mapped snps within the distance threshold,
   my @mapped_snps = ();
  foreach my $s (@m) {
	my ($id,$pos,$allele,$strand) = split(/\:/,$s);
	next if (not defined $id);
	if (( $pos >= $start) and ($pos <= $end)){ push @mapped_snps, $id; }
	elsif ( ( abs ($pos - $start) <= $distance*1_000 ) or ( abs ($pos - $end) <= $distance*1_000 )) { push @mapped_snps, $id; }
  }
   
   next if (scalar @mapped_snps == 0);
   # get the gene position info
   #check if gene was in the list of genes i want to analyze
   unless ( defined $all_genes ) {
      next unless ( ( grep $_ eq $hugo, @genes ) or ( grep $_ eq $ensembl, @genes ) );
   }
   if (defined $analysis_chr){
	next if ($analysis_chr ne $chr);
   }
   # create a pseudo-hash with the gene info
   $gene{$ensembl} = {
      	'hugo'      => $hugo,
      	'chr'       => $chr,
      	'start'     => $start,
      	'end'       => $end,
      	'gene_type' => $gene_type,
      	'snps'      => [],
      	'minp'      => -9,
      	'genotypes' => null,
      	'geno_mat_rows' => [],
	'pvalues' => [],
	'gene_status' => $gene_status,
	'desc' => $description,
   };

   # go over mapped snps and change convert affy ids to rsid.
   # and make a non-redundant set.
   my %nr_snps = ();
   foreach my $s (@mapped_snps) {
      if ( defined $affy_to_rsid ) {
         if ($s !~ m/^rs/){
            if (exists $affy_id{$s}){ $s = $affy_id{$s};}
         }
      }
	  # exclude snps not in the association file nor in the bim file
	  next unless ( exists $assoc_data{$s} );
	  next unless ( exists $bim_ids{$s});
	  $nr_snps{$s} = "";
   }
   @mapped_snps = keys %nr_snps;
   
   # go over the snps mapped to the gene and check if they are in the map
   # and association files. If so, store the min p-value for the gene.
   # if any of the snps is in the files the remove the gene from the analysis.
   foreach my $s (@mapped_snps) {
	  if (defined $v ){ print_OUT("Mapping [ $s ] to [ $ensembl ]");}
	  next if ( grep $_ eq $s, @{ $gene{$ensembl}->{snps} } );
      push @{ $snp_to_gene{$s} }, $ensembl;
      push @{ $gene{$ensembl}->{snps} }, $s;	
      if ( $gene{$ensembl}->{minp} == -9) {
         $gene{$ensembl}->{minp} = $assoc_data{$s}->{pvalue};
      } elsif ( $assoc_data{$s}->{pvalue} < $gene{$ensembl}->{minp} ) {
         $gene{$ensembl}->{minp} = $assoc_data{$s}->{pvalue};
      }
   }
   # remove gene if none of its snps is in the analysis.
   if ( scalar @{ $gene{$ensembl}->{snps} } == 0 ) { delete( $gene{$ensembl} ) }
   else { 
	if (defined $v ){ print_OUT("Gene $ensembl $hugo included in the analysis with [ ", scalar @{ $gene{$ensembl}->{snps} }, " ] mapped SNPs"); }
   }
   
}
close(MAP);

print_OUT("  '->[ " . scalar (keys %gene) . " ] Genes will be analyzed");
print_OUT("  '->[ " . scalar (keys %snp_to_gene) . " ] SNPs mapped to Genes and with association results will be analyzed");


if (scalar keys %gene == 0){
        print_OUT("No genes mapped");
        exit(1);
}

# start a hash to store the SNP-to-SNP correlation values
my %correlation = ();

# if provided get the SNP-to-SNP correlation values from a tab separated file with 3 cols:snp1 snp2 correlatio_value
if (defined $spearman){
   print_OUT("Reading SNP correlation from [ $spearman ]");
   open( SPRMN, $spearman ) or print_OUT("I cannot open [ $spearman ]") and exit(1);
   while (my $ln = <SPRMN>){
      chomp($ln);
      my @a = split(/\s+/,$ln);
	  # take the square of the correlation if they are r2 and user wants to use r values
      $correlation{$a[0]}{$a[1]} = $a[2];
      $correlation{$a[1]}{$a[0]} = $a[2];
	  # set self correlation to 1
      $correlation{$a[0]}{$a[0]} = 1;
      $correlation{$a[1]}{$a[1]} = 1;
   }
}
# start output file and print its header
print_OUT("Output file will be written to [ $out ]");

# if printing sample level scores.
my $out_fh_sample_score = new IO::File if (defined $sample_score);

if (defined $sample_score){
	print_OUT("   '-> Sample scores printed to [ $out.sample_score ]");
	$out_fh_sample_score->open(">$out.sample_score");
	my $iid = "";
	my $fid = "";
	my $pid = "";
	my $mid = "";
	my $sex = "";
	my $pheno = "";
	for (my $i = 0; $i < scalar @fam; $i++){
		$iid .= " $fam[$i]->{iid}";
		$fid .= " $fam[$i]->{fid}";
		$pid .= " $fam[$i]->{pid}";
		$mid .= " $fam[$i]->{mid}";
		$sex .= " $fam[$i]->{sex}";
		$pheno .= " $fam[$i]->{pheno}";
	}
	print $out_fh_sample_score "IID IID$iid\n";
	print $out_fh_sample_score "FID FID$fid\n";
	print $out_fh_sample_score "PID PID$pid\n";
	print $out_fh_sample_score "MID MID$mid\n";
	print $out_fh_sample_score "SEX SEX$sex\n";
	print $out_fh_sample_score "TRAIT BD$pheno\n";
	
# 	my $rng = PDL::GSL::RNG->new('taus');
# 	$rng->set_seed(time()*$$);
# 	my $a=zeroes(5,5,5)
# 	$rng->get_uniform($a);
# 	print $a;
# 	getc;
}

#start count to report advance 
my $count = 0;
print_OUT("Starting to Calculate gene p-values");

# if there are more than 100 genes change the $report variable in order to report every ~ 10 % of genes.
unless (scalar keys %gene < 100){
	  $report = int((scalar keys %gene)/100 + 0.5)*10;
}
# if user defined a genotypes file. read genotypes and store a genotype matrix (rows: samples, cols: genotypes)for each gene 
if (defined $bfile) {
  print_OUT("Reading genotypes from [ $bfile.bed ]");
  # open genotype file
  my $bed = new IO::File;
  $bed->open("<$bfile.bed") or print_OUT("I can not open binary PLINK file [ $bfile ]") and exit(1);
  binmode($bed); # set file type to binary
  # check if the file is a PLINK file in the proper format by checking the first 3 bytes
  my ($buffer,$n_bytes); 
  my $plink_bfile_signature = "";
  read $bed, $plink_bfile_signature, 3;
  if (unpack("B24",$plink_bfile_signature) ne '011011000001101100000001'){
    print_OUT("Binary file is not in SNP-major format, please check you files\n");
    exit(1);
  } else { print_OUT("Binary file is on SNP-major format"); }
  # calculate how many bytes are needed  to encode a SNP
  # each byte has 8 bits with information for 4 genotypes
  my $N_bytes_to_encode_snp = (scalar @fam)/4; # four genotypes per byte
  # if not exact round it up
  if (($N_bytes_to_encode_snp - int($N_bytes_to_encode_snp)) != 0  ){ $N_bytes_to_encode_snp = int($N_bytes_to_encode_snp) + 1;}
  # loop over all genes and extract the genotypes of the SNPs
  foreach my $gn (keys %gene){
    # this will store the genotypes
    my $matrix;
    # loop over the snps mapped to the gene
    foreach my $mapped_snp (@{$gene{$gn}->{snps}}){
      # skip if it does not have association information
      next if (not exists $assoc_data{ $bim[$bim_ids{$mapped_snp}]->{snp_id} } );
      
      if (defined $v){ print_OUT("Adding SNP [  $bim[ $bim_ids{$mapped_snp} ]->{snp_id}  ] to genotypes of $gn"); }
      # because we know the index of the SNP in the genotype file we know on which byte its information starts
      my $snp_byte_start = $N_bytes_to_encode_snp*$bim_ids{$mapped_snp};
      # here i extract the actual genotypes
      my @snp_genotypes = @{ extract_binary_genotypes(scalar @fam,$N_bytes_to_encode_snp,$snp_byte_start,$bed) };
      # store the genotypes.
      # if a snp does not use the 8 bits of a byte the rest of the bits are fill with missing values
      # here i extract the number of genotypes corresponding to the number of samples
      push @{ $matrix}, [@snp_genotypes[0..scalar @fam - 1]];
      # add snp id to matrix row names
      push @{ $gene{$gn}->{geno_mat_rows} }, $bim[ $bim_ids{$mapped_snp} ]->{snp_id};
      # store the p-value of the snp
      push @{ $gene{$gn}->{pvalues} }, $assoc_data{ $bim[ $bim_ids{$mapped_snp} ]->{snp_id} }->{pvalue}; 
    }
    # generate the genotype matrix as a PDL piddle
    $gene{$gn}->{genotypes} = pdl $matrix;
    # calculate gene p-values
    &gene_pvalue($gn);
    # delete the gene's data to keep memory usage low
    delete($gene{$gn});
    $count++;
    &report_advance($count,$report,"Genes");	
 }
} elsif (defined $ped and defined $map){
	my $genotypes = pdl @{ $ped_map_genotypes };
	foreach (my $index_snp = 0; $index_snp <  scalar @bim; $index_snp++){
		# store SNP genotypes only if it has association data
		if (exists $assoc_data{ ${ $bim[$index_snp] }{snp_id} } ){
			# for every gene mapped to this SNP, push inside the genotype matrix this SNP genotypes.
			foreach my $gn (@{ $snp_to_gene{${ $bim[$index_snp] }{snp_id} } }){
				if (defined $v){ print_OUT("Adding SNP [  ${ $bim[$index_snp] }{snp_id}  ] to genotypes of $gn"); }
				next if ( grep $_ eq ${ $bim[$index_snp] }{snp_id} , @{ $gene{$gn}->{geno_mat_rows} } );
				# add the genotypes to the genotype matrix
				$gene{$gn}->{genotypes} = $gene{$gn}->{genotypes}->glue(1,$genotypes(,$index_snp));
				push @{ $gene{$gn}->{geno_mat_rows} }, ${ $bim[$index_snp] }{snp_id};
				push @{ $gene{$gn}->{pvalues} }, $assoc_data{ ${ $bim[$index_snp] }{snp_id} }->{pvalue};
				if (scalar @{ $gene{$gn}->{snps} } == scalar @{ $gene{$gn}->{geno_mat_rows} }){
				  &gene_pvalue($gn);
				  delete($gene{$gn});
				  $count++;# if there are more than 100 genes change the $report variable in order to report every ~ 10 % of genes.
				  unless (scalar keys %gene < 100){
					  my $report = int((scalar keys %gene)/100 + 0.5)*10;
					  &report_advance($count,$report,"GENES");
				  }
				}
			}
		}
	}
}else {
  print_OUT("WARNING: Gene p-values will be calculated with the precomputed correlation only. If correlation for some SNPs pairs are missing you may get wrong results, please check your inputs for completeness");
}
 

# loop over all genes and calculate the gene p-values

$out_fh_sample_score->close() if (defined $sample_score);
# if the user want to get the correlation values print the *.correlation file
if (defined $print_cor){
	  open (COR,">$out.correlation") or print_OUT("Cannot open [ $out.correlation ] to write to") and exit(1);
	  foreach my $snp1 (keys %correlation) {
		 foreach my $snp2 (keys %{$correlation{$snp1}}  ) {
			next if ($snp1 eq $snp2);
			printf COR ("$snp1 $snp2 %.3f\n",abs($correlation{ $snp1 }{ $snp2 }));
			delete($correlation{ $snp1 }{ $snp2 });
		 }
	  }
	  close(COR);
}
print_OUT("Well Done!!");
exit(0);
sub extract_binary_genotypes {
  my $n_genotypes = shift; # number of genotypes per SNP
  my $bytes_per_snp = shift; # number of bytes needed to code the SNP
  my $byte_position = shift; # starting byte position for this SNP
  my $FH = shift; # file handle for the genotype file
  $FH->seek(3 + $byte_position,SEEK_SET); # re-set the file-handle to position start position of the SNP of interest
  my $buffer = ""; # this will store the information read
  my $n_bytes = read $FH, $buffer, $bytes_per_snp; # read the genotypes
  my $data_size = $bytes_per_snp*8; # the amount of data to extract is 8 bits per byte
  my $bin_data = unpack("B$data_size",$buffer); 
  my @bits = ( $bin_data =~ m/\d{8}/g );
  my @genotypes = ();
  foreach my $b (@bits){
    $b = reverse($b); # for some odd reason PLINK stores the genotypes in reverse order
    push @genotypes, @{ get_genotypes($b)};# transform each byte on genotypes
  }
  return(\@genotypes);
}

sub get_genotypes {
  my $b = shift; # a byte
  my @back = ();
  my @genotypes = ( $b =~ m/\d{2}/g ); # extract a pair of number = a genotype
  foreach my $geno (@genotypes){
    if    ( $geno eq '00' ) {  # homozygous 1/1
      push @back, '1';
    } elsif ( $geno eq '11' ) { # -- other homozygous 2/2
      push @back, '3';
    } elsif ( $geno eq '01' ) { # -- heterozygous 1/2
      push @back, '2';
    } elsif ( $geno eq '10' ) { # -- missing genotype 0/0
      push @back, '0';
    } else { print_OUT("This genotype is not recognize [ $geno ]"); }    # genotype not recognize
  }
  return(\@back);
}

sub gene_pvalue {
	my $gn = shift;
   if (defined $v){ print_OUT("____ $gn ____"); }
   if (not defined $bfile){ @{ $gene{$gn}->{geno_mat_rows}} = @{ $gene{$gn}->{snps}}; }
   my $n_snps = scalar @{ $gene{$gn}->{geno_mat_rows}};
   next if ($n_snps == 0);
   # if the gene has just 1 SNP we make that SNP's p value the gene p-value under all methods
	  if ($n_snps == 1){
			if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\tNA\tNA1\t1\n",$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp}); }
			printf OUT ("$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\tNA\tNA1\t1\n",$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp});
			delete($gene{$gn});
			return(); }
   #if the user defined a genotype file then we need to get the number of SNPs in the gene.
   # if the user did not defined neither a genotype nor a file with the SNP-SNP correlations exit the program
   unless (defined $bfile or defined $spearman or defined $ped) {
    print_OUT("You MUST specify file with the SNP-SNP correlations or a genotype file to calculate them by myself\n\n");
    exit(1);
  }
   
   if (defined $v){ print_OUT("Calculating correlation matrix for $n_snps SNPs with association data"); }
   # initialize correlation matrix with zeroes
   my $cor = zeroes $n_snps,$n_snps;
   # set equal weight to all SNPs
   # initialize the second term for the variance of the statistics  for the Makambi method
   # loop over all pairs of SNPs and get their correlation values
   for ( my $i = 0; $i < $n_snps; $i++){
      for ( my $p = $i; $p < $n_snps; $p++){
	 # if the correlation is known
	 if (defined $v){ print_OUT("working in $i:$gene{$gn}->{geno_mat_rows}->[$i]; $p:$gene{$gn}->{geno_mat_rows}->[$p]"); }
	  # if is it a self correlation then don't calculate anything
	 if ($gene{$gn}->{geno_mat_rows}->[$p] eq  $gene{$gn}->{geno_mat_rows}->[$i]){
			set $cor, $i, $p, 1;
            set $cor, $p, $i, 1;
	        $correlation{ ${ $gene{$gn}->{geno_mat_rows} }[$i] }{${ $gene{$gn}->{geno_mat_rows} }[$p]} = 1;
			$correlation{ ${ $gene{$gn}->{geno_mat_rows} }[$p] }{${ $gene{$gn}->{geno_mat_rows} }[$i]} = 1;
			next;
	 }
	 # if the correlation is stored in the %correlation hash then fetch the value
	 if (exists $correlation{ $gene{$gn}->{geno_mat_rows}->[$p] }{ $gene{$gn}->{geno_mat_rows}->[$i] } ){
	    	my $c = $correlation{ $gene{$gn}->{geno_mat_rows}->[$p] }{ $gene{$gn}->{geno_mat_rows}->[$i] };
	    	set $cor, $i, $p, $c;
		set $cor, $p, $i, $c;
      	} else { 
		# if the user did not provided a genotype file and tell that this pair of SNPs does not have a correlation value
		unless ((defined $bfile) or (defined $ped)){print_OUT(" WARNING: this pair of snps does not have a correlation value: [$gn => $gene{$gn}->{snps}->[$p] ] and [ $gene{$gn}->{snps}->[$i] ]"); }
		# compute the correlation
		my $c = "";
		if (defined $pearson_genotypes) { # if requested used the Pearson correlation for genotypes
		  $c = pearson_corr_genotypes([list $gene{$gn}->{genotypes}->(,$i)],[list $gene{$gn}->{genotypes}->(,$p)]);
		} else { # by default use the spearman rank correlation
		  my $d = Statistics::RankCorrelation->new([list $gene{$gn}->{genotypes}->(,$i)],[list $gene{$gn}->{genotypes}->(,$p)]);
		  $c = $d->spearman;
		}
		set $cor, $i,$p, $c;
		set $cor, $p,$i, $c;
		# increase the second term of the Makambi method
		# store the correlation value in case is needed later.
		$correlation{ ${ $gene{$gn}->{geno_mat_rows} }[$i] }{${ $gene{$gn}->{geno_mat_rows} }[$p]} = $c;
		$correlation{ ${ $gene{$gn}->{geno_mat_rows} }[$p] }{${ $gene{$gn}->{geno_mat_rows} }[$i]} = $c;
      }
    }
  }
   if (defined $v){ print_OUT($cor);}
   if (defined $v){ print_OUT("Calculating effective number of tests: "); }
   # calculate number of effective tests by the Gao ($k) and Galwey ($Meff_galwey) method.	  
   my ($k,$Meff_galwey) = number_effective_tests(\$cor);

   # get the log of the SNP p-value
   $gene{$gn}->{pvalues} = pdl @{ $gene{$gn}->{pvalues} };
   # calculate the chi-square statistics for the Makambi method and its p-value
   my $weight = ones $n_snps;
    $weight *= 1/$n_snps;
    $weight /= $weight->sumover;
    if (defined $v){ print_OUT("Weigth = [ $weight ]"); }
    my ($forge_chi_stat,$forge_df) = get_makambi_chi_square_and_df($cor,$weight,$gene{$gn}->{pvalues});

   my $fisher_p_value =  1 - gsl_cdf_chisq_P($forge_chi_stat, $forge_df );
   
   # print out the results
   my $sidak = 1-(1-$gene{$gn}->{minp})**$k;
   if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t%0.3e\t%0.5f\t%0.5f\t%2d\t%3d || $weight\t$gene{$gn}->{pvalues}->log\t@{ $gene{$gn}->{geno_mat_rows} }\n",$gene{$gn}->{minp},$sidak,$fisher_p_value,$forge_chi_stat,$forge_df,$Meff_galwey,$n_snps,$k); }
   printf OUT ("$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t%2d\t%3d\n",$gene{$gn}->{minp},$sidak,$fisher_p_value,$forge_chi_stat,$forge_df,$n_snps,$k);

  if (defined $sample_score){
    my @snp_info = map { $assoc_data{$_}; } @{ $gene{$gn}->{geno_mat_rows} } ;
    # alleles have been coded as 1 : homozygote 1/1, 2 heterozygous, 3: homozygote 2/2 and 0: missing.
    # risk dosage is calculated by: number of risk alleles * odd-ratio. odd ratio for the risk allele
    my ($n_samples,$n_snps) = $gene{$gn}->{genotypes}->dims;
    my $out_line = "$gn $gene{$gn}->{hugo}";
    for (my $person = 0; $person < $n_samples; $person++){
	    my @genotypes = $gene{$gn}->{genotypes}->($person,)->list;
	    my @tmp = ();
	    foreach my $g (@genotypes) {
		    if ($g == 1) { push @tmp, 2; }
		    if ($g == 2) { push @tmp, 1; }
		    if ($g == 3) { push @tmp, 0; }
		    if ($g == 0) { push @tmp, 0; }
	    }
	    my $sample_w = pdl @tmp;
	    my @tmp_or = ();
	    for my $snp (@snp_info){
		    if ($snp->{or} < 1 ){ push @tmp_or, 1/$snp->{or};}
		    else { push @tmp_or, $snp->{or}; }
	    }
	    $sample_w *= pdl @tmp_or;
	    $sample_w += 1/$n_snps;
	    $sample_w /= $sample_w->sumover;
	    my ($person_chi_stat,$person_df) = get_makambi_chi_square_and_df($cor,$sample_w,$gene{$gn}->{pvalues});
	    my $fisher_p_value_galwey = 1 - gsl_cdf_chisq_P( $person_chi_stat, $person_df );
	    $out_line .= " $fisher_p_value_galwey";
    }
    print $out_fh_sample_score "$out_line\n";
  }
}

# this subroutine applies the makambi method to combine p-values from correlated test
# the method receives a correlation matrix, a set of weights for each statistic and a set of probabilies to combine.
# it returns the chi-square and the degrees of freedom of the test

sub get_makambi_chi_square_and_df {
  my $cor = shift; # a pdl matrix with the varoable correlations
  my $w = shift; # a pdl vector with the weights;
  # make sure weiths sum one
  $w /= $w->sumover;
  my $pvalues = shift; # a pdl vector with the p-values to be combined
  # calculate the correlation matrix before the applying the weights
  my $COR_MAT = (3.25*abs($cor) + 0.75*abs($cor)**2);
  my $second = $COR_MAT*$w*$w->transpose; # apply the weights 
  ($second->diagonal(0,1)) .= 0; # set the diagonal to 0
  my $varMf_m = 4*sumover($w**2) + $second->flat->sumover; # calculate the variance of the test statistics
  my $df = 8/$varMf_m; # the degrees of freedom of the test statistic
  my $chi_stat = sumover(-2 * $pvalues->log * $w); # and the chi-square for the combine probability
  $chi_stat = ( $chi_stat/2 ) * $df;
  return ($chi_stat,$df);
}

# this subroutine take an array and return a hash were every element of the line is
# a key and the value is the index in the array
sub get_header {
   my $in = shift;
   my %back = ();
   for (my $i = 0;$i< scalar @$in; $i++){
      $back{$$in[$i]} = $i;
   }
   return(\%back);
}

# this subroutine calculate the number of effective test by the Galwey and Gao method.
sub number_effective_tests {
   my $mat = shift;
   # calculate the eigen value of the correlation matrix
   my $eigens = eigens ${$mat};
   # normalize the values
   my $eigens_norm =  pdl sort { $b <=> $a } list ($eigens/(sumover $eigens) );
   # calculate number of effective test per Gao et al Gao, X., Becker, L. C., Becker, D. M., Starmer, J. D. & Province, M. A. Avoiding the high Bonferroni penalty in genome-wide association studies. Genet. Epidemiol. (2009).
   # the methos simply count how many eigen values are needed to explain 99.5 % of the variance
   my $simpleM = 0;
   for( $simpleM = 0; $simpleM < scalar list $eigens_norm; $simpleM++){
      my $sum = sumover $eigens_norm->slice("0:$simpleM");
      if ($sum >= 0.995){
         $simpleM++;
         last;
      }
   }
   # calculate number of effective test by Galwey, N. W. A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests. Genet. Epidemiol. (2009).
   # this method calculate the square of sum of the square-root of the eigen values divided by the sum.
   my $numerator = 0;
   my $denominator = 0;
   foreach my $e (list $eigens) {
	  if ($e > 0){
			$numerator += sqrt($e);
			$denominator += $e;
	  }
   }
   my $Meff_galwey = ($numerator**2)/$denominator;
   if (defined $v){ print_OUT(" simpleM_Gao = $simpleM; Meff_galwey=$Meff_galwey $numerator $denominator"); }
   return($simpleM,$Meff_galwey);
}

# this subroutine receives a list of SNPs and return the log of the their p-values
sub get_logs {
   my $snps = shift;
   my @stats = ();
   if (defined $v) { print_OUT("____ ASSOCIATIONS ____ \n  [ ",scalar @$snps," ] SNPs"); }
   foreach my $s (@{$snps}){
      my $p = $assoc_data{$s}->{pvalue};
      next if ($p eq "NA");
      if ($p == 0){ push @stats ,1;}
      else {
         push @stats ,log($p);
      }
      if (defined $v) {  print_OUT("\t$s\t$p"); }
   }
   if (defined $v) { print_OUT("______________________"); }
   return(\@stats);
}

# this subroutine print the advance of of something
# it receives 3 parameter: index, rep and tag. index is the advance so far, report how often to report and tag is name for the printing
# it will print a report when index if divisible by rep
sub report_advance {
	  my ($index,$rep,$tag) = @_;
   if (( $index/$rep - int($index/$rep)) == 0) {print_OUT(" '->Done with [ $index ] $tag"); }
}

# this subroutine read the fam file and stores the information in an array
# the elements of the array are pseudo hashes with all the sample's information
sub read_fam {
   my $fam = shift;
   print_OUT("Reading samples info from [ $fam ]");
   open( FAM, $fam ) or print_OUT("Cannot open [ $fam ] file") and exit(1);
   my @back = ();
   while ( my $s = <FAM> ) {
      my @data = split( /\s+/, $s );
      push @back,
        {
         'iid'   => $data[0],
         'fid'   => $data[1],
         'mid'   => $data[2],
         'pid'   => $data[3],
         'sex'   => $data[4],
         'pheno' => $data[5],
        };
   }
   print_OUT("[ " . scalar @back . " ] samples read");
   return ( \@back );
}

# this subroutine read the bim file and store information about the SNPs
# each element of the array returned is a pseudo hash with the SNP information
sub read_bim {
   my $bim = shift;
   print_OUT("Reading SNPs info from [ $bim ]");
   open( BIM, $bim ) or print_OUT("Cannot open [ $bim ] file") and exit(1);
   my @back = ();
   while ( my $snp = <BIM> ) {
      chomp($snp);
      my @data = split( /\t+/, $snp );
	  # if an affy to rsid mapping was provided change the ids
      if ( defined $affy_to_rsid ) {
         if ($data[1] !~ m/^rs/){
            if (exists $affy_id{$data[1]}){ $data[1] = $affy_id{$data[1]};}
         }
      }
      push @back,
        {
         'snp_id' => $data[1],
         'chr'    => $data[0],
         'cm'     => $data[2],
         'pos'    => $data[3],
         'a2'     => $data[4],
         'a1'     => $data[5],
        };
   }
   print_OUT("[ " .  scalar @back . " ] SNPs on BED file");
   return ( \@back );
}

sub read_map {
   my $map = shift;
   print_OUT("Reading SNPs info from [ $map ]");
   open( MAP, $map ) or print_OUT("Cannot open [ $map ] file") and exit(1);
   my @back = ();
   while ( my $snp = <MAP> ) {
      chomp($snp);
      my @data = split( /\t+/, $snp );
	  # if an affy to rsid mapping was provided change the ids
      if ( defined $affy_to_rsid ) {
         if ($data[1] !~ m/^rs/){
            if (exists $affy_id{$data[1]}){ $data[1] = $affy_id{$data[1]};}
         }
      }
      push @back, {
         'snp_id' => $data[1],
         'chr'    => $data[0],
         'cm'     => $data[2],
         'pos'    => $data[3],
         'a2'     => 0,
         'a1'     => 0,
        };
   }
   print_OUT("  '->[ " . scalar @back . " ] SNPs on PED file");
   return ( \@back );
}

sub read_map_and_ped {
	my $file = shift;
	my $map = shift;
	my @bim = @{ read_map($map) };
	print_OUT("Reading Genotypes from [ $file ]");
	open(PED,$file) or print_OUT("Cannot open [ $file ] file") or exit(1);
	my @back_fam = ();
	my @geno_matrix= ();
	while (my $sample = <PED>){
		chomp($sample);
		my ($iid,$fid,$mid,$pid,$sex,$pheno,@genotypes) = split(/\s+/,$sample);
		my $snp_counter = 0;
		for (my $g = 0; $g < scalar @genotypes; $g+=2){
			$geno_matrix[$snp_counter]->{alleles}->{$genotypes[$g]}++;
			$geno_matrix[$snp_counter]->{alleles}->{$genotypes[$g+1]}++;
			push @{ $geno_matrix[$snp_counter]->{genotypes} }, "$genotypes[$g]$genotypes[$g+1]";
			$snp_counter++; 
		}
		push @back_fam,{
			'iid'   => $iid,
			'fid'   => $fid,
			'mid'   => $mid,
			'pid'   => $pid,
			'sex'   => $sex,
			'pheno' => $pheno,
			};
	}
	my @back_genotypes = ();
	my $snp_counter = 0;
	foreach my $snp (@geno_matrix){
		my @alleles = sort { $snp->{alleles}->{$b} <=> $snp->{alleles}->{$a} } keys %{ $snp->{alleles} };
		($snp->{alleles}->{major},$snp->{alleles}->{minor},$snp->{alleles}->{missing}) = 0;
		for (my $i = 0; $i < scalar @alleles; $i++) {
			if ($alleles[$i] == 0) { 
				$snp->{alleles}->{missing} =$alleles[$i]; 
			} else {
				if ( $snp->{alleles}->{major} == 0){
					$snp->{alleles}->{major}= $alleles[$i];
				} else {
					$snp->{alleles}->{minor}= $alleles[$i];
				}
			} 
			
		}
		
		$bim[$snp_counter]->{a1} = $snp->{alleles}->{minor};
		$bim[$snp_counter]->{a2} = $snp->{alleles}->{major};
		foreach my $g (@{ $snp->{genotypes} }){
			my $major_homo = "$snp->{alleles}->{major}$snp->{alleles}->{major}";
			my $minor_homo = "$snp->{alleles}->{minor}$snp->{alleles}->{minor}";
			my $hetero1 = "$snp->{alleles}->{major}$snp->{alleles}->{minor}";
			my $hetero2 = "$snp->{alleles}->{minor}$snp->{alleles}->{major}";
			my $missing = "$snp->{alleles}->{missing}$snp->{alleles}->{missing}";
			my $recoded = 0;
			# homozygous major allele
			if ($g == $major_homo) {
				$recoded = 3;
			} elsif ($g == $minor_homo){ # homozygous minor allele
				$recoded = 1;
			} elsif (($g == $hetero1) or ($g == $hetero2)){ # heterozygous
				$recoded = 2
			} elsif ($g == $missing) { # missing
				$recoded = 0;
			} else { print_OUT(" COULD NOT RECOGNIZE THIS GENOTYPE >$g<\n" . Dumper($snp->{alleles}) . "");
				exit(1);
			}
			push @{$back_genotypes[$snp_counter]}, $recoded;
		}
		$snp_counter++;
	}
	print_OUT("  '-> [ " . scalar @back_fam . " ] samples");
	print_OUT("  '-> [ " . scalar @back_genotypes . " ] SNPs");
	return(\@back_fam, \@back_genotypes,\@bim);
}


sub pearson_corr_genotypes {
  # implemented as in S. Wellek, A. Ziegler, Hum Hered 67, 128 (2009).
  # genotypes must be coded as 1,2 and 3,any other coding will be use. missing genotypes can be set to anything different of 1, 2 or 3.
  # this method will fail is more than 2 out of the four homozygoues haplotypes have counts 0. In such case i recommend to use the default option of the spearman correlation. 
  my $snp1 = shift;
  my $snp2 = shift;
  # initialize to 0 all values for the haplotype cells
  my @cells = (0,0,0,0,0,0,0,0,0);
  # count how many observation of every haplotype are
  for (my $i = 0; $i < scalar @$snp1; $i++){
   $cells[8]++ if ((${$snp1}[$i] == 1) and (${$snp2}[$i] == 1));	
   $cells[7]++ if ((${$snp1}[$i] == 2) and (${$snp2}[$i] == 1));	
   $cells[6]++ if ((${$snp1}[$i] == 3) and (${$snp2}[$i] == 1));	
   $cells[5]++ if ((${$snp1}[$i] == 1) and (${$snp2}[$i] == 2));	
   $cells[4]++ if ((${$snp1}[$i] == 2) and (${$snp2}[$i] == 2));	
   $cells[3]++ if ((${$snp1}[$i] == 3) and (${$snp2}[$i] == 2));	
   $cells[2]++ if ((${$snp1}[$i] == 1) and (${$snp2}[$i] == 3));	
   $cells[1]++ if ((${$snp1}[$i] == 2) and (${$snp2}[$i] == 3));	
   $cells[0]++ if ((${$snp1}[$i] == 3) and (${$snp2}[$i] == 3));	
  }
  if (defined $v) {
	print_OUT("  AA Aa aa");
	print_OUT("BB $cells[0]  $cells[1]  $cells[2]");
	print_OUT("Bb $cells[3]  $cells[4]  $cells[5]");
	print_OUT("bb $cells[6]  $cells[7]  $cells[8]");
  }
  my $h1 = sqrt($cells[0]);
  my $h2 = sqrt($cells[2]);
  my $h3 = sqrt($cells[6]);
  my $h4 = sqrt($cells[8]);
  my $corr = 2*($h1 * $h4 - $h2 * $h3)/sqrt(4*($h1 + $h2)*($h3 + $h4)*($h1 + $h3)*($h2 + $h4));
  return($corr);
}

sub print_OUT {
  my $string = shift;
  print scalar localtime(), "\t$string\n";
  print LOG scalar localtime(), "\t$string\n";
}

__END__

=head1 NAME

 Running network analysis by greedy search

=head1 SYNOPSIS

script [options]

 	-h, --help		print help message
 	-m, --man		print complete documentation
 	-report			how often to report advance
 	-verbose, -v		useful for debugging
 	
 	-ped			Genotype file in PLINK PED format
 	-map			SNP info file in PLINK MAP format
 	-bfile			Files in PLINK binary format, corresping file with extension *.bed, *.bim and *.fam. 
 	-out, -o		Name of the output file
 	-assoc, -a		SNP association file, columns header is necessary and at leat columns with SNP and P names must be present
 	-gene_list, -g		Only analyse genes on the file provided
 	-genes			Provide some gene name in command line, only these genes will be analyzed
 	-all_genes		Analyze all genes in the snp-to-gene mapping file
 	-chr			Anlyze a specific chromosome  
 	-snpmap, -m		Snp-to-gene mapping file
 	-correlation, -cor	SNP-SNP correlation file
 	-affy_to_rsid		Affy id to rs id mapping file
 	-lambda			lambda value to correct SNP statistics, if genomic control is used
 	-print_cor		print out SNP-SNP correlations
 	-pearson_genotypes	Calculate SNP-SNP correlation with a pearson correlation for categorical variables
 	-distance, -d		Max SNP-to-gene distance allowed (in kb) 


=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-report>

how often to report advance. Provide an integer X and the program will report adnvance after X networks are analyzed.

=item B<-verbose, -v>

useful for debugging

=item B<-ped>

Genotype file in PLINK PED format. See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml for details on data format

=item B<-map>

SNP info file in PLINK MAP format

=item B<-bfile>

Files in PLINK binary format, corresping file with extension *.bed, *.bim and *.fam. 

=item B<-out, -o>

Name of the output file. Output file will be tab separated with the following columns: Ensembl_ID,Hugo id,gene_type,chromosome,start,end,min p-value in the gene, Sidak corrected minimum p-value, FORGE p-value, FORGE chi-square, FORGE degrees of freedom, Eigen value ratio method gene p-value, EVR chi-square, EVR degrees freedom, number of snps in gene, number effective tests(Gao et al;PDMID:19434714)
An example would look like:ENSG00000162594	IL23R	protein_coding	1	67632083	67725662	3.491e-02	4.132e-01	2.128e-01	8.42305	6.04941	1.119e-01	28.17326	10.116917	 15

=item B<-assoc, -a>

SNP association file, columns header is necessary and at leat columns with SNP and P names must be present. 

=item B<-gene_list, -g>

Only analyse genes on the file provided. Ids must match either Ensembl Ids or Hugo ids present in the SNP-to-gene mapping file

=item B<-genes>

Provide some gene name in command line, only these genes will be analyzed. Use like -genes IL23R

=item B<-all_genes>

Analyze all genes in the snp-to-gene mapping file

=item B<-chr>

Anlyze a specific chromosome, Chromosome code must match that of the SNP-to-gene mapping file

=item B<-snpmap, -m>

SNP-to-gene mapping file. Format is tab separated with fields:chromosome,start,end,Ensembl id, hugo id, description, SNP1, SNP2, ..., SNPN. 
Each SNP has 4 fields separated by colons: id,position,alleles and strand. This may look overcomplicated but allows to map any kind of variation and store its information without having to use additional files. An example look like:
 
1       67632083        67725662        ENSG00000162594 IL23R   KNOWN   protein_coding  Interleukin-23 receptor Precursor (IL-23R) [Source:UniProtKB/Swiss-Prot;Acc:Q5VWK5]             rs7538938:67132262:T/C:1        rs72669476:67132582:C/T:1       rs72019237:67132863:C/-:1       rs61197134:67132871:C/-:1       rs11208941:67133111:T/C:1 

=item B<-correlation, -cor>

SNP-SNP correlation file. Space separated file with 3 columns, first 2 the SNP ids the the 3th the correlation between them. Like:

rs6660226 rs11209018 0.089

=item B<-affy_to_rsid>

Affy id to rs id mapping file. Tab separated file, looks like:
SNP_A-8389091	rs7593668

=item B<-lambda>

lambda value to correct SNP statistics, if genomic control is used

=item B<-print_cor>

print out SNP-SNP correlations

=item B<-pearson_genotypes>

Calculate SNP-SNP correlation with a pearson correlation for categorical variables

=item B<-distance, -d>

Max SNP-to-gene distance allowed (in kb) 



=back

=head1 DESCRIPTION

TODO


=cut
