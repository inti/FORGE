#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use PDL;
use PDL::Matrix;
use PDL::GSL::CDF;
use PDL::Primitive;
use PDL::NiceSlice;
use PDL::Stats::Basic;
use PDL::Bad;
use IO::File;
use IO::Seekable;
use Fcntl;
use Data::Dumper;
use Pod::Usage;

my $VERSION = "0.9.5.6";

our ( $help, $man, $out, $snpmap, $bfile, $assoc, $gene_list,
    @genes, $all_genes, $analysis_chr, $report, $spearman,
    $affy_to_rsid, @weights_file, $w_header, $v, $lambda,
    $print_cor, $pearson_genotypes,$distance, $sample_score,
    $ped, $map, $ox_gprobs,$sample_score_self, $w_maf,
    $ss_mean, $no_forge, $gc_correction,$g_prob_threshold
);

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
   'gc_correction' => \$gc_correction,
   'print_cor' => \$print_cor,
   'pearson_genotypes' => \$pearson_genotypes,
   'distance|d=i' => \$distance, 
   'sample_score' => \$sample_score,
   'weights|w=s' => \@weights_file,
   'w_header' => \$w_header,
   'ox_gprobs=s' => \$ox_gprobs,
   'g_prob_threshold=f' => \$g_prob_threshold,
   'weight_by_maf|w_maf' => \$w_maf,
   'ss_mean' => \$ss_mean,
   'no_forge' => \$no_forge,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);
pod2usage(0) if (not defined $assoc);
open (LOG,">$out.log") or print_OUT("I can not open [ $out.log ] to write to") and exit(1);
print_OUT("FORGE version [ $VERSION ]. See http://github.com/inti for updates");
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
# defined threshold value for genotype probabilities
defined $g_prob_threshold or $g_prob_threshold = 1.0;
# tell if user wants to print the *.correlation file 
defined $print_cor and print_OUT("Defined -print_cor: I will print the *.correlation file (it is bulky)");
#set output file if not set already.
defined $out or $out = "gene_based_fisher_v041.OUT";
# tell if user wants to correct p-values by genomic control 
if ($lambda != 1){ print_OUT("SNP p-value will be corrected with lambda = [ $lambda ]");}

my $geno_probs = undef;
if (defined $ox_gprobs) {
    $geno_probs = $ox_gprobs ;
    print_OUT("Genotype probabilities in OXFORD format will be read from [ $geno_probs ]");
}

my ($gprobs, $gprobs_index);
if (defined $geno_probs){
    $gprobs = IO::File->new();
    $gprobs_index = IO::File->new();
    
#    if (not -e "$geno_probs.idx") {
    print_OUT("   '-> Making index for genotype probabilities in [ $geno_probs.idx ] file");
    $gprobs->open("<$geno_probs") or print_OUT("I can not open binary PLINK file [ $geno_probs ]") and exit(1);
    $gprobs_index->open("+>$geno_probs.idx") or print_OUT("Can't open $geno_probs.idx for read/write: $!\n");
    build_index(*$gprobs, *$gprobs_index);
        
=head    } else {
        print_OUT("   '-> Found [ $geno_probs.idx ] file, will use it to read the genotype probabilities. Please delete if you do not want to use this file.");
        $gprobs->open("<$geno_probs") or print_OUT("I can not open binary PLINK file [ $geno_probs ]") and exit(1);
        $gprobs_index->open("<$geno_probs.idx") or print_OUT("Can't open $geno_probs.idx for read/write: $!\n") and exit(1);
    }
=cut
} 

unless (defined $no_forge){
	open (OUT,">$out") or print_OUT("I can not open [ $out ] to write to") and exit(1);
	print OUT  "Ensembl_ID\tHugo_id\tgene_type\tchromosome\tstart\tend\tmin_p\tmin_p_SIDAK\tFORGE\tFORGE_chi-square\tFORGE_df\tn_snps\tn_effective_tests\n";
}

# i will read the gene_list and i will load data for just this genes to speed up.
if ( not defined $all_genes and not defined @genes and not defined $gene_list){
  $all_genes = 1;
  print_OUT("Note: You did not provide an option for the set of genes to be analyzed. I will analyze all genes covered by the SNP association file. Check documentation for options -genes and -gene_list otherwise");
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
my $assoc_chis = [] if (defined $gc_correction);

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
  my $A2 = "NA";
  my $A1 = "NA";
  my $OR = 1;
  my $BETA = undef;
  my $R2 = undef;
  my $STAT = undef;
  if (exists $header{"A2"}){ $A2 = $data[$header{"A2"}];}
  if (exists $header{"A1"}){ $A1 = $data[$header{"A1"}];}
  if (exists $header{"OR"}){ $OR = $data[$header{"OR"}];}
  if (exists $header{"BETA"}){ $BETA = $data[$header{"BETA"}];}
  if (exists $header{"R2"}){ $R2 = $data[$header{"R2"}];}
  if (exists $header{"STAT"}){ $STAT = $data[$header{"STAT"}];}
  my $effect = undef;
  if (exists $header{"OR"}){
    $effect = "or";
  } elsif (exists $header{"BETA"}){
        $effect = "beta";
  } elsif (exists $header{"STAT"}){
        $effect = "stat";
  }
  # do some checking
  if (defined $sample_score){
    # if there are not direction of effect defined quite analysis and spit and error.
    if (not defined $effect){
        print_OUT("ERROR: No effect size measure or direction of effect provided. I can not perform Sample Score analysis");
        print_OUT("ERROR: Please provide an odd-ratio, beta or regression coefficient value under the header OR, BETA or STAT, respectively");
        exit(1);
    } 
  }
  $assoc_data{ $data[$header{"SNP"}] } = {
                                        'pvalue' => 1 ,
                                        'id' => $data[$header{"SNP"}],
                                        'a1' => $A1,
                                        'a2' => $A2,
                                        'or' => $OR,
                                        'beta' => $BETA,
                                        'stat' => $STAT,
                                        'r2' => $R2,
                                        'effect_size_measure' => $effect,
        				};
   # correct for genomic control if a lambda > 1 was specified.   
   if ($lambda == 1) {
   	push @{ $assoc_chis }, $data[$header{"P"}];
	$assoc_data{ $data[$header{"SNP"}] }->{'pvalue'} = $data[$header{"P"}];
   } elsif ($lambda > 1) {# transform the p-value on a chi-square, correct it by the inflation factor and transform it again on a p-value 
   	$assoc_data{ $data[$header{"SNP"}] }->{ 'pvalue' }= 1 - gsl_cdf_chisq_P( gsl_cdf_chisq_Pinv( $data[$header{"P"}], 1 )/$lambda, 1 );
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

# Genomic Control adjustment
if (defined $gc_correction){
	print_OUT("Calculating lambda for genomic control correction");
	my $gc_lambda = get_lambda_genomic_control($assoc_chis);
	print_OUT("   '-> lambda (median) of [ $gc_lambda ]");
	if ($gc_lambda > 1){
		print_OUT("   '-> Applying GC correction");
		my $assoc_chis = [];
		foreach my $snp (keys %assoc_data) {
			next if ( $assoc_data{ $snp }->{ 'pvalue' } == 1);
			my $snp_chi = gsl_cdf_chisq_Pinv ( 1 - $assoc_data{ $snp }->{ 'pvalue' }, 1 );
			$snp_chi /=  $gc_lambda;
			$assoc_data{ $snp }->{ 'pvalue' }= 1 - gsl_cdf_chisq_P( $snp_chi, 1 );
			push @{ $assoc_chis }, $assoc_data{ $snp }->{ 'pvalue' };
		}
		$gc_lambda = get_lambda_genomic_control($assoc_chis);
		print_OUT("   '-> After correction the lambda is [ $gc_lambda ]");
	} else {
		print_OUT("   '-> GC correction not applied because lambda is less than 1");
	}
}


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
} elsif (defined $geno_probs){
	print "Getting SNP list from gprobs file\n";
	@bim = @{ get_genotypes_from_ox_format($gprobs, $gprobs_index) };
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
		'ensembl'   => $ensembl,	
      	'chr'       => $chr,
      	'start'     => $start,
      	'end'       => $end,
      	'gene_type' => $gene_type,
      	'snps'      => [],
      	'minp'      => -9,
      	'genotypes' => null,
      	'geno_mat_rows' => [],
		'cor' => null,
        'weights' => null,
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
my $out_fh_sample_score_mat = new IO::File if (defined $sample_score);
my $out_fh_sample_score_dat = new IO::File if (defined $sample_score);

if (defined $sample_score){
  print_OUT("   '-> Sample Scores: Gene scores printed to [ $out.sample_score.mat ]");
  $out_fh_sample_score_mat->open(">$out.sample_score.mat");
  print_OUT("   '-> Sample Scores: Samples data printed to [ $out.sample_score.dat ]");
  $out_fh_sample_score_dat->open(">$out.sample_score.dat");
  print $out_fh_sample_score_dat "IID FID PID MID SEX TRAIT\n";
  for (my $i = 0; $i < scalar @fam; $i++){
    print $out_fh_sample_score_dat  $fam[$i]->{iid}," ",
                                    $fam[$i]->{fid}," ",
                                    $fam[$i]->{pid}," ",
                                    $fam[$i]->{mid}," ",
                                    $fam[$i]->{sex}," ",
                                    $fam[$i]->{pheno},"\n";
  }
  $out_fh_sample_score_dat->close();
}

# create a variable that will store a ref to a hash with the weights
my $weights = {};
if (defined @weights_file){
    print_OUT("Starting to read SNP weigths");
    # create hash refs to store the name of the weight categories
    my $w_classes = {};
    my $w_counter = 0; # category counter, in case the file has not a header.
    # loop over the files and read the weights
    foreach my $w_f (@weights_file){
        ($weights,$w_counter,$w_classes) = read_weight_file($w_f,$w_counter,$w_classes,$weights, $w_header,\%snp_to_gene);
    }
    print_OUT("  '-> [ " . scalar (keys %{$weights}) . " ] weights read");
    # make a single weight vector for each snp
    foreach my $snp (keys %{$weights}){
        # get the snp weights sorted by category name
        my @tmp_all_w = map { $weights->{$snp}{$_}; } sort {$a cmp $b} keys %{ $w_classes };
        # make piddle
        $weights->{$snp} = [@tmp_all_w];
    }
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
	  if (defined $geno_probs){
		  my @lines = map { $bim_ids{$_} + 1 } @{$gene{$gn}->{geno_mat_rows}}; 
		  my ($p_mat,$g_mat) = extract_genotypes_for_snp_list($gene{$gn}->{geno_mat_rows},\@lines,$g_prob_threshold);
		  $gene{$gn}->{genotypes} *=$p_mat;
	  }
	  # Calculate the genotypes correlation matrix
	  $gene{$gn}->{cor} = corr_table($gene{$gn}->{genotypes});
	  
	  
	  # Calculate the weights for the gene
    if (defined @weights_file){
        $gene{$gn}->{weights} = generate_weights_for_a_gene($gene{$gn}->{geno_mat_rows},$weights);
    } else {
        my $n_snps = scalar @{$gene{$gn}->{geno_mat_rows}};
        $gene{$gn}->{weights} = ones $n_snps;
        $gene{$gn}->{weights} *= 1/$n_snps;
        $gene{$gn}->{weights} /= $gene{$gn}->{weights}->sumover;
    }
    # calculate gene p-values
    &gene_pvalue($gn) if (not defined $no_forge);
    &sample_score($gene{$gn},\%assoc_data) if (defined $sample_score);

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
            if (defined @weights_file){
            $gene{$gn}->{weights} = generate_weights_for_a_gene($gene{$gn}->{geno_mat_rows},$weights);
        } else {
            my $n_snps = scalar @{$gene{$gn}->{geno_mat_rows}};
            $gene{$gn}->{weights} = ones $n_snps;
            $gene{$gn}->{weights} *= 1/$n_snps;
            $gene{$gn}->{weights} /= $gene{$gn}->{weights}->sumover;
        }
		if (defined $geno_probs){
			my @lines = map { $bim_ids{$_} + 1 } @{$gene{$gn}->{geno_mat_rows}}; 
			my ($p_mat,$g_mat) = extract_genotypes_for_snp_list($gene{$gn}->{geno_mat_rows},\@lines,$g_prob_threshold);
			$gene{$gn}->{genotypes} *=$p_mat;
		}
		$gene{$gn}->{cor} = corr_table($gene{$gn}->{genotypes});
			
	  &gene_pvalue($gn) if (not defined $no_forge);
      &sample_score($gene{$gn},\%assoc_data) if (defined $sample_score);
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
} elsif (defined $geno_probs) { # in case not plink binary files provided and only a genotype prob file is given
	print_OUT("Reading genotype probabilities from [ $geno_probs ]");
	foreach my $gn (keys %gene){
		my $snp_list = [];
		my $lines = [];
		foreach my $mapped_snp (@{$gene{$gn}->{snps}}){
			next if (not exists $assoc_data{ $bim[$bim_ids{$mapped_snp}]->{snp_id} } );
			if (defined $v){ print_OUT("Adding SNP [  $bim[ $bim_ids{$mapped_snp} ]->{snp_id}  ] to genotypes of $gn"); }
			push @{$snp_list}, $mapped_snp;
			push @{$gene{$gn}->{geno_mat_rows}}, $mapped_snp;
			push @{$lines}, $bim_ids{$mapped_snp} + 1;
		}
		my ($p_mat,$g_mat) = extract_genotypes_for_snp_list($snp_list,$lines,$g_prob_threshold);
		$gene{$gn}->{genotypes} = $g_mat*$p_mat;
		$gene{$gn}->{cor} = corr_table($gene{$gn}->{genotypes});
		# Calculate the weights for the gene
		if (defined @weights_file){
			$gene{$gn}->{weights} = generate_weights_for_a_gene($gene{$gn}->{geno_mat_rows},$weights);
		} else {
			my $n_snps = scalar @{$gene{$gn}->{geno_mat_rows}};
			$gene{$gn}->{weights} = ones $n_snps;
			$gene{$gn}->{weights} *= 1/$n_snps;
			$gene{$gn}->{weights} /= $gene{$gn}->{weights}->sumover;
		}
		
		# calculate gene p-values
		&gene_pvalue($gn) if (not defined $no_forge);
		&sample_score($gene{$gn},\%assoc_data) if (defined $sample_score);
		
		# delete the gene's data to keep memory usage low
		delete($gene{$gn});
		$count++;
		&report_advance($count,$report,"Genes");	
		
	}
}else {
  print_OUT("WARNING: Gene p-values will be calculated with the precomputed correlation only. If correlation for some SNPs pairs are missing you may get wrong results, please check your inputs for completeness");
}
 
$out_fh_sample_score_mat->close() if (defined $sample_score);
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
sub extract_genotypes_for_snp_list{
	my $snp_list = shift;
	my $line_index = shift;
	my $g_prob_threshold = shift;
	my @geno_probs = ();
	my @geno_hard_coded = ();
	# loop over the snps mapped to the gene
	for (my $i = 0; $i < scalar @$snp_list; $i++){
		my $line = line_with_index($gprobs, $gprobs_index, $line_index->[$i]);		
		my @genos = split(/[\t+\s+]/,$line);
		# now loop over all samples for this snps
		my $sample_counter = 0;
		# counter start from 5 because the first columns are chromosome, SNP id, position, minor allele and major allele
		# counter increases by three because each sample has 3 genotype probabilities for the AA, AB and BB, with A the minor allele
		for (my $g = 5; $g < scalar @genos; $g +=3){
			my $snp_prob = pdl @genos[$g..$g+2];
			my $max_index = maximum_ind($snp_prob);
			my $value = undef;
			my $hard_coded = undef;
			if (($snp_prob->dsum == 0) or ($snp_prob->($max_index) < $g_prob_threshold)){
				$value = 0;
				$hard_coded = 0;
			} else {
				$value = $snp_prob->($max_index);
				$hard_coded = (1 + $max_index);
			}
			push @{ $geno_probs[$sample_counter] } , sclr $value;
			push @{ $geno_hard_coded[$sample_counter] }, sclr $hard_coded;
			$sample_counter++;
		}
	}
	my $coded_mat = double mpdl @geno_hard_coded;
	my $prob_mat = double mpdl @geno_probs;	
	return($prob_mat,$coded_mat);
}
sub generate_weights_for_a_gene {
    my $snps = shift;
    my $weights = shift;
    # get the weights for each snp in the gene
    # if there are no weights for an SNP all will be set to 0. Meaning
    # This SNP will get the min weight in all categories present.
    my @w_mat_rows = ();
    foreach my $s (@{$snps}) {
        if (not exists $weights->{$s}) { $weights->{$s} = []; }
        push @w_mat_rows, $weights->{$s};
    }
    # make a matrix with the weigths and get its dimensions
    my $W = pdl @w_mat_rows;
    if (scalar @w_mat_rows == 1){
        $W = ones 1;
    } else {
        $W = abs($W);
        @w_mat_rows = ();
        my @dims = $W->dims();
        # re-scale the weights to be between 0 and 1 (columns)
        for (my $i = 0; $i < $dims[0]; $i++) {
            next if ($W->($i,)->flat->sumover == 0);
            $W->($i,) = ($W->($i,) - $W->($i,)->min) / $W->($i,)->max;
        }
        # sum the weight for each snp
        $W = pdl map { $W->(,$_)->flat->sum; } 0 .. $dims[1] - 1;
        # make the weight for the gene sum 1.
        if ($W->min == 0){
            $W += $W->(which($W > 0))->min/$dims[0];
        }
        $W /= $W->sum;
    }
    return($W);
}

# usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
sub build_index {
    my $data_file  = shift;
    my $index_file = shift;
    my $offset     = 0;

    while (<$data_file>) {
        print $index_file pack("N", $offset);
        $offset = tell($data_file);
    }
}

# usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
# returns line or undef if LINE_NUMBER was out of range
sub line_with_index {
    my $data_file   = shift;
    my $index_file  = shift;
    my $line_number = shift;
    
    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file

    $size = length(pack("N", 0));
    $i_offset = $size * ($line_number - 1);
    
    seek($index_file, $i_offset, 0) or return;
    read($index_file, $entry, $size);
    $d_offset = unpack("N", $entry);
	if (not defined $d_offset){
		return('1');
	} else {
		seek($data_file, $d_offset, 0);
		return scalar(<$data_file>);
	}
}

sub get_genotypes_from_ox_format {
	my $geno_probs = shift;
	my $geno_probs_index = shift;
	my $index = 0;
	my @back = ();
	my $desired_line = 1;
	my $eof = 0;
	while () {
		my $line = line_with_index(*$geno_probs, *$geno_probs_index, $desired_line);
		last if ($line eq '1');
		my ($chr,$snp,$pos,$a1,$a2) = split(/\s+/,$line);
		push @back,{
			'snp_id' => $snp,
			'chr'    => $chr,
			'cm'     => 0,
			'pos'    => $pos,
			'a2'     => $a1,
			'a1'     => $a2,
        };
		$desired_line++;
	}
	print_OUT("[ " .  scalar @back . " ] SNPs on BED file");
	return ( \@back );
}

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
    # 10 indicates missing genotype, otherwise 0 and 1 point to allele 1 (minor) or allele 2 (mayor) in the BIM file, respectively
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
            printf OUT ("$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\tNA\tNA\t1\t1\n",$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp});
            delete($gene{$gn});
            return();
        }
    #if the user defined a genotype file then we need to get the number of SNPs in the gene.
    # if the user did not defined neither a genotype nor a file with the SNP-SNP correlations exit the program
    unless (defined $bfile or defined $spearman or defined $ped) {
        print_OUT("You MUST specify file with the SNP-SNP correlations or a genotype file to calculate them by myself\n\n");
        exit(1);
    }
   
   if (defined $v){ print_OUT($gene{$gn}->{cor});}
   if (defined $v){ print_OUT("Calculating effective number of tests: "); }
   # calculate number of effective tests by the Gao ($k) and Galwey ($Meff_galwey) method.	  
   my ($k,$Meff_galwey) = number_effective_tests(\$gene{$gn}->{cor});

   # get the log of the SNP p-value
   $gene{$gn}->{pvalues} = pdl @{ $gene{$gn}->{pvalues} };
   # calculate the chi-square statistics for the Makambi method and its p-value

    #
    #### WORK ON THE GENE WEIGHTS ####
    #

    if (defined $v){ print_OUT("Weigth = [ $gene{$gn}->{weights} ]"); }

    # if desired weigth by the 1/MAF
    if (defined $w_maf){
	my $MAF_w = [];
	for my $i (0 .. $n_snps - 1) {
	    my $tmp_maf = $gene{$gn}->{genotypes}->(,$i)->flat->sum/$gene{$gn}->{genotypes}->nelem;
	    $tmp_maf = 1 - $tmp_maf if ($tmp_maf > 0.5);
	    push @{$MAF_w}, $tmp_maf;
	}
	$MAF_w = pdl $MAF_w;
	$MAF_w = 1/$MAF_w;
	$gene{$gn}->{weights} *= $MAF_w->transpose;
    }
        
    # Correct the weights by the LD in the gene.
    # the new weight will be the weigthed mean of the gene.
    # the new weight will be the weigthed mean of the gene weights.
    # the weights for the mean are the correlation between the SNP, In that way the
    # weights reflect the correlation pattern of the SNPs 
    # weights reflect the correlation pattern of the SNPs

    my $w_matrix = $gene{$gn}->{weights}*abs($gene{$gn}->{cor}); # multiply the weights by the correaltions
    my @dims = $w_matrix->dims();
    $w_matrix = pdl map { $w_matrix->(,$_)->flat->sum/$gene{$gn}->{weights}->sum; } 0 .. $dims[1] - 1; # sum the rows divided by sum of the weights used
    if ($w_matrix->min == 0){ $w_matrix += $w_matrix->(which($w_matrix == 0))->min/$w_matrix->length; } # make sure NO weights equal 0
    $w_matrix /= $w_matrix->sum; # make sure weights sum 1
    $gene{$gn}->{weights} = $w_matrix/$w_matrix->sum;
   my ($forge_chi_stat,$forge_df) = get_makambi_chi_square_and_df($gene{$gn}->{cor},$gene{$gn}->{weights},$gene{$gn}->{pvalues});

   my $fisher_p_value =  1 - gsl_cdf_chisq_P($forge_chi_stat, $forge_df );
   
   # print out the results
   my $sidak = 1-(1-$gene{$gn}->{minp})**$k;
   if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t%0.3e\t%0.5f\t%0.5f\t%2d\t%3d || $gene{$gn}->{weights}\t$gene{$gn}->{pvalues}->log\t@{ $gene{$gn}->{geno_mat_rows} }\n",$gene{$gn}->{minp},$sidak,$fisher_p_value,$forge_chi_stat,$forge_df,$Meff_galwey,$n_snps,$k); }
   printf OUT ("$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t$n_snps\t$k\n",$gene{$gn}->{minp},$sidak,$fisher_p_value,$forge_chi_stat,$forge_df,$n_snps,$k);

}

# this subroutine calcultes a gene-score for each sample
# inputs are the gene pseudo-hash, the association p-values and the correlation matrix
sub sample_score {
    my $gene = shift; # pseudohash with gene information
    my $assoc = shift; # ref to a hash
	print join " ",$gene->{genotypes}->dims,"\n";
	print $gene->{genotypes}->(1:5,1:5);
	getc;
    # alleles have been coded as 1 : homozygote 1/1 minor allele, 2 heterozygous, 3: homozygote 2/2 major allele and 0: missing.
    
    if (not defined $gene ){ return(0); }
    my ($n_samples,$n_snps) = $gene->{genotypes}->dims;
    my $geno_mat = $gene->{genotypes}->copy;
    my $index = 0;
    # get snps info
    my @snp_info = map { $assoc->{$_}; } @{ $gene->{geno_mat_rows} } ;
   
    # store the SNP effect sizefor the minor allele
    my $snps_effect_size = [];
    foreach my $snp_info (@snp_info){
        if ($snp_info->{effect_size_measure} eq 'or'){
            push @{$snps_effect_size}, log $snp_info->{ $snp_info->{effect_size_measure} } ;            
        } else {
            push @{$snps_effect_size}, $snp_info->{ $snp_info->{effect_size_measure} };
        }
    }
    $snps_effect_size= pdl @{$snps_effect_size}; # make it a piddle
    # multiply the genotypes by the effect size and the weights
    #$geno_mat *= $snps_effect_size->transpose*$gene->{weights}->transpose; 
    
    # calculate the mean of the scores for each SNP
    my $score_means = dsumover $geno_mat/$geno_mat->getdim(0);
    my $sample_z = null;
    if (defined $ss_mean){
	my $mean_over_all_scores = $score_means->davg; # the expected value is the mean of the means
	my $cov = covariance($geno_mat->transpose); # calculate the covariance matrix
	my $var_covMat = $cov->flat->dsum; # this is the variance of the test
	# standarize each sample score by the mean and variance of the distribution
	$sample_z = dsumover $geno_mat->xchg(0,1)/$geno_mat->getdim(1);
	$sample_z -= $mean_over_all_scores;
	$sample_z /= $var_covMat;
    } else {
	my $sum_over_all_scores = $score_means->dsum; # the expected value is the sum of the means
	my $cov = covariance($geno_mat->transpose); # calculate the covariance matrix
	my $var_covMat = sqrt( $cov->flat->dsum ) ; # this is the variance of the test
	$sample_z = dsumover $geno_mat->xchg(0,1);
	$sample_z -= $sum_over_all_scores;
	$sample_z /= $var_covMat;
    }
	
    # add values to output line
    my $out_line = "$gene->{ensembl} $gene->{hugo} " . join " " , $sample_z->list;
    print $out_fh_sample_score_mat "$out_line\n";
}

sub covariance {
    # calculate the covariance matrix by the ML method
    my ( $X ) = shift;
    my $Diff = $X - daverage( $X->xchg(0,1) );
    my $Sigma = ( 1 / ( $X->getdim(1) - 1 ) ) * $Diff->transpose x $Diff;
    return $Sigma;
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
  # I have change this calculation following the results of the paper of Kost et al. this should improve the approximation of the test statistics
  # Kost, J. T. & McDermott, M. P. Combining dependent p-values. Statistics & Probability Letters 2002; 60: 183-190.
  # my $COR_MAT = (3.25*abs($cor) + 0.75*(abs($cor)**2)); # OLD
  my $COR_MAT = (3.263*abs($cor) + 0.710*(abs($cor)**2) + 0.027*(abs($cor)**3)); # NEW 
  my $second = $COR_MAT*$w*($w->transpose); # apply the weights 
  ($second->diagonal(0,1)) .= 0; # set the diagonal to 0
  my $varMf_m = 4*sumover($w**2) + $second->flat->sumover; # calculate the variance of the test statistics
  my $df = 8/$varMf_m; # the degrees of freedom of the test statistic
#  if (defined $sample_score_self){ $df *= 2;}
   
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

sub get_lambda_genomic_control {
	my $p = shift;
	$p = double 1 - pdl $p;
	my $chi = gsl_cdf_chisq_Pinv($p,1);
	return $chi->median/0.456;
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

sub read_weight_file {
    my $file = shift; # file name
    my $w_counter = shift; # weigth counter, to be use if the file has not header
    my $class_names =shift; # header of weights 
    my $weights = shift; # weight for each snp
    my $header_present = shift;
    my $snp_2_gene = shift;
    print_OUT("  '-> Reading weigths from [ $file ]");
    open (WEIGHTS,$file) or print_OUT("I cannot open [ $file ]") and exit(1);
    my $count = 0;
    my $w_header = ();
    while(my $line = <WEIGHTS>){
        chomp($line);
        my ($snp_id,@w) = split(/[\s\t\,]/,$line);
        if ( defined $affy_to_rsid ) {
            if ($snp_id !~ m/^rs/){
                if (exists $affy_id{$snp_id}){ $snp_id = $affy_id{$snp_id};}
            }
        }
        # if it is the first line check if the user declared header in the file
        if ($count == 0){
            if (defined $header_present){
                $w_header = get_header([@w]);
                $count++;
                next;
            } else {
                map {   $w_header->{$w_counter} = $w_counter;
                        $w_counter++;
                    } @w;
            }
            $count++;
        }
        next if (not exists $snp_2_gene->{$snp_id});
        foreach my $category (sort { $w_header->{$a} <=> $w_header->{$b} } keys %{$w_header}){
            my $val = $w[$w_header->{$category}];
            $val = 0 if ($val eq "NA");
            $val = 0 if ($val eq "");
            $weights->{$snp_id}{$category} = $val;
            $class_names->{$category} = "";
        }
    }
    close(WEIGHTS);
    return($weights,$w_counter,$class_names);
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
 	
        Input Files:
 	-ped			Genotype file in PLINK PED format
 	-map			SNP info file in PLINK MAP format
 	-bfile			Files in PLINK binary format, corresping file with extension *.bed, *.bim and *.fam. 
 	-assoc, -a		SNP association file, columns header is necessary and at leat columns with SNP and P names must be present
        -snpmap, -m		Snp-to-gene mapping file
 	-affy_to_rsid		Affy id to rs id mapping file
 	
        Output Files:
 	-out, -o		Name of the output file
 	-print_cor		print out SNP-SNP correlations
 	
        Analsis modifiers:
        -gene_list, -g		Only analyse genes on the file provided
 	-genes			Provide some gene name in command line, only these genes will be analyzed
 	-all_genes		Analyze all genes in the snp-to-gene mapping file
 	-chr			Anlyze a specific chromosome  
 	-distance, -d		Max SNP-to-gene distance allowed (in kb)
	-correlation, -cor	SNP-SNP correlation file
 	-pearson_genotypes	Calculate SNP-SNP correlation with a pearson correlation for categorical variables
        -lambda			lambda value to correct SNP statistics, if genomic control is used
	-gc_correction		Determine the lamda for the genomic control correction from the input p-values
        -no_forge               Do not do a forge analysis. e.g. if only performing sample-scoring
        
        Weigthing
        -weights, -w            File with SNP weights
        -w_header               Indicate if the SNP weight file has a header.
        -weight_by_maf, -w_maf  Weight each SNP its the minor allele frequency (MAF). The actual weight is 1/MAF.
 	
        Sample Score Analysis:
        -sample_score		Generate sample level score (it requieres sample level genotypes).
        -ox_gprobs              Genotype probabilities in OXFORD format.
	-ss_mean		Calculate the sample score as the mean of the sample's risk loading, Default is the sum.
        
	
        
=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-report>

How often to report advance. Provide an integer X and the program will report adnvance after X genes are analyzed.

=item B<-verbose, -v>

useful for debugging

=item B<-ped>

Genotype file in PLINK PED format. See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml for details on data format

=item B<-map>

SNP info file in PLINK MAP format

=item B<-bfile>

Files in PLINK binary format, with extension *.bed, *.bim and *.fam. 

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

=item B<-gc_correction>

Correct the SNP pvalues by the genomic control method. It will calculate the lambda from the data itself. PLEASE MAKE SURE YOU DO NOT FILTER THE SNPS BEFORE RUNNING FORGE OR THE CORRECTION MAY BE ERRONEOUS.

=item B<-print_cor>

print out SNP-SNP correlations

=item B<-pearson_genotypes>

Calculate SNP-SNP correlation with a pearson correlation for categorical variables as described on S. Wellek, A. Ziegler, Hum Hered 67, 128 (2009).

=item B<-distance, -d>

Max SNP-to-gene distance allowed (in kb) 

=item B<-sample_score>

Generate a gene-score for each gene. This method make use of the FORGE method but uses information from risk alleles count and their odd-ratios as weights for each SNP. The calculation is repeated on each sample, generating one gene-score for each sample. Because the p-values used are the same for all samples this score carries information about the combination of risk alleles and their odd-ratios. However, it is different from using simply the count of risk alleles because it account for the LD between the SNPs. This method is still under development and the documentation will be completed later.

=item B<-ss_mean>

Calculate the sample score as the mean of the sample's risk loading, Default is the sum.
Given a genotype matrix G with SxM dimesions, S= samples and M = markers, a genotype probability matrix P also of SxM dimesions, a vector E of M effect sizes we calculate the risk dosage matrix R = G*P*E. Then we calculate the S sample-scores by summing over the columns for each row.
The sample-scores are normalized by FORGE-SS_i = ( sample-score_i - mean)/sd, where mean is the sum of M columns means and sd is the square-root of the sum over the Ms covariance matrix.

=item B<-weights, -w>

File with SNP weights. Multiple files with weights can be used by providing the -weights for each file. the file with weigths can be tab, space or comma separated.
the first columns must be the snp id and the rest the weights. It is possible to have a header in the file, in which case you must use the option -w_header as well. SNP weights are use in the gene analysis and
to calculate the sample scores. How to use different sources of information to make a simple weight for the SNP is a not trivial issue. To simplify things I: a) for each gene's SNPs I make the weight of each category sum 1,
b) then I sum these re-sacel weights for each SNP and c) finally make sure the weight within the gene sum 1.
Please note that SNPs without weights will internally get the weights set to zero. 

=item B<-w_header>

Indicate if the SNP weight file has a header.

=item B<-weight_by_maf, -w_maf>

Weight each SNP its the minor allele frequency (MAF). The actual weight is 1/MAF.

=item B<-ox_gprobs>

File with genotype probabilities in OXFORD format. See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html for details in the file format. An example is
1	rs10489629	67400370 2 4 1 0 0 0 1 0 1 0 0
Columns are chromosome, SNP id, SNP position, minor alelle (A), major alelle (B), prob for AA at sample1,  prob for AB at sample1,  prob for BB at sample1, prob for AA at sample2, etc.


=back

=head1 DESCRIPTION

TODO


=cut
