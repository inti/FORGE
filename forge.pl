#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::GSL::CDF;
use PDL::Primitive;
use PDL::NiceSlice;
use PDL::Stats::Basic;
#use PDL::LinearAlgebra; # commented until re-implement the use of simulation to calculate p-values
use PDL::Bad;
use Data::Dumper;

# Load local functions
use GWAS_IO;
use GWAS_STATS;
use CovMatrix;

our ( $help, $man, $out, $snpmap, $bfile, $assoc, $gene_list,
    @genes, $all_genes, $analysis_chr, $report, $spearman,
    $affy_to_rsid, @weights_file, $w_header, $v, $lambda,
    $print_cor, $pearson_genotypes,$distance, $sample_score,
    $ped, $map, $ox_gprobs,$sample_score_self, $w_maf,
    $ss_mean, $gc_correction,$g_prob_threshold,
	$bgl_gprobs, $flush, $include_gene_type, $exclude_gene_type, $gmt,
	$gmt_min_size,$gmt_max_size, $use_ld_as_corr,
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
   'snpmap|m=s@'      => \$snpmap,
   'report=i'  => \$report,
   'correlation|cor=s' => \$spearman,
   'affy_to_rsid=s' => \$affy_to_rsid,
   'verbose|v' => \$v,
   'lambda=f' => \$lambda,
   'gc_correction' => \$gc_correction,
   'print_cor' => \$print_cor,
   'pearson_genotypes' => \$pearson_genotypes,
	'use_ld' => \$use_ld_as_corr,
   'distance|d=i' => \$distance, 
   'sample_score' => \$sample_score,
   'weights|w=s' => \@weights_file,
   'w_header' => \$w_header,
   'ox_gprobs=s' => \$ox_gprobs,
   'bgl_gprobs=s' => \$bgl_gprobs,
   'g_prob_threshold=f' => \$g_prob_threshold,
   'weight_by_maf|w_maf' => \$w_maf,
   'ss_mean' => \$ss_mean,
   'flush=i' => \$flush,
	'gmt=s@' =>	\$gmt,
	'gmt_min_size=i' =>	\$gmt_min_size,
	'gmt_max_size=i' =>	\$gmt_max_size,
	'gene_type|type=s@' => \$include_gene_type,
	'exclude_gene_type|exclude_type=s@' => \$exclude_gene_type, 
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);
pod2usage(0) if (not defined $assoc);


my $LOG = new IO::File; 
$LOG->open(">$out.log") or print_OUT("I can not open [ $out.log ] to write to",$LOG) and exit(1);

print_OUT("Check http://github.com/inti/FORGE/wiki for updates",$LOG);
print_OUT("LOG file will be written to [ $out.log ]",$LOG);

# defines min and max size for gene-set analyses
defined $gmt_min_size or $gmt_min_size = 2;
defined $gmt_max_size or $gmt_max_size = 999_999_999;

# define distance threshold,
defined $flush or $flush = 1000;

# define distance threshold,
defined $distance or $distance = 20;
print_OUT("Max SNP-to-gene distance allowed [ $distance ] kb",$LOG);

# $report is use to specify how often report the advance when reading input files
defined $report or $report = 50_000;
# tell if analysis is restricted to a specific chromosome
defined $analysis_chr and print_OUT("Restricting analysis to chromosome [ $analysis_chr ]",$LOG);
#  defined lambda value for Genomic Control correction
defined $lambda or $lambda = 1;
# defined threshold value for genotype probabilities
defined $g_prob_threshold or $g_prob_threshold = 1.0;
# tell if user wants to print the *.correlation file 
defined $print_cor and print_OUT("Defined -print_cor: I will print the *.correlation file (it is bulky)",$LOG);
#set output file if not set already.
defined $out or $out = "gene_based_fisher_v041.OUT";
# tell if user wants to correct p-values by genomic control 
if ($lambda != 1){ print_OUT("SNP p-value will be corrected with lambda = [ $lambda ]",$LOG);}

my $geno_probs = undef;
my $geno_probs_format = undef;
if (defined $ox_gprobs) {
    $geno_probs = $ox_gprobs ;
	$geno_probs_format = 'OXFORD';
    print_OUT("Genotype probabilities in OXFORD format will be read from [ $geno_probs ]",$LOG);
} elsif (defined $bgl_gprobs) {
    $geno_probs = $bgl_gprobs ;
	$geno_probs_format = 'BEAGLE';
    print_OUT("Genotype probabilities in BEAGLE format will be read from [ $geno_probs ]",$LOG);
}

my ($gprobs, $gprobs_index);
if (defined $geno_probs){
    $gprobs = IO::File->new();
    $gprobs_index = IO::File->new();
    my $index_name = "$geno_probs.idx";
     $gprobs->open("<$geno_probs") or print_OUT("I can not open genotype probability file [ $geno_probs ]",$LOG) and exit(1);
    if (not -e "$geno_probs.idx") {
    	print_OUT("   '-> Making index for genotype probabilities in [ $index_name ] file",$LOG);
    	$gprobs_index->open("+>$index_name") or print_OUT("Can't open $index_name for read/write: $!\n",$LOG);
    	build_index(*$gprobs, *$gprobs_index);
    } else {
        print_OUT("   '-> Found [ $geno_probs.idx ] file for genotype probabilities",$LOG);
        $gprobs_index->open("<$geno_probs.idx") or print_OUT("Can't open $index_name for read/write: $!\n",$LOG) and exit(1);
	binmode($gprobs_index);
    }
} 

my $OUT = new IO::File; 
$OUT->open(">$out") or print_OUT("I can not open [ $out ] to write to",$LOG) and exit(1);

print $OUT "Ensembl_ID\tHugo_id\tgene_type\tchromosome\tstart\tend";
print $OUT "\tmin_p\tmin_p_SIDAK\tFORGE\tFORGE_chi-square\tFORGE_df";
print $OUT "\tB_fix\tVar_fix\tB_P_fix\tB_random\tVar_random\tB_P_random";
print $OUT "\tI-squared\tQ\tQ_p-value\ttau_squared";
print $OUT "\tn_effect_Galwey\tn_effect_Gao\tn_snps\n";

# i will read the gene_list and i will load data for just this genes to speed up.
if ( not defined $all_genes and not defined @genes and not defined $gene_list){
  $all_genes = 1;
  print_OUT("Note: You did not provide an option for the set of genes to be analyzed. I will analyze all genes covered by the SNP association file. Check documentation for options -genes and -gene_list otherwise",$LOG);
}

if ( not defined $all_genes ) { # in case user want to analyze all genes
  if ( not defined @genes ) { # in case user gave a list of genes in the command line
    print_OUT("Reading Gene List from [ $gene_list ]",$LOG);
    # read file with gene list and store gene names.
    open( GL, $gene_list ) or print_OUT("I can not open [ $gene_list ]",$LOG) and exit(1);
    @genes = <GL>;
    chomp(@genes);
    close(GL);
  } else { 
	print_OUT("Read Gene List command line [ @genes  ]",$LOG); 
  }
} else {
	print_OUT("Going to analyze all genes on [ @$snpmap ] file.",$LOG); 
}

# Now lets going to read the affy id to rsid mapping. This is used to keep all ids in the
# same nomenclature
my %affy_id = ();
if ( defined $affy_to_rsid ) { # if conversion file is defined
   print_OUT("Reading AFFY to rsID mapping from [ $affy_to_rsid ]");
   open( AFFY, $affy_to_rsid ) or print_OUT("I can not open [ $affy_to_rsid ]",$LOG) and exit(1);
   while (my $affy = <AFFY>){
      chomp($affy);
      my @b = split(/\t+/,$affy);
      $affy_id{$b[0]} = $b[1];
   }
   close(AFFY);
}

# Read file with genetic association results.
print_OUT("Reading association file: [ $assoc ]",$LOG);
# create hash to store SNP information
my $assoc_chis = [] if (defined $gc_correction);

my %assoc_data = ();
open( ASSOC, $assoc ) or print_OUT("I can not open [ $assoc ]",$LOG) and exit(1);
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
   exists $header{"SNP"} or print_OUT("Not p-value columns available, these are the headers i found [ " . (keys %header) . " ]",$LOG) and exit(1);
   exists $header{"P"} or print_OUT("Not p-value columns available, these are the headers i found [ " . (keys %header) . " ]",$LOG) and exit(1);
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
	my $SE = undef;
  my $BETA = undef;
  my $R2 = undef;
  my $STAT = undef;
  if (exists $header{"A2"}){ $A2 = $data[$header{"A2"}];}
  if (exists $header{"A1"}){ $A1 = $data[$header{"A1"}];}
  if (exists $header{"OR"}){ 
	  $OR = $data[$header{"OR"}];
	  if ((exists $header{"L95"}) and (exists $header{"U95"})){ 
		  $SE = $header{"U95"} - $header{"L95"};
	  }
  }
  if (exists $header{"BETA"}){ $BETA = $data[$header{"BETA"}]; }
  if (exists $header{"R2"}){ $R2 = $data[$header{"R2"}]; }
  if (exists $header{"STAT"}){ $STAT = $data[$header{"STAT"}]; }
  if ((not defined $SE) and (exists $header{"SE"})){ $SE = $data[$header{"SE"}]; }
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
        print_OUT("ERROR: No effect size measure or direction of effect provided. I can not perform Sample Score analysis",$LOG);
        print_OUT("ERROR: Please provide an odd-ratio, beta or regression coefficient value under the header OR, BETA or STAT, respectively",$LOG);
        exit(1);
    } 
  }
  $assoc_data{ $data[$header{"SNP"}] } = {
                                        'pvalue' => 1 ,
                                        'id' => $data[$header{"SNP"}],
                                        'a1' => $A1,
                                        'a2' => $A2,
                                        'or' => $OR,
										'se' => $SE,
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
	print_OUT("\nPlease check the lambda value is correct\n",$LOG);
	exit(1);
   }  
}
close(ASSOC);

if (scalar keys %assoc_data == 0){ 
	print_OUT("\nNo SNPs with genetic association to used in the analysis\n",$LOG);
	exit(1);
}
print_OUT("[ " . scalar (keys %assoc_data) . " ] SNPs with association data",$LOG);

# Genomic Control adjustment
if (defined $gc_correction){
	print_OUT("Calculating lambda for genomic control correction",$LOG);
	my $gc_lambda = get_lambda_genomic_control($assoc_chis);
	print_OUT("   '-> lambda (median) of [ $gc_lambda ]",$LOG);
	if ($gc_lambda > 1){
		print_OUT("   '-> Applying GC correction",$LOG);
		my $assoc_chis = [];
		foreach my $snp (keys %assoc_data) {
			if ( $assoc_data{ $snp }->{ 'pvalue' } == 1){
				push @{ $assoc_chis }, $assoc_data{ $snp }->{ 'pvalue' };
				next;
			}
			my $snp_chi = gsl_cdf_chisq_Pinv ( 1 - $assoc_data{ $snp }->{ 'pvalue' }, 1 );
			$snp_chi /=  $gc_lambda;
			$assoc_data{ $snp }->{ 'pvalue' }= 1 - gsl_cdf_chisq_P( $snp_chi, 1 );
			push @{ $assoc_chis }, $assoc_data{ $snp }->{ 'pvalue' };
		}
		$gc_lambda = get_lambda_genomic_control($assoc_chis);
		print_OUT("   '-> After correction the lambda is [ $gc_lambda ]",$LOG);
	} else {
		print_OUT("   '-> GC correction not applied because lambda is less than 1",$LOG);
	}
}


#read snp-to-gene mapping and store in a hash with key equal gene name and value
# an array with the snps in the gene.
my @bim = ();
my @fam = ();
my $ped_map_genotypes;
if (defined $bfile) {
	# read the bim file with snp information and fam file with sample information
	@bim = @{ read_bim("$bfile.bim",$affy_to_rsid,\%affy_id) };
	@fam = @{ read_fam("$bfile.fam") };
	print_OUT("[ " . scalar @bim .  " ] SNPs and [ " . scalar @fam .  " ] samples in genotype file",$LOG);
} elsif (defined $ped and defined $map){
	my ($fam_ref,$bim_ref);
	($fam_ref,$ped_map_genotypes,$bim_ref) = read_map_and_ped($ped,$map,$affy_to_rsid,\%affy_id);
	@fam = @$fam_ref;
	@bim = @$bim_ref;
print_OUT("[ " . scalar @bim .  " ] SNPs and [ " . scalar @fam .  " ] samples in genotype file",$LOG);
} elsif (defined $geno_probs){
	print_OUT("Getting list of genotyped SNPs from [ $geno_probs ]",$LOG);
	@bim = @{ get_snp_list_from_ox_format($gprobs, $gprobs_index) } if ($geno_probs_format eq 'OXFORD');
	@bim = @{ get_snp_list_from_bgl_format($gprobs, $gprobs_index) } if ($geno_probs_format eq 'BEAGLE');
	print_OUT("[ " . scalar @bim .  " ] SNPs in genotype file",$LOG);
}

my %bim_ids = ();
my $index = 0;
map {
  $bim_ids{$_->{snp_id}} = $index;
  $index++;
} @bim;

print_OUT("Loading SNP-2-Gene mapping");

for (my $i = 0; $i < scalar @$snpmap; $i++){
	if ($snpmap->[$i] =~ m/\#/){
		print_OUT("   '-> Found [ # ] key on [ $snpmap->[$i] ]. I will generate file names for 26 chromosomes.");
		push @{$snpmap}, @{ make_file_name_array($snpmap->[$i]) };
		splice(@$snpmap,$i,1);
	}
}

my %gene = ();
my %snp_to_gene = ();
my %ids_map = ();
foreach my $snp_gene_mapping_file (@$snpmap){
	if (not -e $snp_gene_mapping_file){
		print_OUT("   '-> File [ $snp_gene_mapping_file ] does not exist, moving on to next file",$LOG);
		next;
	}
	open( MAP, $snp_gene_mapping_file ) or print_OUT("Can not open [ $snp_gene_mapping_file ] file",$LOG) and exit(1);
	print_OUT("   '-> Reading [ $snp_gene_mapping_file ]",$LOG);
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
			'effect_size' => undef,
			'effect_size_se' => undef,
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
		  if (defined $v ){ print_OUT("Mapping [ $s ] to [ $ensembl ]",$LOG);}
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
	   if ( scalar @{ $gene{$ensembl}->{snps} } == 0 ) { 
		   delete( $gene{$ensembl} ); 
	   } else { 
		   if (defined $v ){ print_OUT("Gene $ensembl $hugo included in the analysis with [ ", scalar @{ $gene{$ensembl}->{snps} }, " ] mapped SNPs",$LOG); }
		   $ids_map{$hugo} = $ensembl;	   
	   }
	   
	}
	close(MAP);
}
print_OUT("  '->[ " . scalar (keys %gene) . " ] Genes read from SNP-2-Gene Mapping files",$LOG);
print_OUT("  '->[ " . scalar (keys %snp_to_gene) . " ] SNPs mapped to Genes and with association results will be analyzed",$LOG);


if (scalar keys %gene == 0){
        print_OUT("No genes mapped",$LOG);
        exit(1);
}

if (defined $gmt){
	my $total_gene_sets=0;
	print_OUT("Reading gene-set definitions",$LOG);
	foreach my $gene_set_file (@$gmt){
		print_OUT("  '-> Reading [ $gene_set_file ]",$LOG);
		open (GMT,$gene_set_file) or print_OUT("I can not open [ $gene_set_file ] to read from.",$LOG) and die $!;
		while (my $line = <GMT>){
			chomp($line);
			my ( $p_name, $p_desc, @p_genes ) = split( /\t+/, $line );
			my @gene_with_snp = ();
			# loop over the genes and check if the ids match with any with SNPs
			foreach my $gn (@p_genes) {
				if ( $gn =~ m/\// ) {
					$gn =~ s/\s+//g;
					my @genes = split( /\/{1,}/, $gn );
					map {
						if ( exists $ids_map{$_} ){
							push @gene_with_snp, $ids_map{$_};
						} elsif ( exists $gene{$_} ){
							push @gene_with_snp, $_;
						} else {
							next;
						}
					} @genes;
				} else {
					if ( exists $ids_map{$gn} ){
						push @gene_with_snp, $ids_map{$gn};
					} elsif ( exists $gene{$gn} ){
						push @gene_with_snp, $gn;
					} else {
						next;
					}
				}
			}
			my %tmp = ();
			map { $tmp{$_} = ""; } @gene_with_snp;
			@gene_with_snp = keys %tmp; 
						
			next if (scalar @gene_with_snp < $gmt_min_size);
			next if (scalar @gene_with_snp > $gmt_max_size);
			
			%tmp = ();
			my ($chrs,$starts,$ends) = "";
			
			foreach my $g (@gene_with_snp) { 
				$chrs .= "$gene{$g}->{chr},";
				$starts .= "$gene{$g}->{start},";
				$ends .= "$gene{$g}->{end},";
				foreach my $s (@{ $gene{$g}->{snps} }){
					$tmp{$s} = "";				
				}
			}
			my @p_snps = keys %tmp;
			next if (scalar @p_snps < 2);
			
			map { push  @{ $snp_to_gene{ $_ } }, $p_name; } @p_snps;
			
			$gene{$p_name} = {
				'chr'       => $chrs,
				'start'     => $starts,
				'end'       => $ends,
				'gene_type' => 'gene_set',
				'snps'      => [@p_snps],
				'minp'      => -9,
				'genotypes' => null,
				'geno_mat_rows' => [],
				'cor' => null,
				'weights' => null,
				'pvalues' => [],
				'effect_size' => undef,
				'effect_size_se' => undef,
				'gene_status' => 'KNOWN',
				'ensembl' => $p_name,
				'hugo' => $p_desc,
				'desc' => $p_desc,
				'name' => $p_name,
				'genes' => [@gene_with_snp],
			};
			$total_gene_sets++;
		}
	}

	print_OUT("  '-> Just read [ $total_gene_sets ] gene-sets with mapped genes with size [ $gmt_min_size ] and [ $gmt_max_size ]",$LOG);
}


if (defined $exclude_gene_type or defined $include_gene_type){
	print_OUT("Going to filter genes based on user defined options",$LOG);
	if (defined $exclude_gene_type){
		print_OUT("   '-> Will exclude gene-types [ @$exclude_gene_type ]",$LOG);
	}
	if (defined $include_gene_type){
		print_OUT("   '-> Will include gene-types [ @$include_gene_type ]",$LOG);
	}
	foreach my $gn (keys %gene){
		# apply filter by gene type
		if (defined $exclude_gene_type){
			delete ($gene{$gn}) if ( grep $_ eq $gene{$gn}->{gene_type} , @$exclude_gene_type );
		}
		if (defined $include_gene_type){
			if (exists $gene{$gn}) {
				delete ($gene{$gn}) if (not grep $_ eq $gene{$gn}->{gene_type}, @$include_gene_type );
			}
		}
	}
}

print_OUT("Will analyse [ " . scalar (keys %gene) . " ] genes",$LOG);


# start a hash to store the SNP-to-SNP correlation values
my %correlation = ();

# if provided get the SNP-to-SNP correlation values from a tab separated file with 3 cols:snp1 snp2 correlatio_value
if (defined $spearman){
   print_OUT("Reading SNP correlation from [ $spearman ]",$LOG);
   open( SPRMN, $spearman ) or print_OUT("I cannot open [ $spearman ]",$LOG) and exit(1);
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
print_OUT("Output file will be written to [ $out ]",$LOG);

# create a variable that will store a ref to a hash with the weights
my $weights = {};
if (defined @weights_file){
    print_OUT("Starting to read SNP weigths",$LOG);
    # create hash refs to store the name of the weight categories
    my $w_classes = {};
    my $w_counter = 0; # category counter, in case the file has not a header.
    # loop over the files and read the weights
    foreach my $w_f (@weights_file){
        ($weights,$w_counter,$w_classes) = read_weight_file($w_f,$w_counter,$w_classes,$weights, $w_header,\%snp_to_gene);
    }
    print_OUT("  '-> [ " . scalar (keys %{$weights}) . " ] weights read",$LOG);
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
print_OUT("Starting to Calculate gene p-values",$LOG);

# if there are more than 100 genes change the $report variable in order to report every ~ 10 % of genes.
unless (scalar keys %gene < 100){
  $report = int((scalar keys %gene)/100 + 0.5)*10;
}

my %snp_genotype_stack = ();

if (defined $bfile) {
  print_OUT("Reading genotypes from [ $bfile.bed ]",$LOG);
  # open genotype file
  my $bed = new IO::File;
  $bed->open("<$bfile.bed") or print_OUT("I can not open binary PLINK file [ $bfile ]",$LOG) and exit(1);
  binmode($bed); # set file type to binary
  # check if the file is a PLINK file in the proper format by checking the first 3 bytes
  my ($buffer,$n_bytes); 
  my $plink_bfile_signature = "";
  read $bed, $plink_bfile_signature, 3;
  if (unpack("B24",$plink_bfile_signature) ne '011011000001101100000001'){
    print_OUT("Binary file is not in SNP-major format, please check you files\n",$LOG);
    exit(1);
  } else { print_OUT("Binary file is on SNP-major format",$LOG); }
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
      
		if (defined $v){ print_OUT("Adding SNP [  $bim[ $bim_ids{$mapped_snp} ]->{snp_id}  ] to genotypes of $gn",$LOG); }
		if (exists $snp_genotype_stack{$mapped_snp}) {
		if (defined $v){
				print_OUT("   '-> SNP [ $mapped_snp] already read",$LOG);
			}
			push @{ $matrix }, $snp_genotype_stack{$mapped_snp};
		} else { 
			# because we know the index of the SNP in the genotype file we know on which byte its information starts
			my $snp_byte_start = $N_bytes_to_encode_snp*$bim_ids{$mapped_snp};
			# here i extract the actual genotypes
			my @snp_genotypes = @{ extract_binary_genotypes(scalar @fam,$N_bytes_to_encode_snp,$snp_byte_start,$bed) };
			# store the genotypes.
			# if a snp does not use the 8 bits of a byte the rest of the bits are fill with missing values
			# here i extract the number of genotypes corresponding to the number of samples
			
			my $maf = get_maf([@snp_genotypes[0..scalar @fam - 1]] ); # check the maf of the SNP
			next if ($maf == 0 or $maf ==1);  # go to next if it is monomorphic
			push @{ $matrix }, [@snp_genotypes[0..scalar @fam - 1]];
			$snp_genotype_stack{$mapped_snp} = [@snp_genotypes[0..scalar @fam - 1]];			
		}
		
		# add snp id to matrix row names
		push @{ $gene{$gn}->{geno_mat_rows} }, $mapped_snp;
		# store the p-value of the snp
		push @{ $gene{$gn}->{pvalues} }, $assoc_data{ $mapped_snp }->{pvalue}; 
		my $effect_measure = $assoc_data{ $mapped_snp }->{effect_size_measure};
		if (defined $effect_measure){
			if ($effect_measure eq 'or'){
				push @{ $gene{$gn}->{effect_size} }, log $assoc_data{ $mapped_snp }->{$effect_measure};
			} else {
				push @{ $gene{$gn}->{effect_size} }, $assoc_data{ $mapped_snp }->{$effect_measure};
			}
			if (defined $assoc_data{ $mapped_snp }->{se}){
				if ($effect_measure eq 'or'){
					#push @{ $gene{$gn}->{effect_size_se} }, abs log $assoc_data{ $bim[ $bim_ids{$mapped_snp} ]->{snp_id} }->{se};
					push @{ $gene{$gn}->{effect_size_se} }, $assoc_data{ $mapped_snp }->{se};
				} else {
					push @{ $gene{$gn}->{effect_size_se} }, $assoc_data{ $mapped_snp }->{se};
				}
			}
		}
    }
    # generate the genotype matrix as a PDL piddle
    $gene{$gn}->{genotypes} = pdl $matrix;
	  
	# Calculate the genotypes correlation matrix
	my $more_corrs = "";
	($gene{$gn}->{cor},$gene{$gn}->{cor_ld_r},$more_corrs)  = deal_with_correlations($gene{$gn},\%correlation,$use_ld_as_corr);
	%correlation = (%correlation,%{$more_corrs});
	  
	# Calculate the weights for the gene
	$gene{$gn}->{weights} = deal_with_weights(\@weights_file,$gene{$gn},$w_maf,$weights);
	# calculate gene p-values
	my $z_based_p = z_based_gene_pvalues($gene{$gn});
	if (ref($z_based_p) ne 'HASH' and $z_based_p == -9){
		$z_based_p = {
			'B_stouffer_fix' => "NA",
			'B_stouffer_random' => "NA",
			'B_fix' => "NA",
			'B_random' => "NA",
			'V_fix' => "NA",
			'V_random' => "NA",
			'Chi_fix' => "NA",
			'Chi_random' => "NA",
			'Z_P_fix' => "NA",
			'Z_P_random' => "NA",
			'Q' => "NA",
			'Q_P' => "NA",
			'I2' => "NA",
			'tau_squared' => "NA",
			'N' => scalar @{ $gene{$gn}->{geno_mat_rows} },
		};
	} 
	my $pvalue_based_p = gene_pvalue($gn);

	print $OUT join "\t",($gene{$gn}->{ensembl},$gene{$gn}->{hugo},$gene{$gn}->{gene_type},$gene{$gn}->{chr},$gene{$gn}->{start},$gene{$gn}->{end},
		$gene{$gn}->{pvalues}->min,
		$pvalue_based_p->{sidak_min_p},
		$pvalue_based_p->{fisher},
		$pvalue_based_p->{fisher_chi},
		$pvalue_based_p->{fisher_df},
		$z_based_p->{'B_fix'},
		$z_based_p->{'V_fix'},
		$z_based_p->{'Z_P_fix'},
		$z_based_p->{'B_random'},
		$z_based_p->{'V_random'},
		$z_based_p->{'Z_P_random'},
		$z_based_p->{'I2'},
		$z_based_p->{'Q'},
		$z_based_p->{'Q_P'},
		$z_based_p->{'tau_squared'},
		$pvalue_based_p->{Meff_Galwey},
		$pvalue_based_p->{Meff_gao},
		scalar @{ $gene{$gn}->{geno_mat_rows} });
	print $OUT "\n";
	# remove LD measures and genotypes from the stack that will not use again to free memory
	foreach my $snp (  @{ $gene{$gn}->{geno_mat_rows} } ){
	  if (scalar @{ $snp_to_gene{ $snp } } == 1){
		  foreach my $pair ( keys %{$correlation{$snp}}){
			  delete($correlation{$snp}{$pair});
			  delete($correlation{$pair}{$snp});
		  }
		  delete($snp_genotype_stack{ $snp });
		  delete($snp_to_gene{ $snp });
	  } else {
		  splice(@{ $snp_to_gene{ $snp } },0,1);
	  }
	}
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
	if (defined $v){ print_OUT("Adding SNP [  ${ $bim[$index_snp] }{snp_id}  ] to genotypes of $gn",$LOG); }
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
			my ($p_mat,$d_mat) = extract_genotypes_for_snp_list($gene{$gn}->{geno_mat_rows},\@lines,$g_prob_threshold,$geno_probs_format,$gprobs,$gprobs_index);
			$gene{$gn}->{genotypes} *=$p_mat;
		}
		$gene{$gn}->{cor} = corr_table($gene{$gn}->{genotypes});
			
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
  print_OUT("WARNING: Gene p-values will be calculated with the precomputed correlation only. If correlation for some SNPs pairs are missing you may get wrong results, please check your inputs for completeness",$LOG);
}



# if the user want to get the correlation values print the *.correlation file
if (defined $print_cor){
  open (COR,">$out.correlation") or print_OUT("Cannot open [ $out.correlation ] to write to",$LOG) and exit(1);
  foreach my $snp1 (keys %correlation) {
    foreach my $snp2 (keys %{$correlation{$snp1}}  ) {
      next if ($snp1 eq $snp2);
      printf COR ("$snp1 $snp2 %.3f\n",abs($correlation{ $snp1 }{ $snp2 }));
      delete($correlation{ $snp1 }{ $snp2 });
    }
  }
  close(COR);
}

print_OUT("Well Done!!",$LOG);

$LOG->close();
$OUT->close();

exit(0);

sub deal_with_weights {
	my $w_files = shift;
	my $gn = shift;
	my $w_by_maf =shift;
	my $weights = shift;
	my $back = null;
	my $N = scalar @{ $gn->{geno_mat_rows} };
	if (defined @$w_files){
        $back = generate_weights_for_a_gene($gn->{geno_mat_rows},$weights);
    } else {
        $back = ones $N;
        $back *= 1/$N;
        $back /= $back->dsum;
    }
	
	# if desired weigth by the 1/MAF
	if (defined $w_by_maf){
		my $MAF_w = get_maf_weights($gn->{genotypes});
		$back *= $MAF_w->transpose;
	}
	
	# Correct the weights by the LD in the gene.
	# the new weight will be the weigthed mean of the gene.
	# the new weight will be the weigthed mean of the gene weights.
	# the weights for the mean are the correlation between the SNP, In that way the
	# weights reflect the correlation pattern of the SNPs 
	# weights reflect the correlation pattern of the SNPs
	my $C;
	if (defined $gn->{cor_ld_r}){
		$C = $gn->{cor_ld_r};
	} else {
		$C = $gn->{cor};
	} 
	my $w_matrix = $back * abs($C); # multiply the weights by the correaltions
	my @dims = $w_matrix->dims();
	$w_matrix = pdl map { $w_matrix->(,$_)->flat->sum/$back->dsum; } 0 .. $dims[1] - 1; # sum the rows divided by sum of the weights used
	if ($w_matrix->min == 0){ $w_matrix += $w_matrix->(which($w_matrix == 0))->min/$w_matrix->length; } # make sure NO weights equal 0
	$w_matrix /= $w_matrix->sum; # make sure weights sum 1
	
	$back = $w_matrix/$w_matrix->sum;
	return($back);
	
}
	

sub deal_with_correlations {
	my $gn = shift;
	my $correlation = shift;
	my $use_ld = shift;
	my %new_corrs = ();
	my $shrunken_matrix = cov_shrink($gn->{genotypes}->transpose);
	my $pearson_cor = $shrunken_matrix->{cor};
	my $n = scalar @{ $gn->{'geno_mat_rows'} };
	my $cor_ld_r = undef; zeroes $n,$n; 
	if (defined $use_ld){
		$cor_ld_r = zeroes $n,$n; 
		for (my $i = 0; $i < $n; $i++){
			for (my $j = $i; $j < $n; $j++){
				next if ($j == $i);
				if (exists $correlation->{$gn->{'geno_mat_rows'}->[$i]}{$gn->{'geno_mat_rows'}->[$j]}){
					set $cor_ld_r, $i, $j, $correlation->{$gn->{'geno_mat_rows'}->[$i]}{$gn->{'geno_mat_rows'}->[$j]};
					set $cor_ld_r, $j, $i, $correlation->{$gn->{'geno_mat_rows'}->[$i]}{$gn->{'geno_mat_rows'}->[$j]};
				} else {
					my $ld = calculate_LD_stats([ $gn->{'genotypes'}->(,$i)->list ],[ $gn->{'genotypes'}->(,$j)->list ]);	
					set $cor_ld_r, $i, $j, $ld->{r};
					set $cor_ld_r, $j, $i, $ld->{r};
					
					$new_corrs{ $gn->{'geno_mat_rows'}->[$i] }{ $gn->{'geno_mat_rows'}->[$j] } = $ld->{r};
					$new_corrs{ $gn->{'geno_mat_rows'}->[$j] }{ $gn->{'geno_mat_rows'}->[$i] } = $ld->{r};
				}
			}
		}

	} else {
		$cor_ld_r = undef;
	}	
	return($pearson_cor,$cor_ld_r,\%new_corrs);
}


sub get_maf_weights {
	my $genotypes = shift;
	my $MAF_w = $genotypes->sumover/(2*$genotypes->getdim(0));
	my $maf_to_flip = which($MAF_w >0.5);
	$MAF_w->index($maf_to_flip) .= 1 - $MAF_w->index($maf_to_flip);
	$MAF_w /= $MAF_w->sum;
	$MAF_w = $MAF_w->flat; 
	return($MAF_w);
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

sub make_file_name_array {
	my $file = shift;
	my @back = ();
	my @body = split(/\#/,$file);
	for my $chr (1..26){
		push @back, join "$chr", @body;
	}
	return([@back]);
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

sub gene_pvalue {
    my $gn = shift;
    if (defined $v){ print_OUT("____ $gn ____",$LOG); }
    my $n_snps = scalar @{ $gene{$gn}->{geno_mat_rows}};
    # if the gene has just 1 SNP we make that SNP's p value the gene p-value under all methods
	if ($n_snps < 2){
		if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\tNA\tNA1\t1\n",$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp}); }
		return({
			'fisher' => -1,
			'fisher_df' => -1,
			'fisher_chi' => -1,
			'sidak_min_p' => -1,
			'Meff_gao' => -1,
			'Meff_Galwey' => -1,
		});
	}

   if (defined $v){ print_OUT("Calculating effective number of tests: ",$LOG); }
   # calculate number of effective tests by the Gao ($k) and Galwey ($Meff_galwey) method.	  
   my ($k,$Meff_galwey) = number_effective_tests(\$gene{$gn}->{cor});

   # get the log of the SNP p-value
	if (ref($gene{$gn}->{pvalues}) eq "ARRAY"){
		$gene{$gn}->{pvalues} = pdl @{ $gene{$gn}->{pvalues} };
	}
   # calculate the chi-square statistics for the Makambi method and its p-value

    #
    #### WORK ON THE GENE WEIGHTS ####
    #

    if (defined $v){ print_OUT("Weigth = [ $gene{$gn}->{weights} ]",$LOG); }


   my ($forge_chi_stat,$forge_df) = get_makambi_chi_square_and_df($gene{$gn}->{cor},$gene{$gn}->{weights},$gene{$gn}->{pvalues});
   my $fisher_p_value = sclr double  1 - gsl_cdf_chisq_P($forge_chi_stat, $forge_df );
   
	
	my $sidak = sclr double 1-(1- $gene{$gn}->{pvalues}->min)**$k;
	$forge_df = sclr $forge_df;
	$forge_chi_stat = sclr $forge_chi_stat;
	my $back = {
		'fisher' => $fisher_p_value,
		'fisher_df' => $forge_df,
		'fisher_chi' => $forge_chi_stat,
		'sidak_min_p' => $sidak,
		'Meff_gao' => $k,
		'Meff_Galwey' => $Meff_galwey,
	};
   # print out the results
  
   if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t%0.3e\t%0.5f\t%0.5f\t%2d\t%3d || $gene{$gn}->{weights}\t$gene{$gn}->{pvalues}->log\t@{ $gene{$gn}->{geno_mat_rows} }\n",$gene{$gn}->{minp},$sidak,$fisher_p_value,$forge_chi_stat,$forge_df,$Meff_galwey,$n_snps,$k); }
	return($back);
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



# this subroutine print the advance of of something
# it receives 3 parameter: index, rep and tag. index is the advance so far, report how often to report and tag is name for the printing
# it will print a report when index if divisible by rep
sub report_advance {
	  my ($index,$rep,$tag) = @_;
   if (( $index/$rep - int($index/$rep)) == 0) {print_OUT(" '->Done with [ $index ] $tag",$LOG); }
}


sub print_OUT {
	my $string = shift;
	my @file_handles = @_; 	
	print scalar localtime(), "\t$string\n";
	unless (scalar @file_handles == 0){
		foreach my $fh (@file_handles){
			print $fh scalar localtime(), "\t$string\n";
		}
	}
}

sub read_weight_file {
    my $file = shift; # file name
    my $w_counter = shift; # weigth counter, to be use if the file has not header
    my $class_names =shift; # header of weights 
    my $weights = shift; # weight for each snp
    my $header_present = shift;
    my $snp_2_gene = shift;
    print_OUT("  '-> Reading weigths from [ $file ]",$LOG);
    open (WEIGHTS,$file) or print_OUT("I cannot open [ $file ]",$LOG) and exit(1);
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
	-gmt			Gene-set definition file

	Output Files:
	-out, -o		Name of the output file
	-print_cor		print out SNP-SNP correlations

	Analsis modifiers:
	-gene_list, -g		Only analyse genes on the file provided
	-genes			Provide some gene name in command line, only these genes will be analyzed
	-all_genes		Analyze all genes in the snp-to-gene mapping file
	-chr			Anlyze a specific chromosome
	-gene_type, -type	Only analyses genes of this type, e.g. protein_coding
	-exclude_gene_type, -exclude_type	Exclude genes of this type from the analysis, e.g. pseudogenes
	-distance, -d		Max SNP-to-gene distance allowed (in kb)
	-correlation, -cor	SNP-SNP correlation file
	-pearson_genotypes	Calculate SNP-SNP correlation with a pearson correlation for categorical variables
	-use_ld				Use Linkage Disequilibrium as measure of SNP-SNP correlation
	-lambda			lambda value to correct SNP statistics, if genomic control is used
	-gc_correction		Determine the lamda for the genomic control correction from the input p-values
	-gmt_min_size		Min number of genes in gene-sets to be analysed. default = 2
	-gmt_max_size		Max number of genes in gene-sets to be analysed. default = 999999999

	Weigthing
	-weights, -w            File with SNP weights
	-w_header               Indicate if the SNP weight file has a header.
	-weight_by_maf, -w_maf  Weight each SNP its the minor allele frequency (MAF). The actual weight is 1/MAF.
	
        
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
 
=item B<-assoc, -a>
 
SNP association file, columns header is necessary and at leat columns with SNP and P names must be present. 
 
=item B<-snpmap, -m>

SNP-to-gene mapping file. Format is tab separated with fields:chromosome,start,end,Ensembl id, hugo id, description, SNP1, SNP2, ..., SNPN. 
Each SNP has 4 fields separated by colons: id,position,alleles and strand. This may look overcomplicated but allows to map any kind of variation and store its information without having to use additional files. An example look like:

1       67632083        67725662        ENSG00000162594 IL23R   KNOWN   protein_coding  Interleukin-23 receptor Precursor (IL-23R) [Source:UniProtKB/Swiss-Prot;Acc:Q5VWK5]             rs7538938:67132262:T/C:1        rs72669476:67132582:C/T:1       rs72019237:67132863:C/-:1       rs61197134:67132871:C/-:1       rs11208941:67133111:T/C:1 

=item B<-affy_to_rsid>

Affy id to rs id mapping file. Tab separated file, looks like:
SNP_A-8389091	rs7593668

=item B<-gmt>
 
Gene-set definition file. Format is (tab separated) 
NAME	Description	GENE1	GENE2	GENEN
 
=item B<-out, -o>

Name of the output file. Output file will be tab separated with the following columns: Ensembl_ID,Hugo id,gene_type,chromosome,start,end,min p-value in the gene, Sidak corrected minimum p-value, FORGE p-value, FORGE chi-square, FORGE degrees of freedom, Eigen value ratio method gene p-value, EVR chi-square, EVR degrees freedom, number of snps in gene, number effective tests(Gao et al;PDMID:19434714)
An example would look like:ENSG00000162594	IL23R	protein_coding	1	67632083	67725662	3.491e-02	4.132e-01	2.128e-01	8.42305	6.04941	1.119e-01	28.17326	10.116917	 15

=item B<-print_cor>
 
print out SNP-SNP correlations
 
=item B<-gene_list, -g>

Only analyse genes on the file provided. Ids must match either Ensembl Ids or Hugo ids present in the SNP-to-gene mapping file

=item B<-genes>

Provide some gene name in command line, only these genes will be analyzed. Use like -genes IL23R

=item B<-all_genes>

Analyze all genes in the snp-to-gene mapping file

=item B<-chr>

Anlyze a specific chromosome, Chromosome code must match that of the SNP-to-gene mapping file

=item B<-gene_type, -type>

Only analyses genes of this type, e.g. protein_coding

=item B<-exclude_gene_type, -exclude_type>
 
Exclude genes of this type from the analysis, e.g. pseudogenes

=item B<-distance, -d>
 
Max SNP-to-gene distance allowed (in kb) 
 
=item B<-correlation, -cor>

SNP-SNP correlation file. Space separated file with 3 columns, first 2 the SNP ids the the 3th the correlation between them. Like:

rs6660226 rs11209018 0.089

=item B<-pearson_genotypes>
 
Calculate SNP-SNP correlation with a pearson correlation for categorical variables as described on S. Wellek, A. Ziegler, Hum Hered 67, 128 (2009).

=item B<-use_ld>
 
Use Linkage Disequilibrium as measure of SNP-SNP correlation. Default is Pearson's correlation.

=item B<-lambda>

lambda value to correct SNP statistics, if genomic control is used

=item B<-gc_correction>

Correct the SNP pvalues by the genomic control method. It will calculate the lambda from the data itself. PLEASE MAKE SURE YOU DO NOT FILTER THE SNPS BEFORE RUNNING FORGE OR THE CORRECTION MAY BE ERRONEOUS.

=item B<-gmt_min_size>
 
Min number of genes in gene-sets to be analysed. default = 2

=item B<-gmt_max_size>
 
Max number of genes in gene-sets to be analysed. default = 999999999

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
 

=head2 DESCRIPTION
 
TODO

=head3 Examples
 
=item B<1. Basic Analysis>
 
To perform a basic gene-based testing with the example files run:
 
>perl forge.pl -bfile example/example -assoc example/example.assoc -snpmap example/example.snpmap -out test

=item B<2. Adding gene-sets to the analysis>
 
To perform a basic gene-based testing with the example files run:

>perl forge.pl -bfile example/example -assoc example/example.assoc -snpmap example/example.snpmap -out test -gmt example/example.gmt
 

