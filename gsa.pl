#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
use PDL::GSL::CDF;
use PDL::Stats::Basic;
use PDL::Ufunc;
use PDL::LinearAlgebra qw(mchol);
use Data::Dumper;
use IO::File;

# Load local functions
use GWAS_IO;
use GWAS_STATS;
use CovMatrix;

our (	$help, $man, $gmt, $pval,
	$perm, $out, $max_size, 
	$min_size, $recomb_intervals,
	$report,$gene_sets,$ref_list,
	$set_stat, $input_z, $verbose_output,
	$append, $add_file_name,$all_are_background,
	$gs_coverage,$snp_assoc, $bfile,
	$snpmap, $distance, $affy_to_rsid,$gene_set_list,
	$quick_gene_cor,$gene_cor_max_dist,$print_ref_list,
	$mnd,$mnd_N,$gene_p_type
);

GetOptions(
	'help' => \$help,
	'man' => \$man,
	'gmt=s@' => \$gmt,
	'file=s' => \$pval,
	'perm=i' => \$perm,
	'out|o=s' => \$out,
	'max_size=i' => \$max_size,
	'min_size=i' => \$min_size,
	'recomb_intervals=s' => \$recomb_intervals,
	'report=i' => \$report,
	'gene_sets=s@' => \$gene_sets,
	'gene_set_list=s'   => \$gene_set_list,
	'ref_list=s' => \$ref_list,
	'set_stat=s' => \$set_stat, # stat to calculate over the sub-networks
	'z_score' => \$input_z,
	'verbose_output|v' => \$verbose_output,
	'append' => \$append,
	'add_file_name' => \$add_file_name,
	'backgroung_all_genes' => \$all_are_background,
	'gs_coverage=f' => \$gs_coverage,
	'snp_assoc=s@' => \$snp_assoc,
	'bfile=s'	=> \$bfile,
	'snpmap|m=s@' => \$snpmap,
	'distance|d=i' => \$distance,
	'affy_to_rsid=s' => \$affy_to_rsid,
	'quick_gene_cor' => \$quick_gene_cor,
	'gene_cor_max_dist=i' => \$gene_cor_max_dist,
	'print_ref_list' => \$print_ref_list,
	'mnd' => \$mnd,
	'mnd_n=i' => \$mnd_N,
	'gene_p_type' => \$gene_p_type,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);
pod2usage(0) if (not defined $pval);
pod2usage(0) if (not defined $out);
pod2usage(0) if (not defined $gmt);

my $LOG = new IO::File; 
$LOG->open(">$out.log") or print_OUT("I can not open [ $out.log ] to write to",$LOG) and exit(1);

print_OUT("Check http://github.com/inti/FORGE/wiki for updates",$LOG);
print_OUT("LOG file will be written to [ $out.log ]",$LOG);

if (defined $gs_coverage){
	print_OUT("Removing pathway with less than [ " . $gs_coverage*100 . " % ] coverage",$LOG);
}
defined $gs_coverage or $gs_coverage = 0;
if (defined $perm){
	print_OUT("Will run [ $perm ] permutations to estimate the mean and std deviation of the statistics undel the null",$LOG);
}
if (defined $input_z and not defined $perm){ $perm = 10_000; }

defined $distance or $distance = 20;
defined $report or $report = 250;
defined $max_size or $max_size = 99_999_999; 
defined $min_size or $min_size = 10;

if (defined $mnd){
	defined $mnd_N or $mnd_N = 1_000_000;
	use Pareto_Distr_Fit qw( Pgpd );
	if (defined $gene_p_type) {
		unless (grep $_ eq $gene_p_type, ('sidak','fisher','z_fix','z_random')){
			print_OUT("I do not recognise the gene p-valyes type entered [ $gene_p_type ]\n");
			print_OUT("Option are sidak, fisher, z_fix and z_random");
			print_OUT("bye!");
			exit(1);
		}
	} else { $gene_p_type = 'z_fix'; }
	print_OUT("Will use multivariate normal distribution sampling to estimate gene-gene correlations. Using [ $mnd_N ] simulations with the [ $gene_p_type] gene p-value");
}	
print_OUT("Will analyses Gene-set between [ $min_size ] and [ $max_size ] in size",$LOG);

defined $set_stat or $set_stat = 'mean';
my $network_stat = defined_set_stat($set_stat);
print_OUT("Gene-set statistics set to [ $set_stat ]",$LOG);

my %ref_genes = ();
if (defined $ref_list){
	open (REF,$ref_list) or die $!; 
	my @tmp = <REF>;
	chomp(@tmp);
	map {$ref_genes{lc($_)} = ""; } @tmp;
	print_OUT("[ " . scalar (keys %ref_genes) . " ] read from [ $ref_list ]",$LOG);
}

my %gene_data = ();
print_OUT("Reading gene stats [ $pval ]",$LOG);
open( PVAL, $pval ) or die $!;
while (my $line = <PVAL>) {
	chomp($line);
	my ( $gn, $p ) = split( /\t/, $line );
	$gn = lc($gn);
	if (defined $ref_list){
		next if (not exists $ref_genes{$gn});
	}
	next if ( $p eq "NA" );
	unless (defined  $input_z){
		next if ( $p !~ m/\d/);
		if ( $p == 0 ) {
			print_OUT("WARNING: This gene [ " . uc($gn) . " ] has p-value equal [ $p ], I do not know how to transformit to Z score",$LOG);
			next;
		}
	}
	my $z = 0;
	if (defined $input_z){
		$z = abs($p);
	} else {
		if ($p == 1){ $p = 0.99999999; }
		$z = -1 * gsl_cdf_ugaussian_Pinv($p);
	}
	$z = eval($z);
	if ( exists $gene_data{$gn} ) {
		if ( $gene_data{$gn}->{stat} > $z ) { 
			$gene_data{$gn}->{stat} = $z; 
			$gene_data{$gn}->{pvalue} = $p;
		}
		
	} else { 
		$gene_data{$gn} = {
					'stat' => $z,
					'pvalue'=> $p,
					'recomb_int' => [],
					'node_similarity'=> null,
					'genotypes' => null,
					'snps' => undef,
					'chr' => undef,
					'start' => undef,
					'end' => undef,
		};
	}
}
close(PVAL);
print_OUT("   '-> [ " . scalar (keys %gene_data) . " ] genes will be analysed",$LOG);

if (defined $gene_sets){
	print_OUT("Reading Gene-sets from command line [ " . @{$gene_sets} . " ]",$LOG);
}

if ( defined $gene_set_list ) { # in case user gave a list of genes in the command line
	print_OUT("Reading Gene-set List from [ $gene_set_list ]",$LOG);
	# read file with gene list and store gene names.
	open( GL, $gene_set_list ) or print_OUT("I can not open [ $gene_set_list ]",$LOG) and exit(1);
	my @gs = <GL>;
	chomp(@gs);
	push  @{$gene_sets},@gs; 
	close(GL);
} 


# stats for all genes in pathways
my %genes_in_paths = ();
my %gene_2_paths_map = ();
# store all pathways and it information 
my @pathways        = ();
foreach my $gene_set_file (@$gmt){
	print_OUT("Reading gene-set definitions from [ $gene_set_file ]",$LOG);
	open( GMT, $gene_set_file ) or print_OUT("Cannot open [ $gene_set_file ]") and die $!;
	while ( my $path = <GMT> ) {
		my ( $p_name, $p_desc, @p_genes ) = split( /\t/, $path );
		if (defined $gene_sets){
			next unless (grep $_ eq $p_name,@$gene_sets); 
		}
		my @p_stats = ();
		my @gene_with_p = ();
		my $total = 0;
		foreach my $gn (@p_genes) {
			$gn = lc($gn);
			if ( $gn =~ m/\// ) {
				$gn =~ s/\s+//g;
				my @genes = split( /\/{1,}/, $gn );
				map {
					$total++;
					next if ( not exists $gene_data{$_} );
					push @p_stats, $gene_data{$_}->{stat};
					push @gene_with_p, $_;
				} @genes;
			} else {
				$total++;
				next if ( not exists $gene_data{$gn} );
				push @p_stats, $gene_data{$gn}->{stat};
				push @gene_with_p, $gn;
			}
		}
		next if ( scalar @p_stats < $min_size );
		next if ( scalar @p_stats > $max_size );
		next if ( $gs_coverage > (scalar @p_stats)/$total);
		my $p = { 'name' => $p_name,
				'desc' => $p_desc ,
				'stats' => [],
				'N_all' => $total,
				'N_in' => undef,
				'genes' => \@gene_with_p,
				'gene_recomb_inter_redundant' => [],
				'intervals' => undef,
				'node_weigths' => null,
				'z_stat_raw' => undef,
				'z_stat_empi' => undef,
				'snps' => undef,
				};
		
		# this section reduces the gene sets to sets of genes of non-overlapping recombination intervals
		# it first check in which recombination intervals the genes are
		# then it will keep those that do not share recomb inter with other genes.
		# the cases where more than one gene belong to the same recombination interval
		# are solve by choosing the best p-value per recombination interval.
		$p->{N_in} = scalar @{$p->{'genes'}} if (not defined $p->{N_in});
		if (scalar @{ $p->{'stats'} } == 0){  
			@{ $p->{'stats'} } = map { $gene_data{$_}->{stat} } @{$p->{'genes'}};  
		}
		map { 
			$genes_in_paths{$_} = $gene_data{$_}->{stat}; 
			push @{ $gene_2_paths_map{$_} }, $p->{name};
		} @{ $p->{genes} };
		map { $genes_in_paths{$_} = $gene_data{$_}->{stat}; } @{ $p->{"gene_recomb_inter_redundant"} };
		push @pathways, $p;
	}
	close(GMT);
}

if (defined $ref_list){
	foreach my $gn (keys %ref_genes){
		if (exists $gene_data{lc($gn)}->{stat}){
			$genes_in_paths{lc($gn)} = $gene_data{lc($gn)}->{stat};
		}
	}
}
# if analysing only 1 pathway then reference list are all genes in pathways
if (scalar @pathways == 1){
	foreach my $gn (keys %gene_data){
		if (exists $gene_data{lc($gn)}->{stat}){
			$genes_in_paths{lc($gn)} = $gene_data{lc($gn)}->{stat};
		}
	}
}

print_OUT("   '-> [ " . scalar @pathways . " ] gene-sets will be analysed",$LOG);



# Now lets going to read the affy id to rsid mapping. This is used to keep all ids in the
# same nomenclature
my %affy_id = ();
if ( defined $affy_to_rsid ) { # if conversion file is defined
	print_OUT("Reading AFFY to rsID mapping from [ $affy_to_rsid ]",$LOG);
	open( AFFY, $affy_to_rsid ) or print_OUT("I can not open [ $affy_to_rsid ]") and exit(1);
	while (my $affy = <AFFY>){
		chomp($affy);
		my @b = split(/\t+/,$affy);
		$affy_id{$b[0]} = $b[1];
	}
	close(AFFY);
}

my %snps_covered = ();
if (defined $snp_assoc){
	print_OUT("Reading SNPs included in analysis",$LOG);
	foreach my $file (@$snp_assoc){
		open (ASSOC,$file) or print_OUT("Cannot open [ $file ]") and die $!;
		while(my $line = <ASSOC>){
			my @d = split(/[\t+\s+]/,$line);
			$snps_covered{$d[0]} = "";
		}
		close(ASSOC);
	}
	print_OUT("   '->[ " . scalar (keys %snps_covered) . " ] SNP read",$LOG);
}
my $N_bytes_to_encode_snp;
my %bim_ids = ();
my @bim = ();
my @fam = ();
my $bed = new IO::File;

if (defined $bfile){
	print_OUT("Checking genotypes on [ $bfile.bed ]",$LOG);
	# open genotype file
	$bed->open("<$bfile.bed") or print_OUT("I can not open binary PLINK file [ $bfile ]") and exit(1);
	binmode($bed); # set file type to binary
	my $plink_bfile_signature = "";
	read $bed, $plink_bfile_signature, 3;
	if (unpack("B24",$plink_bfile_signature) ne '011011000001101100000001'){
		print_OUT("Binary file is not in SNP-major format, please check you files",$LOG);
		exit(1);
	} else { 
		print_OUT("Binary file is on SNP-major format",$LOG); 
	}

	# read the bim file with snp information and fam file with sample information
	@bim = @{ read_bim("$bfile.bim",$affy_to_rsid,\%affy_id) };
	@fam = @{ read_fam("$bfile.fam") };
	print_OUT("[ " . scalar @bim .  " ] SNPs and [ " . scalar @fam .  " ] samples in genotype file",$LOG);
	my $index = 0;
	map {
		$bim_ids{$_->{snp_id}} = $index;
		$index++;
	} @bim;
	# calculate how many bytes are needed  to encode a SNP
	# each byte has 8 bits with information for 4 genotypes
	$N_bytes_to_encode_snp = (scalar @fam)/4; # four genotypes per byte
	# if not exact round it up
	if (($N_bytes_to_encode_snp - int($N_bytes_to_encode_snp)) != 0  ){ $N_bytes_to_encode_snp = int($N_bytes_to_encode_snp) + 1;}
}

if (defined $snpmap){
	print_OUT("Loading SNP-2-Gene mapping",$LOG);
	for (my $i = 0; $i < scalar @$snpmap; $i++){
		if ($snpmap->[$i] =~ m/\#/){
			print_OUT("   '-> Found [ # ] key on [ $snpmap->[$i] ]. I will generate file names for 26 chromosomes",$LOG);
			push @{$snpmap}, @{ make_file_name_array($snpmap->[$i]) };
			splice(@$snpmap,$i,1);
		}
	}
	foreach my $snp_gene_mapping_file (@$snpmap){
		if (not -e $snp_gene_mapping_file){
			print_OUT("   '-> File [ $snp_gene_mapping_file ] does not exist, moving on to next file",$LOG);
			next;
		}
		my $MAP = new IO::Handle;
		print_OUT ("   '-> Reading [ $snp_gene_mapping_file ]",$LOG);
		$MAP->open("$snp_gene_mapping_file");
		while (my $read = $MAP->getline()) {
			chomp($read);
			# the line is separate in gene info and snps. the section are separated by a tab.
			my ($chr,$start,$end,$ensembl,$hugo,$gene_status,$gene_type,$description,@m) = split(/\t+/,$read);
			#check if gene was in the list of genes i want to analyze
			my $gene_id = undef;
			if ( exists $gene_data{lc($hugo)}){
				$gene_id= lc($hugo);
			} elsif (exists $gene_data{lc($ensembl)}){
				$gene_id= lc($ensembl);
			} else {
				next;
			}		
			# exlude genes that are not in gene-sets
			next unless (exists $genes_in_paths{lc($hugo)} or exists $genes_in_paths{lc($ensembl)});

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
				next unless ( exists $snps_covered{$s} );
				next unless ( exists $bim_ids{$s} );
				$nr_snps{$s} = "";
			}
			@mapped_snps = keys %nr_snps;
			if ( scalar @mapped_snps == 0 ){
				delete( $genes_in_paths{ $gene_id } );
				next;
			}
			$gene_data{$gene_id}->{snps} = [@mapped_snps];
			$gene_data{$gene_id}->{chr} = $chr;
			$gene_data{$gene_id}->{start} = $start;
			$gene_data{$gene_id}->{end} = $end;
		}
		$MAP->close;
	}
	
	print_OUT("Reading genotypes for genes",$LOG);
}


# mean of all stats across the genome
my $background_values;
if (defined $all_are_background) {
	my @values = ();
	foreach my $gn (keys %gene_data){
		push @values, $gene_data{$gn}->{stat};
	}
	$background_values = pdl @values;
} else {
	$background_values = pdl values %genes_in_paths;
}
print_OUT("  '-> [ " . scalar $background_values->list() . " ] genes will be used to calculate the patameters of null",$LOG);
$background_values->inplace->setvaltobad( "inf" );
$background_values->inplace->setvaltobad( "-inf" );
my $mu = $background_values->davg;
# sd of all stats across the genome
my $sd = $background_values->stdv;

print_OUT("  '-> Background mean = [ $mu ] and stdv [ $sd ]",$LOG);

if (defined $print_ref_list){
	print_OUT("Printing reference list to [ $out.refList ]",$LOG);
	my @ids= ();
	if (defined $all_are_background) {
		@ids = keys %gene_data;
	} else {
		@ids = keys %genes_in_paths;
	}
	open (REFLIST,">$out.refList") or print_OUT("I cannot write to [ $out.refList ]") and exit(1);
	print REFLIST join "\n", @ids;
	close(REFLIST);
}

my %snp_genotype_stack = ();
my %correlation_stack  = ();
my %gene_genotype_var  = ();
my $i = 0;

if (defined $bfile){
	my @new_pathways = ();
	foreach my $p (@pathways) {
		report_advance($i++,$report,"Gene-Sets Genotypes",$LOG);
		next if ( scalar @{ $p->{genes} } < $min_size );
		
		my @new_genes = ();
		foreach my $gn ( @{$p->{ genes } } ){
			next if ( not exists $genes_in_paths{$gn} );
			next if ( not exists $gene_data{$gn} );
			next if ( not defined $gene_data{$gn}->{snps} );
			next if ( $gene_data{$gn}->{genotypes}->isempty == 0 );
			# this will store the genotypes
			my $matrix = [];
			my @gn_snps_with_genotypes = ();
			foreach my $snp ( @{ $gene_data{$gn}->{snps} } ){
				if (exists $snp_genotype_stack{$snp}) {
					push @{ $matrix }, $snp_genotype_stack{$snp};
					push @gn_snps_with_genotypes, $snp;
				} else { 
					# because we know the index of the SNP in the genotype file we know on which byte its information starts
					my $snp_byte_start = $N_bytes_to_encode_snp*$bim_ids{$snp};
					# here i extract the actual genotypes
					my @snp_genotypes = @{ extract_binary_genotypes(scalar @fam,$N_bytes_to_encode_snp,$snp_byte_start,$bed) };
					# store the genotypes.
					# if a snp does not use the 8 bits of a byte the rest of the bits are fill with missing values
					# here i extract the number of genotypes corresponding to the number of samples
					
					my $maf = get_maf([@snp_genotypes[0..scalar @fam - 1]] ); # check the maf of the SNP
					next if ($maf == 0 or $maf ==1);  # go to next if it is monomorphic
					push @{ $matrix }, [@snp_genotypes[0..scalar @fam - 1]];
					push @gn_snps_with_genotypes, $snp;
					$snp_genotype_stack{$snp} = [@snp_genotypes[0..scalar @fam - 1]];			
				}
			}
			next if (scalar @gn_snps_with_genotypes == 0);
			$gene_data{$gn}->{snps} = [@gn_snps_with_genotypes];
			$gene_data{$gn}->{genotypes} = pdl $matrix;
			push @new_genes, $gn;
		}
		
		next if (scalar @new_genes < $min_size);
		
		$p->{genes} = [ @new_genes ];
		my @stats = map { $gene_data{$_}->{stat} } @{ $p->{genes} }; 
		$p->{stats} = [ @stats ];
		$p->{N_in} = scalar @{ $p->{genes} };
		my ($more_corr,$more_var,$G_cor) = "";
		($p->{gene_cor_mat},$more_corr,$more_var,$G_cor) = calculate_gene_corr_mat($p->{genes},\%gene_data,\%correlation_stack,\%gene_genotype_var,$quick_gene_cor,$gene_cor_max_dist,$mnd,1000,$gene_p_type);
		$p->{ genotypes_corr } = ${ $G_cor };
		%correlation_stack = ( %correlation_stack, %{$more_corr} );
		%gene_genotype_var = ( %gene_genotype_var, %{$more_var} );
		
		# store the pathway information
		push @new_pathways, $p;
		
		foreach my $gn (  @{ $p->{genes} } ){
			if (scalar @{ $gene_2_paths_map{$gn} } == 1){
				delete($gene_data{$gn}->{genotypes});
			} else {
				splice( @{ $gene_2_paths_map{$gn} },0,1);
			}
		}
	}
	@pathways = @new_pathways;
	
	%snp_genotype_stack = ();
	%correlation_stack = ();
	%gene_genotype_var  = ();
}

print_OUT("  '->[ " . scalar (keys %gene_data) . " ] Genes read from SNP-2-Gene Mapping files with genotpe data",$LOG);
if (scalar (keys %gene_data) == 0){
	print_OUT("No gene-sets to analyze\n");
	exit(1);
}


# if permutation are performed store null distribution in here
my %null_size_dist = ();

print_OUT("Writing output to [ $out ]",$LOG);
if (defined $append){
	open (OUT,">>$out") or die $!;
} else {
	open (OUT,">$out") or die $!;
	print OUT "name\traw_p\traw_z";
	
	if (defined $bfile and defined $snpmap and defined $snp_assoc){
		print OUT "\tZ_fix\tV_fix\tZ_P_fix\tZ_random\tV_random\tZ_P_random\tI2\tQ\tQ_P\ttau_squared";
		if (defined $mnd){
			print OUT "\tSIM_Z_FIX\tSIM_Z_RANDOM\tSEEN_FIX\tSEEN_RANDOM\tN\tpareto_fix_Phat\tpareto_fix_Phatci_low\tpareto_fix_Phatci_up\tpareto_random_Phat\tpareto_random_Phatci_low\tpareto_random_Phatci_up";			
		}
	}
	if (defined $perm){
		print OUT ("\tempi_p:$set_stat\tempi_z:$set_stat\tmean_set\tmean_null\tsd_null");
	} 
	print OUT ("\tN_in\tN_all\tdesc");
	if (defined $add_file_name){
		print OUT "\tdataset";
	}
	print OUT "\n";
}

print_OUT("Starting to Analyse the [ " . scalar @pathways . " ] gene-sets",$LOG);
my $c = 0;
while (my $p = shift @pathways) {
	report_advance($c,$report,"Gene-Sets");
	$c++;
	next if ( $p->{N_in} < $min_size);
	my $m = $p->{N_in};
	my $Sm = mean( $p->{stats} );
	my $z_score = ( ( $Sm - $mu ) * (sqrt($m) ) ) / $sd;
	$p->{z_stat_raw} = $z_score;
	
	if ( $z_score == -0 ) {
		print OUT "$p->{name}\tNA\t$z_score";
	} else {
		printf OUT ("$p->{name}\t%.3e\t%.3f", 1 - gsl_cdf_ugaussian_P($z_score),$z_score);
	}

	if (defined $bfile and defined $snpmap and defined $snp_assoc){
		my $z = pdl $p->{stats};
		my $var = dsum($p->{gene_cor_mat});
		my $se = $z->nelem * ones $z->nelem;
		my $Z_STATS = get_fix_and_radom_meta_analysis($z,$se,undef,$p->{gene_cor_mat});
		my $num = 100;
		my $simulated_set_stats = simulate_mnd_gene_set($p,$Z_STATS,10,$mnd_N);
		printf OUT ("\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.2f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e", 
			$Z_STATS->{'B_fix'},
			$Z_STATS->{'V_fix'},
			$Z_STATS->{'Z_P_fix'},
			$Z_STATS->{'B_random'},
			$Z_STATS->{'V_random'},
			$Z_STATS->{'Z_P_random'},
			$Z_STATS->{'I2'},
			$Z_STATS->{'Q'},
			$Z_STATS->{'Q_P'},
			$Z_STATS->{'tau_squared'},
			$simulated_set_stats->{'z_fix'},
			$simulated_set_stats->{'z_random'},
			$simulated_set_stats->{'seen_fix'},
			$simulated_set_stats->{'seen_random'},
			$simulated_set_stats->{'N'},
			$simulated_set_stats->{'pareto_fix_Phat'},
			$simulated_set_stats->{'pareto_fix_Phatci_low'},
			$simulated_set_stats->{'pareto_fix_Phatci_up'},
			$simulated_set_stats->{'pareto_random_Phat'},
			$simulated_set_stats->{'pareto_random_Phatci_low'},
			$simulated_set_stats->{'pareto_random_Phatci_up'},
		);
		
	}
	
	if (defined $perm){
		$Sm = &$network_stat( $p->{stats} );
		if (not defined $null_size_dist{$m}){
			my $rand_stats = pdl @{ get_null_distribution($perm,[keys %genes_in_paths],$p,{}, \%gene_data)};
			$null_size_dist{$m} = { 'mean' => $rand_stats->average, 'sd' => $rand_stats->stdv };
		}
		my $z_score_empirical = ( ( $Sm - $null_size_dist{$m}->{mean}) )/$null_size_dist{$m}->{sd};
		
		$p->{z_stat_empi} = $z_score_empirical;
		
		if ( $z_score == -0 ) {
			print OUT "\tNA\tNA\t";
		} else {
			printf OUT ("\t%.3e\t%.3f\t%.3f\t%.3f\t%.3f", 1 - gsl_cdf_ugaussian_P($z_score_empirical),$z_score_empirical, $Sm,$null_size_dist{$m}->{mean},$null_size_dist{$m}->{sd});
		}
		print OUT "\t$m\t$p->{N_all}\t$p->{desc}";
		if (defined $add_file_name){
			print OUT "\t$pval";
		}
		if (defined  $verbose_output){
			map { print OUT "\t$_:",$gene_data{$_}->{stat}; } sort { $gene_data{$b}->{stat} <=> $gene_data{$a}->{stat} } @{$p->{genes}};
		}
		print OUT "\n";
	} else {
		print OUT "\t$m\t$p->{N_all}\t$p->{desc}";
		if (defined $add_file_name){
			print OUT "\t$pval";
		}
		if (defined  $verbose_output){
			map { print OUT "\t$_:",$gene_data{$_}->{stat}; } sort { $gene_data{$b}->{stat} <=> $gene_data{$a}->{stat} } @{$p->{genes}};
		}
		print OUT "\n";
	}   
}

print_OUT("Well Done",$LOG);

exit;

sub simulate_mnd_gene_set {
	my $gene_set_data = shift; # HASH ref
	my $set_stats = shift; # HASH ref
	my $target = shift;
	my $MAX = shift;
	
	my $max_step_size = 10_000;
	my $total = 0;
	my $step=100;
	
	my $se = $gene_set_data->{gene_cor_mat}->getdim(0) * ones $gene_set_data->{gene_cor_mat}->getdim(0);
	my $SEEN = zeroes 2;
	
	my $cholesky = mchol($gene_set_data->{gene_cor_mat});
	my $fix_stat = [];
	my $random_stat = [];
	while ($SEEN->min < $target){
		
		my ($sim,$c) = rmnorm($step,0,$gene_set_data->{gene_cor_mat},$cholesky);
		my $sim_chi_df1 = $sim**2;
		my $sim_p =  1 - gsl_cdf_chisq_P($sim_chi_df1,1);
		for (my $sim_n = 0; $sim_n < $sim_p->getdim(0); $sim_n++){
			my $sim_z = flat -1*gsl_cdf_ugaussian_P($sim_p->($sim_n,));
			my $sim_set_stats = get_fix_and_radom_meta_analysis($sim_z,$se,undef,$gene_set_data->{gene_cor_mat});
			
			#print $set_stats->{'Z_P_fix'}," >= ", $sim_set_stats->{'Z_P_fix'}," || ", $set_stats->{'Z_P_random'} ," >= ", $sim_set_stats->{'Z_P_random'},"\n";
			$SEEN->(0)++ if ( $set_stats->{'Z_P_fix'} >= $sim_set_stats->{'Z_P_fix'} );
			$SEEN->(1)++ if ( $set_stats->{'Z_P_random'} >= $sim_set_stats->{'Z_P_random'} );
			push @{$fix_stat}, -1*gsl_cdf_ugaussian_Pinv($sim_set_stats->{'Z_P_fix'});
			push @{$random_stat}, -1*gsl_cdf_ugaussian_Pinv($sim_set_stats->{'Z_P_random'});
		}
		$total += $step;
		
		if ($SEEN->min != 0){
			$step = 1.1*(10*($total)/$SEEN->min);
		} elsif ($step < $max_step_size){ 
			$step *=10; 
		}
		if ($step > $MAX){ $step = $MAX; }
		last if ($total > $MAX);
	}
	
	my $fix_null_stats = pdl $fix_stat;
	my $random_null_stats = pdl $random_stat;
	my $fix_observed = -1*gsl_cdf_ugaussian_Pinv($set_stats->{'Z_P_fix'});
	my $random_observed = -1*gsl_cdf_ugaussian_Pinv($set_stats->{'Z_P_random'});
	
	my ($pareto_fix_Phat,$pareto_fix_Phatci_low,$pareto_fix_Phatci_up) = Pareto_Distr_Fit::Pgpd($fix_observed,$fix_null_stats,250,0.05);
	my ($pareto_random_Phat,$pareto_random_Phatci_low,$pareto_random_Phatci_up) = Pareto_Distr_Fit::Pgpd($random_observed,$random_null_stats,250,0.05);
	
	my $back = {
		'z_fix' => sclr ($SEEN->(0)+1)/($total +1),
		'z_random' => sclr ($SEEN->(1)+1)/($total +1),
		'N' => $total,
	};
	
	$back->{'seen_fix'} = sclr $SEEN->(0);
	$back->{'seen_random'} = sclr $SEEN->(1);
	$back->{'pareto_fix_Phat'} = sclr $pareto_fix_Phat;
	$back->{'pareto_fix_Phatci_low'} = $pareto_fix_Phatci_low;
	$back->{'pareto_fix_Phatci_up'} =  $pareto_fix_Phatci_up;
	$back->{'pareto_random_Phat'} = sclr $pareto_random_Phat;
	$back->{'pareto_random_Phatci_low'} =  $pareto_random_Phatci_low;
	$back->{'pareto_random_Phatci_up'} =  $pareto_random_Phatci_up;
	
	return($back);
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

sub get_genotype_matrix_var {
	my $G = shift;
	my $C = cov_shrink($G->transpose);	
	my $V = $C->{cor}->dsum();
	return($V);
}
sub check_overlap {
	my $start_uno = shift;
	my $end_uno = shift;   
	my $start_dos = shift;   
	my $end_dos = shift;
	my $back = 0;

=h
	              start_1 .......................... end_1
	 start_2 .......................... end_2
	 
=cut
	if (($start_uno > $start_dos) and ($start_uno < $end_dos) ){$back = 1; }
	
=h
	 start_1 .......................... end_1
	                  start_2 .......................... end_2	 
=cut
	if (($start_dos > $start_uno) and ($start_dos < $end_uno) ){$back = 1; }
	
=h
	start_1 .......................... end_1
	start_2 .......................... end_2

=cut
	if (($start_uno == $start_dos) or ($end_uno == $end_dos) ){$back = 1; }
	
	# else there is no overlap
	return($back);
}


sub calculate_gene_corr_mat {
	my $genes = shift; # ARRAY ref
	my $gene_data = shift; # HASH ref
	my $gene_gene_corr = shift; # HASH reaf
	my $var = shift; # HASH ref
	my $quick_cor = shift; # 0,1
	my $max_gene_dist = shift; # integer	
	my $mnd= shift; # 0,1
	my $mnd_N = shift; # integer
	my $gene_p_type = shift; # 0,1
	
	my %new_corrs = ();
	my %new_vars = ();
	my $G_cov = ""; # will store the correlation between genotypes.
	
	# define correlation matrix
	my $corr = stretcher(ones scalar @$genes); 
	if (defined $mnd){
		my $G = null;
		my $genotypes_stack = [];
		my @mat_gene_idx = ();
		my $old_end = -1;
		for (my $i = 0; $i < scalar @$genes; $i++){
			push @{$genotypes_stack}, $gene_data->{ $genes->[$i] }->{genotypes};		
			$mat_gene_idx[ $i ]= { 
				'id' => $genes->[$i], 
				'start' => $old_end + 1, 
				'end' =>  $old_end + 1 - 1 + $gene_data->{ $genes->[$i] }->{genotypes}->getdim(1),
			};
			$old_end = $mat_gene_idx[ $i ]->{end};
		}
		$G = $G->glue(1,@$genotypes_stack);
		my $G_corr = cov_shrink($G->transpose);
		my $status = undef;
		($G_cov,$status) = check_positive_definite($G_corr->{cor},1e-8);
		if ($status == 1){
			print "Matrix never positive definite";
			exit(1);
		}
		my $cholesky = mchol($G_cov);
		my ($sim,$c) = rmnorm($mnd_N,0,$G_cov,$cholesky);
		
		my $sim_chi_df1 = $sim**2;
		my $sim_p =  1 - gsl_cdf_chisq_P($sim_chi_df1,1);
		my @null_stats = (); 
		my %effective_tests = ();
		my $stats_bag = [];
		for (my $stat_sim = 0; $stat_sim < $sim_p->getdim(0); $stat_sim++){
			my $r_p = $sim_p->($stat_sim,)->flat;
			my  @null_stats = ();
			foreach my $g (@mat_gene_idx){
				my $fake_gene = {
					'pvalues' => $r_p->($g->{start}:$g->{end}),
					'effect_size' => undef,
					'effect_size_se' => undef,
					'cor' => $G_cov->($g->{start}:$g->{end},$g->{start}:$g->{end}),
					'weights' => ((ones $r_p->($g->{start}:$g->{end})->nelem)/$r_p->($g->{start}:$g->{end})->nelem),
					'geno_mat_rows' => [ $r_p->($g->{start}:$g->{end})->list ],
				};
				my $n_stat = undef;
				if ($gene_p_type eq 'sidak' or $r_p->($g->{start}:$g->{end})->nelem == 1){
					
					if (not defined $effective_tests{ $g->{id} }){
						$effective_tests{ $g->{id} } = number_effective_tests(\$G_cov->($g->{start}:$g->{end},$g->{start}:$g->{end}));
					}
					$n_stat = 1 - (1 - $fake_gene->{'pvalues'}->min)**$effective_tests{ $g->{id} };
					$n_stat = gsl_cdf_ugaussian_Pinv($n_stat );
					
				} elsif ($gene_p_type eq 'fisher'){
					my ($sim_fisher_chi_stat,$sim_fisher_df) = get_makambi_chi_square_and_df($fake_gene->{cor},$fake_gene->{weights},$fake_gene->{'pvalues'} );
					my $sim_fisher_p_value = sclr double  1 - gsl_cdf_chisq_P($sim_fisher_chi_stat, $sim_fisher_df );
					$n_stat = gsl_cdf_ugaussian_Pinv($sim_fisher_p_value);
				} elsif (($gene_p_type eq 'z_fix') or ($gene_p_type eq 'z_random')){
					my $sim_z_gene = z_based_gene_pvalues($fake_gene,$mnd);
					$n_stat = gsl_cdf_ugaussian_Pinv($sim_z_gene->{'Z_P_fix'}) if ($gene_p_type eq 'z_fix');
					$n_stat = gsl_cdf_ugaussian_Pinv($sim_z_gene->{'Z_P_random'}) if ($gene_p_type eq 'z_random');
				}
				push @null_stats, $n_stat;
			}
			push @{$stats_bag}, [@null_stats];
		}
		$stats_bag = pdl $stats_bag;
		my $gene_stats_cor = cov_shrink($stats_bag);
		$corr = $gene_stats_cor->{cor};
	} else {
		for (my $i = 0; $i < scalar @$genes; $i++){
			# get name of gene i
			my $gn_i = $genes->[$i];
			
			# get indexes for its SNPs in the snp correlation matrix
			my $idx_i = sequence $gene_data->{ $gn_i }->{genotypes}->getdim(1);
			# get its variance if it has not been calculated already
			if (not exists $var->{ $gn_i }){
				$var->{ $gn_i } =  get_genotype_matrix_var($gene_data->{ $gn_i }->{genotypes});
				$new_vars{ $gn_i } = $var->{ $gn_i };
			}
			for (my $j = $i; $j < scalar @$genes; $j++){
				next if ($j == $i); 
				# get name of gene j
				my $gn_j = $genes->[$j];
				
				# if defined a quick gene-gene correlation the correlation will only be calculate between genes in 
				# the same chromosome
				if (not exists $gene_gene_corr->{ $gn_i }{ $gn_j } and defined $quick_cor){
					# if chromosomes are different set the correlation to 0
					if ($gene_data->{ $gn_i }->{chr} ne $gene_data->{ $gn_j }->{chr} ) {
						$gene_gene_corr->{ $gn_i }{ $gn_j } = 0;
						$new_corrs{$gn_j}{$gn_i} = 0;
					} else { # if are in the same chromosome
						# if user provided a maximum distance to evaluate correlations	
						if (defined $max_gene_dist){
							# if they overlap we will need to calculate the correlation
							# return 0 from check_overlap mean there is no overlap
							if (check_overlap($gene_data->{ $gn_i }->{start},$gene_data->{ $gn_i }->{end},check_overlap($gene_data->{ $gn_j }->{start},$gene_data->{ $gn_j }->{end}) == 0 )){
								
								# check the distance between the gene coordinates
								if ( $gene_data->{ $gn_i }->{start} > $gene_data->{ $gn_j }->{end}   ){
									#		                    start_i.....end_i
									#		start_j.....end_j
									if ( ($gene_data->{ $gn_i }->{start} - $gene_data->{ $gn_j }->{end}) > $max_gene_dist * 1000 ){
										$gene_gene_corr->{ $gn_i }{ $gn_j } = 0;
										$new_corrs{$gn_j}{$gn_i} = 0;
									}
								} elsif ( $gene_data->{ $gn_j }->{start} > $gene_data->{ $gn_i }->{end}  ){
									#		start_i.....end_i
									#							start_j.....end_j
									if ( ($gene_data->{ $gn_j }->{start} - $gene_data->{ $gn_i }->{end}) > $max_gene_dist * 1000 ){
										$gene_gene_corr->{ $gn_i }{ $gn_j } = 0;
										$new_corrs{$gn_j}{$gn_i} = 0;
									}
								}
								
							}
							
						} 
					}
				}
				
				# next if this correlation was already calculated
				if (exists $gene_gene_corr->{ $gn_i }{ $gn_j }){
					set $corr, $i ,$j, $gene_gene_corr->{ $gn_i }{ $gn_j };
					set $corr, $j ,$i, $gene_gene_corr->{ $gn_i }{ $gn_j };
					next;
				}
				# get indexes for its SNPs in the snp correlation matrix
				my $idx_j = $gene_data->{ $gn_i }->{genotypes}->getdim(1) + sequence $gene_data->{ $gn_j }->{genotypes}->getdim(1);
				# get its variance if it has not been calculated already
				if (not exists $var->{ $gn_j }){
					$var->{ $gn_j } =  get_genotype_matrix_var($gene_data->{ $gn_j }->{genotypes});
					$new_vars{ $gn_j } = $var->{ $gn_j };
				}
				
				# combine the genotype data information
				my $r_i_j = zeroes $gene_data->{ $gn_i }->{genotypes}->getdim(1), $gene_data->{ $gn_j }->{genotypes}->getdim(1);
				for ( my $i_snp = 0; $i_snp < $gene_data->{ $gn_i }->{genotypes}->getdim(1); $i_snp++ ){
					my $genotype_i = $gene_data->{ $gn_i }->{genotypes}->(,$i_snp);
					
					for ( my $j_snp = 0; $j_snp < $gene_data->{ $gn_j }->{genotypes}->getdim(1); $j_snp++ ){
						
						my $genotype_j = $gene_data->{ $gn_j }->{genotypes}->(,$j_snp);
						my $c = corr( $genotype_j, $genotype_i );
						set $r_i_j, $i_snp, $j_snp, $c->sclr;
						
					}
				}
				my $c_i_j = double $r_i_j->dsum()/sqrt( $var->{ $gn_i } * $var->{ $gn_j }  );
				set $corr, $i ,$j, $c_i_j;
				set $corr, $j ,$i, $c_i_j;
				
				$new_corrs{$gn_j}{$gn_i} = $new_corrs{$gn_i}{$gn_j} = $c_i_j;
			}
		}
	}
	return($corr,\%new_corrs,\%new_vars,\$G_cov);
}





sub get_makambi_df {
  my $cor = shift; # a pdl matrix with the varoable correlations
  my $w = shift; # a pdl vector with the weights;
  # make sure weiths sum one
  $w = $w*abs($cor); # multiply the weights by the correaltions
  my @dims = $w->dims();
  $w = pdl map { $w->(,$_)->flat->sum/$w->sum; } 0 .. $dims[1] - 1; # sum the rows divided by sum of the weights used
  if ($w->min == 0){ $w += $w->(which($w == 0))->min/$w->length; } # make sure NO weights equal 0
  $w /= $w->sum; # make sure weights sum 1
  
  # calculate the correlation matrix before the applying the weights
  # I have change this calculation following the results of the paper of Kost et al. this should improve the approximation of the test statistics
  # Kost, J. T. & McDermott, M. P. Combining dependent p-values. Statistics & Probability Letters 2002; 60: 183-190.
  # my $COR_MAT = (3.25*abs($cor) + 0.75*(abs($cor)**2));
  my $COR_MAT = (3.263*abs($cor) + 0.710*(abs($cor)**2) + 0.027*(abs($cor)**3)); 
  my $second = $COR_MAT*$w*($w->transpose); # apply the weights 
  ($second->diagonal(0,1)) .= 0; # set the diagonal to 0
  my $varMf_m = 4*sumover($w**2) + $second->flat->sumover; # calculate the variance of the test statistics
  my $df = 8/$varMf_m; # the degrees of freedom of the test statistic
  return ($df);
}


sub check_if_exist {
	my $bait = shift;
	my $array = shift;
	return( grep $_ eq $bait, @$array);
}
sub check_if_overlap {
	my $array1 = shift;
	my $array2 = shift;
	my $match = 0;

	if (scalar @$array1 < scalar @$array2){
		foreach my $g (@$array1){ $match += check_if_exist($g, $array2); }
	} else {
		foreach my $g (@$array2){ $match += check_if_exist($g, $array1); }
	}
	return($match);
}
sub std_dev {
	my $ar       = shift;
	my $elements = scalar @$ar;
	my $sum      = 0;
	my $sumsq    = 0;
	foreach (@$ar) {
		$sum   += $_;
		$sumsq += ( $_**2 );
	}
	return sqrt( $sumsq / $elements - ( ( $sum / $elements )**2 ) );
}


sub mean {
	my $ar = pdl @_;
	$ar->inplace->setvaltobad( "inf" );
	$ar->inplace->setvaltobad( "-inf" );
	return $ar->average;
}


sub get_null_distribution {
        my $N = shift; # number of statisitics to be generated
	my $genes = shift;
	my $path_info = shift;
        my $sub_net_genes = $path_info->{genes}; # genes in the subnetwork
        my $size = scalar @$sub_net_genes; # size of the subnetwork
        my $comb_int = shift; # recombination interval information
	my $gene_stats = shift;
        my @rand_stat = (); # array to store the statistics from each permutation
        for (my $i = 0; $i < $N; $i++){
		if (defined $verbose_output){
			&report_advance($i,100,"Permutations");
		}
		# if sampling conditional to recombination intervals, first check is any group of genes  
		# from the subnetwotk are in the same recombination interval. If not then sample normally.
		my %intervals = ();
		
		my @normal_sampling = ();
		my @interval_sampling = ();
		if (defined $recomb_intervals){
			if ( defined $path_info->{intervals}){
				%intervals = %{$path_info->{intervals}};
			} else {
				foreach my $g (@$sub_net_genes){
					foreach my $int (@{$gene_stats->{$g}->{recomb_int}}){
						push @{ $intervals{$int} }, $g;
					}
				}
			}
			while (my ($int,$int_genes) = each %intervals){
				if (scalar @$int_genes > 1){ push @interval_sampling,$int_genes; } # genes that are in the same recombination interval
				else { push @normal_sampling, @$int_genes;} # all other genes
			}
		} else {
			@normal_sampling = @$sub_net_genes; # if not recombination interval information has been provided sampling is normal
		}
		my @sampled_genes_normal = (); # to store genes from the sampling
		my @sampled_genes_interval = (); # to store genes from the sampling with recomb intervals
		# get genes from separated recombination intervals by sampling unconditional to the recombination intervals
		if (scalar @normal_sampling > 0){
			# if sampling is not conditional on the node degree
			my @index = @{ get_rand_index(scalar @$genes, scalar @normal_sampling)};
			@sampled_genes_normal = @$genes[@index];			
		}
		# if some genes need to be sampled conditional on recombination intervals
		
		if (scalar @interval_sampling > 0){
			# loop over all group of genes
			for (my $interval = 0; $interval < scalar @interval_sampling; $interval++){
				# get the genes ids
				my @interval_genes = @{ $interval_sampling[$interval] };
				my @seed_intervals = map { @{ $gene_stats->{ $_ }->{recomb_int} } } @interval_genes;
				my %tmp = ();
				foreach my $int (@seed_intervals) {
					map { $tmp{$_} = "" }  @{ $comb_int->{$int}->{genes} };
				}
				my @interval_seed = keys %tmp;
				my @all_intervals = keys %{ $comb_int };
				# sample one gene and then sample the rest from its interval
				my $idx = int(rand(scalar @all_intervals));
				my @interval_sampling_universe = @{ $comb_int->{$all_intervals[$idx]}->{genes} };
				# skip the interval if has less genes than need to be sampled
				if (scalar @interval_genes > @interval_sampling_universe){
					$interval--;
					next;
				}
				my $ratio = scalar @interval_seed/scalar @interval_sampling_universe;
				$ratio = 1/$ratio if ($ratio >1);
				my $rand = rand();
				if ($ratio < $rand){
					$interval--;
					next;
				}
				my @index = @{ get_rand_index(scalar @interval_sampling_universe, scalar @interval_genes)};
				my @selected_genes = sort { $gene_stats->{ $b }->{stat} <=> $gene_stats->{ $a }->{stat}} @interval_sampling_universe[@index];
				foreach my $s (@selected_genes){
					push @sampled_genes_interval, { 'id'=>$s, 'n' => scalar @interval_sampling_universe };
				}
			}
		}
		my @sampled_gene_stats = map { $gene_stats->{$_}->{stat}; } @sampled_genes_normal;
		foreach my $sampled_gene ( @sampled_genes_interval ){
			push @sampled_gene_stats, $gene_stats->{ $sampled_gene->{id} }->{stat};
		}
		my $s = &$network_stat(\@sampled_gene_stats);
		push @rand_stat, $s;
        }
        return(\@rand_stat);
}
sub defined_set_stat {
        my $stat = shift;
        if ($stat eq 'stouffer_z_score') {  
		return(\&stat_set_stouffer_z_score);
	} elsif ($stat eq 'mean'){ 
		return(\&stat_set_mean_z_score);
	} else { die("The statistics [ $stat ] you want to apply to the sub-networks is not recognized\n\n"); }

}
sub stat_set_mean_z_score {
	my $stat = shift;
	my $back = 0;
	map { $back += $_ } @$stat;
	$back /= scalar @$stat;
	return($back);
}
sub stat_set_stouffer_z_score {
	my $stat = shift;
	$stat = pdl @$stat;
	my $back = $stat->sum/sqrt(length $stat->list);
	return($back);
}

sub get_rand_index {
	my $max = shift;
	my $N = shift;
        my @universe = (0..$max-1);
        my @index = ();
	for (1 .. $N){
                my $i = int(rand(scalar @universe));
        	push @index, splice (@universe,$i,1);
	}
	return(\@index);
}


sub report_advance {
	my ($index,$rep,$tag) = @_;
	if (( $index/$rep - int($index/$rep)) == 0) {
		print_OUT("   '->Done with [ $index ] $tag",$LOG);
	}
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


__END__

=head1 NAME

 Perl implementation of PAGE: parametric analysis of gene set enrichment. As bonus includes strategies to correct for gene clusters.

=head1 DESCRIPTION

B<This program> will read a file with gene symbols and statistics and performed the PAGE gene-set analysis. It reads a at least two files: a file with gene id's and p-values and a second file with the gene-set definitions. Please check the original paper Kim SY, Volsky DJ: PAGE: parametric analysis of gene set enrichment. BMC Bioinformatics 2005, 6:144. for details of the method. In addition it implements correction for Linkage Disequilibrium which is useful when analysing results from genome-wide association studies. Please check http://github.com/inti/ for updates and documentation.  

=head1 SYNOPSIS

script [options]

	General options
 	-h, --help		print help message
 	-m, --man		print complete documentation
 	-report			how often to report advance
 	-gmt			gene-set definitions on GMT format
 	-file			gene p-value file
 	-out, -o		output file
 	-max_size		max gene-set size
 	-min_size		min gene-set size
	-ref_list		set reference list for the analysis
	
	Analysis modifiers
 	-set_stat		statistics to calculate over the sub-networks
	-z_score		input values are z_scores (the absolute values will be used)
	-cgnets			File with gene ids. Only gene-sets that contain them will be analyzed.
	-cgnets_all		Set all genes in gene-sets as CGNets seed. it basically runs the CGNets analysis for all genes.
	-Neff_gene_sets		Calculate the number effective gene-sets being analysed.
	-gs_coverage		number [0,1]. Fraction of the gene-set that must be covered by the
				experiment for the gene set to be considered in the analysis
	-recomb_intervals	correct for genes that are in the same recombination interval
	-interval_merge		Integer. Recombination intervals closer than this number will be merge. It is helpfull to
				assess the effect of residual long range LD. It can be very strict if you use a large number.
	-interval_merge_by_chr	Same as before but will define a whole chromosome as the interval. It is the extreme of the
				previous option and overrides it.
	-best_per_interval	Used by default if recomb_intervals is used. If two or more genes are in the
				same recombination interval, the one with the best statistic will be selected. This behavior
				will be mantained on the permutation sampling to calculate the null. No permutation will be run
				with this option, I may fix this in the feature.
	
	Output modifiers:
 	-append			Append results to output file rather than overwrite it
 	-add_file_name		add the input file name to the result

	Permutations:
	-perm			number of permutations
	-complete_interval_sampling	Implements a method to correct for recombination intervals but using all
					statistics in the recombination interval.

=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-report>

how often to report advance. Provide an integer X and the program will report adnvance after X networks are analyzed.

=item B<-gmt>

ggene-set definitions on GMT format

=item B<-file>

gene p-value file. Tab separated file with at least two columns: gene_id and p-value. please make sure that not p-values equal to 0 are included, those genes will be excluded from the analysis.

=item B<-perm>

number of permutations

=item B<-out>

output file: ithe output file looks like 
GO0007156	9.032e-06	4.288	4.474e-06	4.441
GO0016339	1.656e-04	3.590	2.306e-04	3.502

columns are: 
1) gene-set id
2) PAGE asymtotic p-value, z_score, 
3) p-value calculated with the null distribution calculated via sampling
4) z score calculated with a null distribution calculated via sampling.

the last 2 columsn will only be printed if permutatins are run.

=item B<-max_size>

max gene-set size

=item B<-min_size>

min gene-set size

=item B<-recomb_intervals>

correct for genes that are in the same recombination interval. On GWAS analysis one would often derive one p-value per genes and these will be correlated if the genetic variants on different genes are in Linkage Desequilibrium, as for example when genes lay on the same recombination interval (pice of genome between two recombination hot-spots). The proble is that the statistics calculated across the sub-network asssume that the gene p-values are independent. With this option is possible to provide a set of genomic interval that group genes, e.g. recombination interval. This information will be use during the MC sampling. For example, if we have a network with 10 genes, 5 of which lay on the same recombination interval. With the -recomb_intervals option set on the montecarlo sampling 5 of the genes will be obtain from a single recombination interval (elsewhere in the genome), thus assuring the genetic structure in the subnetwotk is somehow preserved. 

=item B<-interval_merge>

Integer. Recombination intervals closer than this number will be merge. It is helpfull to assess the effect of residual long range LD. It can be very strict if you use a large number.

=item B<-interval_merge_by_chr>

Same as before but will define a whole chromosome as the interval. It is the extreme of the previous option and overrides it.

=item B<-z_score>

Input values are z_score intead of p-values, for calculations the absolute value of the z-score will be taken. Permutations must be performed together with this option.

=item B<-append>

Append results to output file rather than overwrite it

=item B<-add_file_name>

Add the input file name to the result. Usefull is analysing different data sets that wil be concatenated in on single file

=item B<-ref_list>

set reference list for the analysis

=item B<-set_stat>

statistics to calculate over the sub-networks. options are stouffer_z_score and mean. mean if the default


=item B<-best_per_interval>

Used by default if recomb_intervals is used. If two or more genes are in the same recombination interval, the one with the best statistic will be selected. This behavior will be mantainedon the permutation sampling to calculate the null.

=item B<-complete_interval_sampling>

Implements a method to correct for recombination intervals but using all statistics in the recombination interval.

=item B<-cgnets>

Implements the CGNet analysis. It receives a list of genes and will restric the analysis to gene-sets that contain them.

=item B<-cgnets_all>

Set all genes in gene-sets as CGNets seed. it basically runs the CGNets analysis for all genes.

=item B<-Neff_gene_sets>

Calculate the number effective gene-sets being analysed. It can be quite slow is many pathways are under analysis but it will finish in a reasinable time. It is quite useful when performing the CGNets analysis.

=back



=cut
