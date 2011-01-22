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
use Data::Dumper;

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
	$cgnets,$complete_interval_sampling,
	$best_x_interval, $gs_coverage, $interval_merge,
	$interval_merge_by_chr, $node_similarity, $cgnets_all,
	$Neff_gene_sets, $max_processes, $snp_assoc, $bfile,
	$snpmap, $distance, $affy_to_rsid,$gene_set_list,
	$quick_gene_cor,$gene_cor_max_dist
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
	'cgnets=s' => \$cgnets,
	'best_per_interval' => \$best_x_interval,
	'complete_interval_sampling' => \$complete_interval_sampling,
	'gs_coverage=f' => \$gs_coverage,
	'interval_merge=i' => \$interval_merge,
	'interval_merge_by_chr' => \$interval_merge_by_chr,
	'node_similarity=s'=> \$node_similarity, # NOT documented
	'cgnets_all' => \$cgnets_all,
	'Neff_gene_sets' => \$Neff_gene_sets,
	'n_runs=i' => \$max_processes,
	'snp_assoc=s@' => \$snp_assoc,
	'bfile=s'	=> \$bfile,
	'snpmap|m=s@' => \$snpmap,
	'distance|d=i' => \$distance,
	'affy_to_rsid=s' => \$affy_to_rsid,
	'quick_gene_cor' => \$quick_gene_cor,
	'gene_cor_max_dist=i' =>$gene_cor_max_dist,
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
defined $interval_merge_by_chr and $interval_merge = "inf";
if (defined $interval_merge and defined $best_x_interval){
	print_OUT("Using segment distance of [ $interval_merge ] to join recombination intervals during sampling",$LOG);
}
if (defined $perm){
	print_OUT("Will run [ $perm ] permutations to estimate the mean and std deviation of the statistics undel the null",$LOG);
	if (defined $best_x_interval){
		print_OUT("Permutations with bext-per-interval option not supported at the moment.",$LOG);
		exit(0);
	}
}
if (defined $input_z and not defined $perm){ $perm = 10_000; }

defined $interval_merge or $interval_merge = 0;
defined $distance or $distance = 20;
defined $report or $report = 100;
defined $max_size or $max_size = 99_999_999; 
defined $min_size or $min_size = 10;
print_OUT("Will analyses Gene-set between [ $min_size ] and [ $max_size ] in size",$LOG);

if (defined $recomb_intervals) {
	if ((defined $best_x_interval) and (defined $complete_interval_sampling)){
		print_OUT("Please choose one of the two sampling scheems -best_per_interval or -complete_interval_samplig\nFor more information type gsa.pl -man",$LOG);
		exit(1);
	}
	if (not defined $complete_interval_sampling){
		$best_x_interval = 1;
		print_OUT("Using best gene per recombination interval for permutation sampling",$LOG);
	} elsif (defined $complete_interval_sampling) {
		print_OUT("Using all genes in recombination interval for permutation sampling",$LOG);
	}
}

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
		if ( $gene_data{$gn}->{stat} < $z ) { $gene_data{$gn}->{stat} = $z; }
		$gene_data{$gn}->{pvalue} = $p;
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

# if a CGNets analysis is performed.
# read the list of genes for which gene-sets are to be analysed
# discard the ones without gene-statistic and report how many are left

my @candidate_genes = ();
if (defined $cgnets){
	open (LIST, $cgnets) or die $!;
	print_OUT("Reading CGNets seed genes from [ $cgnets ]",$LOG);
	@candidate_genes = <LIST>;
	close(LIST);
	chomp(@candidate_genes);
	my %tmp = ();
	map { $tmp{lc($_)} = ""; } @candidate_genes;
	@candidate_genes = keys %tmp;
	print_OUT("  '-> [ " . scalar @candidate_genes . " ] unique CGNets seeds read",$LOG);
	for (my $i = 0; $i < scalar @candidate_genes; $i++){
		$candidate_genes[$i] = lc($candidate_genes[$i]);
		if (not exists $gene_data{$candidate_genes[$i]}){ splice(@candidate_genes,$i,1); }
	}
	if (scalar @candidate_genes == 0){
		print_OUT("None of the CGNets seed genes has statistics",$LOG);
		print_OUT("Bye for now",$LOG);
		exit(1);
	}
	print_OUT("  '-> of them [ " . scalar @candidate_genes . " ] have statistics",$LOG);
}



my %recomb_int = ();
if ( defined $recomb_intervals) {
	%recomb_int = %{ read_recombination_intervals($recomb_intervals,\%gene_data) };
}

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
		if (defined $best_x_interval){
			my @reduced_set = (); # store here the genes of the reduced set
			my %intervals = (); # store here the comb inter of the gene-set genes
			# loop over the gene of the gene set
			foreach my $gene ( @{ $p->{genes} } ){
				# if the gene is not map to any recombiantion interval. e.g. it is in a recomb hotspot, keep it
				if (scalar @{ $gene_data{$gene}->{recomb_int}} == 0){
					push @reduced_set, $gene;
				} else { # if the gene is in a recomb interval. store the recomb int ID and all genes map to it from the gene set
					foreach my $int (@{ $gene_data{$gene}->{ recomb_int }}){
						push @{ $intervals{$int} }, $gene;
					}
				}
			}
			# now for each of the recombination intervals with genes from the gene-set
			# check if it has more than one gene, if so choose the best p-value and keep it.
			my %nr_4_reduced_set = ();
			my %redundant_genes = ();
			if ($interval_merge > 0){
				my @all_int = sort {$a <=>$b } map { $_ =~ m/interval.(\d+)/; } keys %intervals; 
				my %int_chr = ();
				map { push @{ $int_chr{ $recomb_int{"interval.$_"}->{chr} }  }, $_; } @all_int;
				my %ri_groups = ();
				my $group_counter = 0;
				my %int_membership = ();
				foreach my $chr (keys %int_chr){
					my @ri = @{ $int_chr{$chr} };
					if (scalar @ri == 1){
						push @{ $ri_groups{$group_counter} }, $ri[0];
						$int_membership{$ri[0]} = $group_counter;
						$group_counter++;
						next;
					}
					for (my $i = 0; $i < scalar @ri; $i++){
						my $int1 = $ri[$i];
						for (my $p = $i; $p < scalar @ri; $p++){
							next if ($p == $i);
							my $int2 = $ri[$p];
							if (abs ($int1 - $int2) > $interval_merge){
								if (not exists $int_membership{$int1} and not exists $int_membership{$int2}){
									$int_membership{$int1} = $group_counter;
									push @{ $ri_groups{$group_counter} }, $int1;
									$group_counter++;
									$int_membership{$int2} = $group_counter;
									push @{ $ri_groups{$group_counter} }, $int2;
									$group_counter++;
								}
								if ( not exists $int_membership{$int1} ) {
									$int_membership{$int1} = $group_counter;
									push @{ $ri_groups{$group_counter} }, $int1;
									$group_counter++;
								}
								if ( not exists $int_membership{$int2} ) {
									$int_membership{$int2} = $group_counter;
									push @{ $ri_groups{$group_counter} }, $int2;
									$group_counter++;
								}
							} else {
								# if none has a group
								if (not exists $int_membership{$int1} and not exists $int_membership{$int2}){
									$int_membership{$int1} = $group_counter;
									$int_membership{$int2} = $group_counter;
									push @{ $ri_groups{$group_counter} }, $int1, $int2;
									$group_counter++;
									next;
								}
								if (exists $int_membership{$int1} and exists $int_membership{$int2}){
									map {
										$int_membership{$_} = $int_membership{$int2};
									} @{ $ri_groups{$int_membership{$int1}} };
									
									push @{ $ri_groups{$int_membership{$int2}} }, @{ $ri_groups{$int_membership{$int1}} };
									delete $ri_groups{$int_membership{$int1} };
								}
								if ( not exists $int_membership{$int1} ) {
									$int_membership{$int1} = $int_membership{$int2};
									push @{ $ri_groups{$int_membership{$int2}} }, $int1;
								}
								if ( not exists $int_membership{$int2} ) {
									$int_membership{$int2} = $int_membership{$int1};
									push @{ $ri_groups{$int_membership{$int1}} }, $int2;
								}
							}
						}
					}
				}
				# Merge intervals that are close to each other.
				foreach my $group (keys %ri_groups){
					next if (scalar  @{ $ri_groups{$group} } == 1);
					for (my $c = 1; $c < scalar  @{ $ri_groups{$group} }; $c++){
						push @{ $intervals{ "interval.$ri_groups{$group}->[0]" } }, @{ $intervals{ "interval.$ri_groups{$group}->[$c]" } };
						delete($intervals{ "interval.$ri_groups{$group}->[$c]" } );
					}
				}
			}
			foreach my $int (keys %intervals){
				# if only one gene in the recombination interval, keep it
				if (scalar @{ $intervals{$int} } == 1){
					map { $nr_4_reduced_set{$_} = 1; } @{ $intervals{$int} };
				} else { # if there are more than one gene, keep the one with the best p-value.
					@{ $intervals{$int} } = sort {  $gene_data{$b}->{stat} <=> $gene_data{$a}->{stat} } @{ $intervals{$int} };
					my $top = shift @{ $intervals{$int} };
					$nr_4_reduced_set{ $top } = 1 + scalar @{ $intervals{$int} }; 
					# store the ids of the genes dropped for record
					map { $redundant_genes{$_} = ""; } @{ $intervals{$int} };
				} 
			}
			# finally redefine the gene set with the independent genes and their statistics.
			@{ $p->{"gene_recomb_inter_redundant"} } = keys %redundant_genes;
			foreach my $g (keys %nr_4_reduced_set){
				push @{ $p->{'genes'} }, $g;
				if ( $nr_4_reduced_set{ $g } > 1){
					my $sidak_p = 1 - ( 1 - $gene_data{$g}->{pvalue})**$nr_4_reduced_set{ $g };
					$sidak_p -= 1e-15 if ($sidak_p ==1);
					$sidak_p =  $gene_data{$g}->{pvalue} * $nr_4_reduced_set{ $g } if ($sidak_p ==0);
					push @{ $p->{'stats'} } , -1 * gsl_cdf_ugaussian_Pinv($sidak_p);
				} else {
					push @{ $p->{'stats'} } , $gene_data{$g}->{stat};
				}
			}
			$p->{intervals} = \%intervals;
			@{$p->{'genes'}} = keys %nr_4_reduced_set;
			$p->{N_in} = scalar @{$p->{'genes'}}; #+  scalar @{$p->{"gene_recomb_inter_redundant"}};
		}
		$p->{N_in} = scalar @{$p->{'genes'}} if (not defined $p->{N_in});
		if (scalar @{ $p->{'stats'} } == 0){  @{ $p->{'stats'} } = map {$gene_data{$_}->{stat}} @{$p->{'genes'}};  }
		map { 
			$genes_in_paths{$_} = $gene_data{$_}->{stat}; 
			push @{ $gene_2_paths_map{$_} }, $p->{name};
		} @{ $p->{genes} };
		map { $genes_in_paths{$_} = $gene_data{$_}->{stat}; } @{ $p->{"gene_recomb_inter_redundant"} };
		push @pathways, $p;
	}
	close(GMT);
}
# remove redundant pathways. those that have the same genes with p-values
my %nr = ();
for (my $p = 0; $p < scalar @pathways; $p++){
	if (exists $nr{ join "",@{ $pathways[$p]->{genes} } }){ splice(@pathways,$p,1); }
	else { $nr{ join "",@{ $pathways[$p]->{genes} } } = ""; }
}

if (defined $ref_list){
	foreach my $gn (keys %ref_genes){
		if (exists $gene_data{lc($gn)}->{stat}){
			$genes_in_paths{lc($gn)} = $gene_data{lc($gn)}->{stat};
		}
	}
}
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
	my $bed = new IO::File;
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
		open( MAP, $snp_gene_mapping_file ) or print_OUT ("Can not open [ $snp_gene_mapping_file ] file") and exit(1);
		print_OUT ("   '-> Reading [ $snp_gene_mapping_file ]",$LOG);
		while ( my $read = <MAP> ) {
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
			next unless (exists $genes_in_paths{lc($hugo)} or  exists $genes_in_paths{lc($ensembl)});

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
				delete($gene_data{$gene_id});
				next;
			}
			$gene_data{$gene_id}->{snps} = [@mapped_snps];
			$gene_data{$gene_id}->{chr} = $chr;
			$gene_data{$gene_id}->{start} = $start;
			$gene_data{$gene_id}->{end} = $end;
		}
		close(MAP);
	}

	print_OUT("Reading genotypes for genes",$LOG);
}
my %snp_genotype_stack = ();
my %correlation_stack  = ();
my %gene_genotype_var  = ();
my $i = 0;

if (defined $bfile){
	foreach my $p (@pathways) {
		report_advance($i++,$report,"Gene-Sets Genotypes",$LOG);
		if ( scalar @{ $p->{genes} } < $min_size ) {
			$p->{N_in} = -1;
			next;
		}
		my @new_genes = ();
		foreach my $gn ( @{$p->{ genes } } ){
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
			$gene_data{$gn}->{genotypes} = int pdl $matrix;
			push @new_genes, $gn;
		}
		
		if (scalar @new_genes < $min_size){
			$p->{N_in} = -1;
			next;
		}
		
		$p->{genes} = [ @new_genes ];
		my @stats = map { $gene_data{$_}->{stat} } @{ $p->{genes} }; 
		$p->{stats} = [ @stats ];
		$p->{N_in} = scalar @{ $p->{genes} };
		my ($more_corr,$more_var) = "";
		($p->{gene_cor_mat},$more_corr,$more_var) = calculate_gene_corr_mat($p->{genes},\%gene_data,\%correlation_stack,\%gene_genotype_var,$quick_gene_cor,$gene_cor_max_dist);	
		%correlation_stack = ( %correlation_stack, %{$more_corr} );
		%gene_genotype_var = ( %gene_genotype_var, %{$more_var} );
		
		foreach my $gn (  @{ $p->{genes} } ){
			if (scalar @{ $gene_2_paths_map{$gn} } == 1){
				delete($gene_data{$gn}->{genotypes});
			} else {
				splice( @{ $gene_2_paths_map{$gn} },0,1);
			}
		}
	}
	%snp_genotype_stack = ();
	%correlation_stack = ();
	%gene_genotype_var  = ();
}

print_OUT("  '->[ " . scalar (keys %gene_data) . " ] Genes read from SNP-2-Gene Mapping files with genotpe data",$LOG);
if (scalar (keys %gene_data) == 0){
	print_OUT("No gene-sets to analyze\n");
	exit(1);
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
my $mu = $background_values->average;
# sd of all stats across the genome
my $sd = $background_values->stdv;

# if permutation are performed store null distribution in here
my %null_size_dist = ();

print_OUT("Writing output to [ $out ]",$LOG);
if (defined $append){
	open (OUT,">>$out") or die $!;
} else {
	open (OUT,">$out") or die $!;
	print OUT "name\traw_p\traw_z";
	
	if (defined $bfile and defined $snpmap and defined $snp_assoc){
		print OUT "\tstouffer_z_p\tstouffer_z\tstouffer_z_var";
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
		$p->{z_stouffer} = dsum($z)/sqrt($var);
		printf OUT ("\t%.3e\t%.3f\t%.3e", gsl_cdf_ugaussian_P(-1*$p->{z_stouffer}),$p->{z_stouffer},$var);
	}
	
	if (defined $perm){
		$Sm = &$network_stat( $p->{stats} );
		if (defined $recomb_intervals){
			my $rand_stats = pdl @{ get_null_distribution($perm,[keys %genes_in_paths],$p,\%recomb_int, \%gene_data) };
			$null_size_dist{$m} = { 'mean' => $rand_stats->average, 'sd' => $rand_stats->stdv };
		} else {
			if (not defined $null_size_dist{$m}){
				my $rand_stats = pdl @{ get_null_distribution($perm,[keys %genes_in_paths],$p,\%recomb_int, \%gene_data)};
				$null_size_dist{$m} = { 'mean' => $rand_stats->average, 'sd' => $rand_stats->stdv };
			}
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
# Performe CGNet analysis
if (defined $cgnets or defined $cgnets_all){
	# print out header
	open (OUTCGNETS, ">$out.cgnets") or die $!;
	print OUTCGNETS join "\t", ("gene_symbol","gene_pvalue","combined_z","N_gene_sets","Neff_gs","gene_set_name","best_p","sidak_best_p","sidak_experiment_wise_best_p","gs_name_BF","log10_best_BF","log10_global_BF","Neff_gs_exp_wise","all_gs_names","\n");
	print scalar localtime(), "\t", "Writting CGNets per gene results to [ $out.cgnets ]\n";
	# now calculate the significance for the CGNets.
	# significance will be calculated for the CGNets of each seed node.
	my @candidate_paths = ();
	# If -cgnets_all then all gene in gene-sets and all gene-sets will be analysed.
	if (defined $cgnets_all){ @candidate_genes = keys %genes_in_paths; } 
	print_OUT("   '-> [ " . scalar @candidate_genes . " ] CGNets seed genes",$LOG);
	my $index = 0;
	my %gene_index = ();
	my $p_index = 0; # this is the index of the gene-set in the @candidate_paths array
	my $g_index = 0; # this will be the index of the gene dimension in the $p_cor piddle matrix.
	my %genes_paths_map = ();
	print_OUT("Cheking overlap between CGNets seeds and gene-sets",$LOG);
	# generate and index for the gene-sets and the genes.
	# the index will be used to fill up the correlation matrix among gene-sets
	foreach my $p (@pathways){
		next if ($p->{N_in} < 3);
		# consider the genes the were previouslt eliminated due to redundancy in the recombination intervals
		my $all_p_genes = [@{$p->{genes}},@{ $p->{"gene_recomb_inter_redundant"} }];
		# check if the gene-set genes overlap with the CGNets seeds
		my $overlap = 0;
		if (defined $cgnets_all) {$overlap = 1;} # by definition all genes and gene-sets will be included so always they will overlap
		else { $overlap = check_if_overlap($all_p_genes,\@candidate_genes)}
		if ($overlap == 1){
			push @candidate_paths, $p;
			foreach my $g (@{ $all_p_genes }){
				# if this is the first time we find this gene
				# give it an index number.
				if (not exists $gene_index{ $g }){
					$gene_index{ $g } = $g_index;
					$g_index++;
				}
				# record the link between this pathway and the gene 
				push @{ $genes_paths_map{$g} }, { 'path_index' => $p_index, 'gene_index' => $gene_index{$g} };
			}
			$p_index++;
		}
	}

	print_OUT("   '-> [ " . scalar @candidate_paths . " ] gene-sets with CGNets seed genes",$LOG);
	# check if there are CGNets seed without a gene-set a report it
	my @candidates_no_genesets = ();
	map {  if (not exists $genes_paths_map{$_}){ push @candidates_no_genesets, $_; }  } @candidate_genes;
	if (scalar @candidates_no_genesets != 0){
		print_OUT("PLEASE NOTE that [ " . scalar @candidates_no_genesets . " ] CGNets seed have not been mapped to gene-sets",$LOG);
		print_OUT("   '-> Writting their ids to [ $out.cgnets.without_genesets ]",$LOG);
		open (NOMAPPED, ">$out.cgnets.without_genesets") or die $!;
		print NOMAPPED join "\n",@candidates_no_genesets;
		print NOMAPPED "\n";
		close(NOMAPPED);
		if (scalar @candidates_no_genesets == scalar @candidate_genes){
			print_OUT("Hello!! There are no CGNets seed left for analysis. Bye now.",$LOG);
			exit(1);
		}
	}
	
	my $overall_cgnet_Meff = 1;
	$overall_cgnet_Meff = scalar @candidate_paths if (not defined $Neff_gene_sets);
	my $p_cor = null;
	if (defined $Neff_gene_sets){
		# generare the matrix with gene-sets X gene statistics
		my $path_gene_content = zeroes scalar @candidate_paths, scalar keys %genes_paths_map;
		foreach my $g (keys %genes_paths_map){
			foreach my $index_pair (@{ $genes_paths_map{$g} }) {
				set $path_gene_content, $index_pair->{path_index}, $index_pair->{gene_index}, $gene_data{$g}->{stat};
			}
			
		}
		print_OUT("Calculating effective number of gene-set for experiment wise correction",$LOG);
		print_OUT("Calculating correlation matrix among gene-sets. This may take a while, go for a coffe is using -cgnets_all",$LOG);
		# this is to ensure piddle of > than 1 Gb can be generated
#		if (scalar @candidate_paths > 1_000){
#			$p_cor = float $PDL::BIGPDL=zeroes( scalar @candidate_paths, scalar @candidate_paths);
#		} else {
			$p_cor = float zeroes scalar @candidate_paths, scalar @candidate_paths;
#		}
		$p_cor->diagonal(1,0)++;
		$p_cor = $path_gene_content->transpose->corr_table;
		# calculate effective number of gene-sets for the experiment wise correction
		$overall_cgnet_Meff = number_effective_tests($p_cor);
		print_OUT("\t   '-> [ $overall_cgnet_Meff ] effective number of gene-sets",$LOG);
	} else {
		print_OUT("Will use the number of pathway as number of tests. Be aware this is a conservative estimate",$LOG);
	}
	
	my $counter  = 0;
	my $rep = int((scalar @candidate_genes)/25);
	$rep = 1 if ($rep == 0);
	# now calculate the CGNets statistics
	foreach my $seed (@candidate_genes){
		&report_advance($counter,$rep," CGNets seed");
		$counter++;
		$seed = lc($seed);
		#print "SEED $seed\n";
		#unless (not defined $node_similarity){
		#	if (not exists $row_number{ $seed }){
		#		print scalar localtime(), "\t", "Gene [ $seed ] does not have similarity measures\n";
		#		next;
		#	}
		#}
		
		# get the indexes in correlation matrix of the pathways of this CGNet seed 
		my $seed_paths_index = pdl map { $_->{path_index} } @{ $genes_paths_map{$seed}};
		# initialize the statistics overall CGNets 
		my $cgnet_z = 0;
		# number of gene-sets containing this seed
		my $N_paths = @{ $genes_paths_map{$seed} };
		if ($N_paths > 1){
			# initialize number of effective CGNets as number of gene-sets
			my $seed_Meff_paths = $N_paths;
			my @seed_paths_names = (); # names of these gene-sets
			my $ps = []; # store here the statistics of each gene-set
			my $Zs = [];
			my $Ms = [];
			my $seed_per_path_Z = [];
			foreach (@candidate_paths[$seed_paths_index->list]) {
				my $seed_stat_for_gene_set = undef;
				if (not defined $best_x_interval){
					$seed_stat_for_gene_set = $gene_data{$seed}->{stat};					
				} elsif (grep $_ eq $seed, $_->{genes}){
					$seed_stat_for_gene_set = $gene_data{$seed}->{stat};
				} else {
					my @seed_intervals = @{ $gene_data{$seed}->{recomb_int} };
					foreach my $g (@{$_->{genes}}){
						if (check_if_overlap($gene_data{$g}->{recomb_int},$gene_data{$seed}->{recomb_int})){
							$seed_stat_for_gene_set = $gene_data{$g}->{stat};
							last;
						}
					}
				}
				if (not defined $seed_stat_for_gene_set ) {
					print_OUT("I could not match the interval statistics for [ $seed ]",$LOG);
					exit;
				}
				push @{ $Zs }, $_->{z_stat_raw};
				push @{ $Ms }, $_->{N_in};
				push @{ $seed_per_path_Z }, $seed_stat_for_gene_set;
				push @seed_paths_names,$_->{name};
			}
			my $P_Gplus = double $gene_data{$seed}->{pvalue};
			my $P_Gminus = 1 - double $gene_data{$seed}->{pvalue};
			$Zs = double pdl $Zs; # the network z-scores
			$Ms = double pdl $Ms; # the network sizes
			$seed_per_path_Z = pdl @{ $seed_per_path_Z }; # the contribution of the seed node to each gene-set
			# these are the p-values for the gene-sets including the seed.
			my $P_ni_Gplus = double 1-gsl_cdf_ugaussian_P( $Zs ); # now as a piddle
			
			# Calculate the corrected statistic for each gene-set by using the stats of the seed node as if it were under the oposite hypothesis
			my $P_ni_Gminus = double ($Zs*$sd)/sqrt($Ms)  + $mu  - $seed_per_path_Z/$Ms +-$seed_per_path_Z/$Ms; 
			$P_ni_Gminus = ( ( $P_ni_Gminus - $mu) * (sqrt($Ms - 1 ) ) )/$sd; # here are the gene-set z-scores minus the seed contribution
			$P_ni_Gminus = double 1-gsl_cdf_ugaussian_P($P_ni_Gminus); # now as p-values
			# calculate the gene-sets stats without the seed
			my $P_ni_noG = double ($Zs*$sd)/sqrt($Ms)  + $mu  - $seed_per_path_Z/$Ms; 
			$P_ni_noG = ( ( $P_ni_noG - $mu) * (sqrt($Ms - 1 ) ) )/$sd; # here are the gene-set z-scores minus the seed contribution
			$P_ni_noG = double 1-gsl_cdf_ugaussian_P($P_ni_noG); # now as p-values
			# the overall BF has information for the contribution of all gene-sets
			my ($P_N_Gplus,$P_N_Gminus) = undef;
			
			if (defined $Neff_gene_sets){
				# extract the correlations from the gene-set correlation matrix
				my $path_cor = $p_cor->($seed_paths_index,$seed_paths_index);
				# calculate number of effective tests
				$seed_Meff_paths = number_effective_tests($path_cor);
				$path_cor = abs($path_cor);
				# calculate the makambi fisher p-value over the gene-sets of this seed
				my $w = ones $N_paths; # the weights
				my ($chi_stat_P_ni_Gplus,$df_stat_P_ni_Gplus) = get_makambi_chi_square_and_df($path_cor, $w, $P_ni_Gplus);
				# p-value over all networks
				$P_N_Gplus = double 1 - gsl_cdf_chisq_P($chi_stat_P_ni_Gplus,$df_stat_P_ni_Gplus);
				# change the to other tail
				# the z-score over all netwotks
				$cgnet_z = double -1 * gsl_cdf_ugaussian_Pinv( $P_N_Gplus );
				# the -log10 of the p-values over all networks
				
				my ($chi_stat_P_ni_Gminus,$df_stat_P_ni_Gminus) = get_makambi_chi_square_and_df($path_cor, $w, $P_ni_Gminus);
				my $P_N_Gminus = double 1 - gsl_cdf_chisq_P($chi_stat_P_ni_Gminus,$df_stat_P_ni_Gminus);
				
			} else {
				$cgnet_z = average pdl map { $_->{z_stat_raw}; } @candidate_paths[$seed_paths_index->list];
				$P_N_Gminus = double 1 - gsl_cdf_ugaussian_P(davg (-1 * gsl_cdf_ugaussian_Pinv($P_ni_Gminus))) ;
				$P_N_Gplus = double 1 - gsl_cdf_ugaussian_P($cgnet_z);
			}
			# now lets calculate the p(G|N) = p(G)*p(N|G)/p(N)
			# correct P(Ni|G) for number of gene-sets tested
			my $log10_P_G_Ni = double log10(1-$P_Gplus) + log10(1-$P_ni_Gplus) - log10( double((1-$P_ni_Gplus)*$P_Gplus) + double((1-$P_ni_Gminus)*$P_Gminus));
			# this are the bayes factors for each gene-set: p(G|N)/p(G)
			my $log10_BF = double $log10_P_G_Ni - log10(1-$P_Gplus);
			
			my $log10_P_G_N = double log10( $P_Gplus ) + log10( $P_N_Gplus ) - log10( double($P_N_Gplus*$P_Gplus) + double($P_N_Gminus*$P_Gminus) );
			my $log10_overall_BF = double $log10_P_G_N - log10( $P_Gplus );
			my $best_bf_path_name = $seed_paths_names[$log10_BF->maximum_ind];
			# the sidak corrected tests for the best gene-set
			my $sidak_p = 1 - (1 - $P_ni_noG->min)**$seed_Meff_paths;
			my $sidak_experiment_wise = 1 - (1 - $P_ni_noG->min)**$overall_cgnet_Meff;
			$seed = uc($seed);
			my $best_path_name = $seed_paths_names[$P_ni_noG->minimum_ind];
			my $path_names = join ",", @seed_paths_names;
			printf OUTCGNETS ("$seed\t%.3e\t%.3f\t$N_paths\t%.3f\t$best_path_name\t%.2e\t%.2e\t%.2e\t$best_bf_path_name\t%.3f\t%.3f\t%.3f\t$path_names\n",
					$P_Gplus,
					$cgnet_z,$seed_Meff_paths,
					$P_ni_noG->min,$sidak_p,
					$sidak_experiment_wise,
					$log10_BF->max,
					$log10_overall_BF,
					$overall_cgnet_Meff);		
		} elsif (scalar @{ $genes_paths_map{$seed}} == 1){
			my $seed_stat_for_gene_set = undef;
			if (not defined $best_x_interval){
				$seed_stat_for_gene_set = $gene_data{$seed}->{stat};					
			} elsif (grep $_ eq $seed, $_->{genes}){
				$seed_stat_for_gene_set = $gene_data{$seed}->{stat};
			} else {
				my @seed_intervals = @{ $gene_data{$seed}->{recomb_int} };
				foreach my $g (@{$_->{genes}}){
					if (check_if_overlap($gene_data{$g}->{recomb_int},$gene_data{$seed}->{recomb_int})){
						$seed_stat_for_gene_set = $gene_data{$g}->{stat};
						last;
					}
				}
			}
			$cgnet_z = $candidate_paths[$seed_paths_index->list]->{z_stat_raw};
			my $P_Ni_G = double ($cgnet_z*$sd)/sqrt($candidate_paths[$seed_paths_index->list]->{N_in})  + $mu  - 1/$candidate_paths[$seed_paths_index->list]->{N_in};#$seed_stat_for_gene_set/$candidate_paths[$seed_paths_index->list]->{N_in}; 
			$P_Ni_G = ( ( $P_Ni_G - $mu) * (sqrt($candidate_paths[$seed_paths_index->list]->{N_in} - 1 ) ) )/$sd; # here are the gene-set z-scores minus the seed contribution
			$P_Ni_G = double 1-gsl_cdf_ugaussian_P($P_Ni_G); # now as p-values
			my $log10_P_G = double log10($gene_data{$seed}->{pvalue});
			# now lets calculate the p(G|N) = p(G)*p(N|G)/p(N)
			my $P_Ni = double 1-gsl_cdf_ugaussian_P($cgnet_z);
			my $log10_P_G_Ni = double $log10_P_G + $P_Ni_G->log10 - $P_Ni->log10;
			# this are the bayes factors for each gene-set: p(G|N)/p(G)
			my $log10_BF = double $log10_P_G_Ni - $log10_P_G;
			
			my $p = 1-gsl_cdf_ugaussian_P($cgnet_z);
			my $sidak_experiment_wise = 1 - (1 - $p)**$overall_cgnet_Meff;
			$seed = uc($seed);
			my $sidak_experiment_wise10_BF10_BF = "FOO";
			printf OUTCGNETS ("$seed\t%.3f\t1\t1\t$candidate_paths[$seed_paths_index->list]->{name}\t%.2e\t%.2e\t%.2e\t$candidate_paths[$seed_paths_index->list]->{name}\t%.3f\t%.3f\t%.3f\t$candidate_paths[$seed_paths_index->list]->{name}\n", $cgnet_z,$p,$p,$sidak_experiment_wise10_BF10_BF,$overall_cgnet_Meff);
		} elsif (scalar @{ $genes_paths_map{$seed}} == 0){ next; }
	}


}

print_OUT("Well Done",$LOG);
exit;

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
	
=h
	              start_1 .......................... end_1
	 start_2 .......................... end_2
	 
=cut
	
	return(1) if (($start_uno => $start_dos) and ($start_uno <= $end_dos) );
	
=h
	 start_1 .......................... end_1
	                  start_2 .......................... end_2	 
=cut
	
	return(1) if (($start_dos => $start_uno) and ($start_dos <= $end_uno) );
	# else there is no overlap
	return(0);
}

sub calculate_gene_corr_mat {
	my $genes = shift; # ARRAY ref
	my $gene_data = shift; # HASH ref
	my $gene_gene_corr = shift; # HASH reaf
	my $var = shift; # HASH ref
	my $quick_cor = shift; # 0,1
	my $max_gene_dist = shift; # integer

	my %new_corrs = ();
	my %new_vars = ();

	# define correlation matrix
	my $corr = stretcher(ones scalar @$genes); 
	
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
			if (defined $quick_cor){
				# if chromosomes are different set the correlation to 0
				if ($gene_data->{ $gn_i }->{chr} ne $gene_data->{ $gn_j }->{chr} ) {
						$gene_gene_corr->{ $gn_i }{ $gn_j } = 0;
				} else { # if are in the same chromosome
					# if user provided a maximum distance to evaluate correlations	
					if (defined $max_gene_dist){
						# if they overlap we will need to calculate the correlation
						# return 0 from check_overlap mean there is no overlap
						if (check_overlap($gene_data->{ $gn_i }->{start},$gene_data->{ $gn_i }->{end},check_overlap($gene_data->{ $gn_j }->{start},$gene_data->{ $gn_j }->{end}) == 0 ){
							
							# check the distance between the gene coordinates
							if ( $gene_data->{ $gn_i }->{start} > $gene_data->{ $gn_j }->{end}   ){
								#		                    start_i.....end_i
								#		start_j.....end_j
								if ( ($gene_data->{ $gn_i }->{start} - $gene_data->{ $gn_j }->{end}) > $max_gene_dist * 1000 ){
									$gene_gene_corr->{ $gn_i }{ $gn_j } = 0;
								}
							} elsif ( $gene_data->{ $gn_j }->{start} > $gene_data->{ $gn_i }->{end}  ){
								#		start_i.....end_i
								#							start_j.....end_j
								if ( ($gene_data->{ $gn_j }->{start} - $gene_data->{ $gn_i }->{end}) > $max_gene_dist * 1000 ){
									$gene_gene_corr->{ $gn_i }{ $gn_j } = 0;
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
	return($corr,\%new_corrs,\%new_vars);
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
				if (defined $best_x_interval){
					my $top_per_rand_invertal = shift @selected_genes;
					push @sampled_genes_interval , { 'id'=>$top_per_rand_invertal, 'n' => scalar @interval_sampling_universe };
				} else {
					foreach my $s (@selected_genes){
						push @sampled_genes_interval, { 'id'=>$s, 'n' => scalar @interval_sampling_universe };
					}
				}
			}
		}
		my @sampled_gene_stats = map { $gene_stats->{$_}->{stat}; } @sampled_genes_normal;
		foreach my $sampled_gene ( @sampled_genes_interval ){
			if (defined $best_x_interval){
				my $sidak_p = 1 - ( 1 - $gene_data{ $sampled_gene->{id} }->{pvalue} )**$sampled_gene->{n};
				$sidak_p -= 1e-15 if ($sidak_p ==1);
				$sidak_p =  $gene_data{$sampled_gene->{id}}->{pvalue} * $sampled_gene->{n} if ($sidak_p ==0);
				push @sampled_gene_stats, -1 * gsl_cdf_ugaussian_Pinv($sidak_p);
			} else {
				push @sampled_gene_stats, $gene_stats->{ $sampled_gene->{id} }->{stat};
			}
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

sub read_recombination_intervals {
        my $file = shift;
	my $gene_info = shift;
        print_OUT("Reading mapping between genes and recombination intervals from [ $file ]",$LOG);
        my %back;
        open (FH,$file) or die $!;
        while (my $line = <FH>){
                my ($chr,$start,$end,$id,$all_genes) = split(/\t/,$line);
                my @genes = map { $_=~ m/__(.*)/; } split(/;/,$all_genes);
		@genes = map { lc($_); } @genes;
                $back{$id} = {  'id'=>$id,
                                'chr'=>$chr,
                                'start'=>$start,
                                'end'=>$end,
                                'genes'=> [],
                                'N'=> scalar @genes };
                foreach my $g (@genes){
                        next if (not exists $gene_info->{$g});
                        push @{ $gene_info->{$g}->{recomb_int} }, $id;
                        push @{ $back{$id}->{genes} }, $g;
                }
        }
        return(\%back);
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
