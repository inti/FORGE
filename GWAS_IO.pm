package GWAS_IO;
use strict;
use warnings;
use Carp;
use Exporter qw (import);


# define verbose = undef;
my $v = undef;

# check the modules needed

eval { 
	use PDL;
	use PDL::Matrix;
	use PDL::NiceSlice;
	use PDL::GSL::CDF;
	use IO::File;
	use IO::Seekable;
	use Fcntl;
};
if ($@) { 
	print "Some libraries does not seem to be in you system. quitting\n";
	exit(1);
}

our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);

@EXPORT = qw( extract_stats_from_mperm_dump_all_files build_index line_with_index extract_binary_genotypes extract_genotypes_for_snp_list get_snp_list_from_bgl_format get_snp_list_from_ox_format read_bim read_fam read_map_and_ped );				# symbols to export by default
@EXPORT_OK = qw( extract_stats_from_mperm_dump_all_files build_index line_with_index extract_binary_genotypes extract_genotypes_for_snp_list get_snp_list_from_bgl_format get_snp_list_from_ox_format read_bim read_fam read_map_and_ped);			# symbols to export on request




# extract stats from plink *.mperm.dump.all files
sub extract_stats_from_mperm_dump_all_files {
	my $file = shift;
	open (DUMP,$file) or die$!;
	my $data = [];
	my $c = 0;
	while (my $line = <DUMP>){
		if ($c == 0){
			$c++;
			next;
		}
		chomp($line);
		my @stats = split(/\s+/,$line);
		my $chi = pdl @stats[1.. scalar @stats -1];
		my $z = gsl_cdf_ugaussian_Pinv( gsl_cdf_chisq_P($chi,1) );
		push @{$data},$z;
		$c++;
	}
	$data = pdl $data;
	return($data);
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

# this subroutine read the fam file and stores the information in an array
# the elements of the array are pseudo hashes with all the sample's information
sub read_fam {
	my $fam = shift;
	print "Reading samples info from [ $fam ]\n";
	open( FAM, $fam ) or print "Cannot open [ $fam ] file\n" and exit(1);
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
	#print_OUT("[ " . scalar @back . " ] samples read");
	return ( \@back );
}

# this subroutine read the bim file and store information about the SNPs
# each element of the array returned is a pseudo hash with the SNP information
sub read_bim {
	my $bim = shift;
	my $affy_to_rsid = shift;
	my $affy_id = shift;
	print "Reading SNPs info from [ $bim ]\n";
	open( BIM, $bim ) or print "Cannot open [ $bim ] file\n" and exit(1);
	my @back = ();
	while ( my $snp = <BIM> ) {
		chomp($snp);
		my @data = split( /\t+/, $snp );
		# if an affy to rsid mapping was provided change the ids
		if ( defined $affy_to_rsid ) {
			if ($data[1] !~ m/^rs/){
				if (exists $affy_id->{$data[1]}){ $data[1] = $affy_id->{$data[1]};}
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
	#print_OUT("[ " .  scalar @back . " ] SNPs on BED file");
	return ( \@back );
}

sub read_map {
	my $map = shift;
	my $affy_to_rsid = shift;
	my $affy_id = shift;
	#print_OUT("Reading SNPs info from [ $map ]");
	open( MAP, $map ) or #print_OUT("Cannot open [ $map ] file") and exit(1);
	my @back = ();
	while ( my $snp = <MAP> ) {
		chomp($snp);
		my @data = split( /\t+/, $snp );
		# if an affy to rsid mapping was provided change the ids
		if ( defined $affy_to_rsid ) {
			if ($data[1] !~ m/^rs/){
				if (exists $affy_id->{$data[1]}){ $data[1] = $affy_id->{$data[1]};}
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
	#print_OUT("  '->[ " . scalar @back . " ] SNPs on PED file");
	return ( \@back );
}

sub read_map_and_ped {
	my $file = shift;
	my $map = shift;
	my $affy_to_rsid = shift;
	my $affy_id = shift;
	
	my @bim = @{ read_map($map,$affy_to_rsid,$affy_id) };
	#print_OUT("Reading Genotypes from [ $file ]");
	open(PED,$file) or #print_OUT("Cannot open [ $file ] file") or exit(1);
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
			} else { #print_OUT(" COULD NOT RECOGNIZE THIS GENOTYPE >$g<\n" . Dumper($snp->{alleles}) . "");
				exit(1);
			}
			push @{$back_genotypes[$snp_counter]}, $recoded;
		}
		$snp_counter++;
	}
	#print_OUT("  '-> [ " . scalar @back_fam . " ] samples");
	#print_OUT("  '-> [ " . scalar @back_genotypes . " ] SNPs");
	return(\@back_fam, \@back_genotypes,\@bim);
}

sub extract_genotypes_for_snp_list{
	my $snp_list = shift;
	my $line_index = shift;
	my $g_prob_threshold = shift;
	my $geno_probs_format = shift; 
	my $gprobs = shift;
	my $gprobs_index = shift;
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
		my $start_index = 0;
		$start_index = 5 if ($geno_probs_format eq 'OXFORD');
		$start_index = 3 if ($geno_probs_format eq 'BEAGLE');
		for (my $g = $start_index; $g < scalar @genos; $g +=3){
			my $snp_prob = pdl @genos[$g..$g+2];
			my $max_index = maximum_ind($snp_prob);
			my $value = undef;
			
			if ($snp_prob->dsum == 0){ 
				$value = 0;
			} else {
				$value = $snp_prob->($max_index);
			}
			if ($value < $g_prob_threshold) {$value = 0;}
			push @{ $geno_probs[$sample_counter] } , sclr $value;
			
			my $dossage = 0*$snp_prob->(0) + 1*$snp_prob->(1) + 2*$snp_prob->(2);
			push @{ $geno_hard_coded[$sample_counter] }, sclr $dossage;
			$sample_counter++;
		}
	}
	my $coded_mat = transpose double pdl @geno_hard_coded;
	my $prob_mat = transpose double pdl @geno_probs;	
	return($prob_mat,$coded_mat);
}
sub get_snp_list_from_bgl_format {
	my $geno_probs = shift;
	my $geno_probs_index = shift;
	my $affy_to_rsid = shift;
	my $affy_id = shift;
	
	my $index = 0;
	my @back = ();
	my $desired_line = 1;
	my $eof = 0;
	while () {
		my $line = line_with_index(*$geno_probs, *$geno_probs_index, $desired_line);
		last if ($line eq '1');
		my ($snp,$a1,$a2) = split(/\s+/,$line);
		if ( defined $affy_to_rsid ) {
			if ($snp !~ m/^rs/){
				if (exists $affy_id->{$snp}){ $snp = $affy_id->{$snp};}
			}
		}		
		push @back,{
			'snp_id' => $snp,
			'chr'    => 0,
			'cm'     => 0,
			'pos'    => 0,
			'a2'     => $a1,
			'a1'     => $a2,
        };
		$desired_line++;
	}
	#print_OUT("[ " .  scalar @back . " ] SNPs on BEAGLE format genotype probability file");
	return ( \@back );
}

sub get_snp_list_from_ox_format {
	my $geno_probs = shift;
	my $geno_probs_index = shift;
	my $affy_to_rsid = shift;
	my $affy_id = shift;
	
	my $index = 0;
	my @back = ();
	my $desired_line = 1;
	my $eof = 0;
	while () {
		my $line = line_with_index(*$geno_probs, *$geno_probs_index, $desired_line);
		last if ($line eq '1');
		my ($chr,$snp,$pos,$a1,$a2) = split(/\s+/,$line);
		if ( defined $affy_to_rsid ) {
			if ($snp !~ m/^rs/){
				if (exists $affy_id->{$snp}){ $snp = $affy_id->{$snp};}
			}
		}
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
	#print_OUT("[ " .  scalar @back . " ] SNPs on OXFORD format genotype probability file");
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
		} else { 
			#print_OUT("This genotype is not recognize [ $geno ]"); 
		}    # genotype not recognize
		
	}
	return(\@back);
}
