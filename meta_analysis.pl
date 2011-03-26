#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::GSL::CDF;
use Data::Dumper;

use GWAS_STATS qw( get_lambda_genomic_control );

our (	$help,$man, $out, $w, $in_files,
		$stouffer,$w_col,$random,$stat_or,
		$stat_col,$id_col, $stat_pvalue,$header, 
		$min_studies, $gc_correction,$lambda,$sdtze,
		$gc_correction_all,$filter_col,$filter
);


GetOptions( 
	'help|h'			=> \$help,
	'man'				=> \$man,
	'out|o=s'			=> \$out, #name of the output file
	'file|f=s@'			=> \$in_files,
	'weights|w=f@'		=> \$w,
	'w_col=i@'			=> \$w_col,
	'id_col=i@'			=> \$id_col,
	'stat_col=i@'		=> \$stat_col,
	'stat_pvalue'		=> \$stat_pvalue,
	'stat_or'			=> \$stat_or,
	'header'			=> \$header,
	'min_studies|min_st=i' => \$min_studies,
	'gc_correction|gc=i@'		=> \$gc_correction,
	'gc_correction_all|gc_all'		=> \$gc_correction_all,
	'lambda=i@'			=> \$lambda,
	'sdtze=i@'			=> \$sdtze,
    'filter=s@'  => \$filter,
    'filter_col=i' => \$filter_col,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);
pod2usage(0) if (not defined $in_files);


defined $min_studies or $min_studies = -1;

if (not defined $stat_pvalue and not defined $stat_or){
	print_OUT("Please state wheter is a p-value [ -stat_pvalue ] or odd-ratio [ -stat_or ] based analysis");
	exit(1);
}

if (scalar @$in_files < 2){
	print_OUT("Please provide at least two input files");
	exit(1);
}

my %data     = ();
my @studies  = ();
my %total_seen = (); # records how many times and ids has been seen
# Defined weights
# weights can be either variable specific is provided with the option -w_col 
# or study specific if provided with the option -w
# if the option -w_col is defined then for each file the variable weight will be extracted
my @study_weights = list ones scalar @$in_files;
if (defined $w){
	@study_weights = @$w;
}

# defined gc_correction for all studies
if (defined $gc_correction_all){
	$gc_correction = [ list ones scalar @$in_files ];
}

if (defined $filter_col or  defined $filter){
    print_OUT("Applying filters [ @$filter ] to columns [ $filter_col ]");
} elsif ((defined $filter_col and not defined $filter) or (defined $filter_col and not defined $filter)) {
    print_OUT("You must privide both a filter and columns to filter on\n");
    exit(1);
}

open(OUT,">$out") or die $!;
print_OUT("Writting output to [ $out ]");


foreach my $file (@$in_files) {
	# get name of input file 
	my $st = $file;
	($st) = ( $file =~ m/[\/\w+]{1,}\/(.*)/ ) if ($file =~ m/\//);
	# get the study weights. if not provided will be 1
	my $st_w = shift @study_weights;
	# get the column number of variable specific weight 
	my $st_w_col = shift @$w_col if (defined $w_col);
	
	# get the stat_col
	my $st_stat_col = 2;
	$st_stat_col = shift @$stat_col if (defined $stat_col);
	
	# get the stat_col
	my $st_id_col = 1;
	$st_id_col = shift @$id_col if (defined $id_col);
	
	
	push @studies, $st; # store the name of this study for the print out.
	# print out some infor about the information read in the file
	print_OUT("Reading [ $file ]");
	print_OUT("   '-> Variables weights in col [ $st_w_col ]") if (defined $w_col);
	print_OUT("   '-> Variables stat in col [ $st_stat_col ]") if (defined $stat_col);
	print_OUT("   '-> Study weight [ $st_w ]") if (defined $w);

	open( IN, $file ) or print_OUT("I cannot open file [ $file ]\n") and die $!;
	my $counter = 0;
	while ( my $line = <IN> ) {
		if (defined $header and $counter == 0){
			$counter++;
			next;
		}
		chomp($line);
		# split line
		my @d = split( /[\s+\t+]/, $line );
        if (defined $filter){
            next unless (grep $_ eq $d[ $filter_col - 1], @{$filter});
        }
		# set the variable weight to the var of the study
		my $var_w = double $st_w;
		# modify the vartiable weight if a columns with its variance is given
		if ( defined $w_col){
			if ( $st_w_col > scalar @d) {
				print_OUT("Variable weight columns number is greater than number of columns: [ $st_w_col ] > [ " . scalar @d . " ]");
				print_OUT("Finishing execution\n");
				print_OUT("Problem found at line [ $. ] : >> $line <<");
				exit(1);
			}
			if ($d[ $st_w_col -1 ] != 0){
				$var_w /= double $d[ $st_w_col -1 ];
			}
		}
		$var_w = sclr $var_w; # make a perl scalar to simply calculations later
		# get the variable stat
		my $var_stat = $d[ $st_stat_col - 1 ]; 
		if (defined $stat_pvalue){ 
			if ( $var_stat == 1 ){
				$var_stat = double 1-2.2e-16;
			}
			$var_stat = -1 * gsl_cdf_ugaussian_Pinv( $var_stat );
		} elsif (defined $stat_or){
			$var_stat = log $var_stat;
		}
		# store the data as a pseudohash
		$data{ $d[ $st_id_col -1 ] }{$st} = { 'w' => $var_w, 'stat' => $var_stat };
		$total_seen{ $d[ $st_id_col -1 ] }++;
	}
	close(IN);
}

if (defined $gc_correction){
	print_OUT("Calculating lambda for genomic control correction");
	
	my %studies_p = ();
	foreach my $var (keys %data){
		my $f_counter = 0;
		foreach my $file ( @$in_files ){
			my $st = $file;
			($st) = ( $file =~ m/[\/\w+]{1,}\/(.*)/ ) if ($file =~ m/\//);
			if (not exists $data{$var}{$st}){
				$f_counter++; 
				next;
			}
			next unless ($gc_correction->[$f_counter++] == 1);
			push @{ $studies_p{$st} }, 1 - gsl_cdf_ugaussian_P($data{$var}{$st}->{stat});
		} 
	}
	foreach my $st ( keys %studies_p ){
		print_OUT(" '-> [ $st ]");
		my $gc_lambda = get_lambda_genomic_control($studies_p{$st});
		print_OUT("   '-> lambda (median) of [ $gc_lambda ]");
		if ($gc_lambda > 1){
			print_OUT("   '-> Applying GC correction in [ $st ]");
			my $assoc_p = [];
			foreach my $var (keys %data) {
				next unless (exists $data{$var}{$st});
				my $p = 1 - gsl_cdf_ugaussian_P($data{$var}{$st}->{stat});
				if ( $p == 1){
					push @{ $assoc_p }, $p;
					next;
				}
				my $var_chi = gsl_cdf_chisq_Pinv ( 1 - $p , 1 );
				$var_chi /=  $gc_lambda;
				$p = gsl_cdf_chisq_P( $var_chi, 1 );
				$data{$var}{$st}->{stat} = gsl_cdf_ugaussian_Pinv( $p );
				push @{ $assoc_p }, 1 - $p;
			}
			$gc_lambda = get_lambda_genomic_control($assoc_p);
			print_OUT("   '-> After correction the lambda is [ $gc_lambda ]");
		} else {
			print_OUT("   '-> GC correction not applied because lambda is less than 1");
		}
	}
}

print_OUT("Calculating Stats");
print OUT join "\t", ("id","stouffer_fix","B_fix","Var_fix","Chi-square_df1_fix","B_fix_P","stouffer_fix","B_random","Var_random","Chi-square_df1_random","B_random_P","Q","Q_P","tau_squared","I2","N","binomial_rank_pvalue",@studies);
print OUT "\n";
foreach my $gn ( keys %data ) {
	my $count = scalar keys %{$data{ $gn } };
	next if ($count < 2); # next if only one study looked at this gene
	next if ($count < $min_studies);
	next if ($total_seen{ $gn } > scalar @studies);

	my $all_w = [];
	my $all_b = [];
	# collect stat and weights for each study
	foreach my $st ( keys %{ $data{ $gn } } ) {
		push @{$all_w}, $data{ $gn }{$st}->{'w'};
		push @{$all_b}, $data{ $gn }{$st}->{'stat'};
	}
	$all_w = pdl $all_w;
	if ($all_w->minimum() == 1 and $all_w->maximum() ==1){
		$all_w *= $all_w->nelem;
	} else {
		$all_w = 1/$all_w;
	}
	my $metaResults = get_fix_and_radom_meta_analysis($all_b,$all_w);
	
	my $max_b = maximum pdl $all_b;
	my $min_p = gsl_cdf_ugaussian_P( -1 * $max_b );
	my $binomial_rank_p = 1 - gsl_cdf_binomial_P(1,$min_p,$metaResults->{'N'});

	printf OUT ("$gn\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6f\t%.6e\t%.6f\t%.6f\t$metaResults->{'N'}\t%.6e", 
				$metaResults->{'B_stouffer_fix'},
				$metaResults->{'B_fix'},
				$metaResults->{'V_fix'}, 
				$metaResults->{'Chi_fix'},
				$metaResults->{'Z_P_fix'},
				$metaResults->{'B_stouffer_random'},
				$metaResults->{'B_random'},
				$metaResults->{'V_random'}, 
				$metaResults->{'Chi_random'},
				$metaResults->{'Z_P_random'},
				$metaResults->{'Q'},
				$metaResults->{'Q_P'},
				$metaResults->{'tau_squared'},
				$metaResults->{'I2'},
				$binomial_rank_p,
				);
	map {
		if (defined $w_col){
			if   ( exists $data{ $gn }{$_} ) { printf OUT ("\t%.6f(%.3f)",$data{ $gn }{$_}->{'stat'},1/$data{ $gn }{$_}->{'w'}); }
			else { print OUT "\tNA"; }			
		} else {
			if   ( exists $data{ $gn }{$_} ) { printf OUT ("\t%.6f",$data{ $gn }{$_}->{'stat'}); 
				#print Dumper($data{ $gn }),"\n";
			}
			else { print OUT "\tNA"; }
		}
	} @studies;
	print OUT "\n";
	delete($data{ $gn });
}

print_OUT("Done baby");

exit;
sub get_fix_and_radom_meta_analysis {
	my $B = shift;
	my $SE = shift;
	my $external_w = shift;
	my $VarCov = shift;
	if (ref($B) eq 'ARRAY') { $B = pdl $B; }
	if (ref($SE) eq 'ARRAY') { $SE = pdl $SE; }
	if (ref($external_w) eq 'ARRAY') { $external_w = pdl $external_w; }
	my $N = $B->nelem;
	if (not defined $VarCov) {
		$VarCov = stretcher(ones $N);
	}
	my $W = null;
	if ( (ref($external_w) eq 'PDL') and (nelem($external_w > 0) > 0) and (dsum( ($external_w - $external_w->davg)**2 ) > 0 )){
		$W = 1/($SE + 1/$external_w); 
	} else {
		$W = 1/$SE; 
	}

	# calculate fix effect estimate
	my $norm_w_fix = $W/$W->dsum;
	#my $B_fix = dsum($B*$W)/$W->dsum;
	#my $V_fix = dsum($W*$W->transpose*$VarCov);
	my $B_fix = dsum($B*$norm_w_fix)/$norm_w_fix->dsum;
	my $V_fix = dsum($norm_w_fix*$norm_w_fix->transpose*$VarCov);
	my $fix_chi_square_df1 = ($B_fix**2)/$V_fix;
	
	my $stouffer_w_fix = $W/dsum($W**2);
	my $B_stouffer_fix = dsum($B*$stouffer_w_fix)/sqrt(dsum($stouffer_w_fix*$stouffer_w_fix->transpose*$VarCov));
	
	# calculate heteroogeneity parameter Q
	my $Q = 0.0;
	{
		my $cor = $VarCov->copy();
		$cor->diagonal(0,1) .=0;
		$cor = (3.263*abs($cor) + 0.710*(abs($cor)**2) + 0.027*(abs($cor)**3));
		my $df = 8/(dsum($cor*$W*$W->transpose) + 4*dsum($W**2));
		my $Q_naive = dsum ($W * (($B_fix - $B)**2));
		$Q = sclr gsl_cdf_chisq_Pinv( gsl_cdf_chisq_P(  $df*0.5*$Q_naive, $df), $N - 1);
	}	
	# calculate heteroogeneity parameter I-squared
	my $I_squared = 0.0;
	if ($Q > ($N - 1)){
		eval { $I_squared = 100*($Q - ($N - 1))/$Q; };
		if ($@){
			$I_squared = 0.0;
		}
	}
	# calculate tau-squared
	my $tau_squared = 0.0;
	if ($Q > ($N - 1)){
		eval { $tau_squared = ($Q - ($N - 1))/(dsum($W) - dsum($W**2)/$W->dsum); };
		if ($@){
			# tau-squared equal 0 if was < 0
			$tau_squared = 0.0;
		}
	}
	# calculate the random effect estimate
	my $w_star = null;
	if ( (ref($external_w) eq 'PDL') and (nelem($external_w > 0) > 0) and (dsum( ($external_w - $external_w->davg)**2 ) > 0 )){
		$w_star = 1/( $tau_squared + $SE + 1/$external_w);
	} else {
		$w_star = 1/( $tau_squared + $SE);
	}
	my $norm_w_random = $w_star/$w_star->dsum;
	my $B_random = dsum( $B*$norm_w_random)/$norm_w_random->dsum;
	#my $B_random = dsum( $B*$w_star)/$w_star->dsum;
	my $V_random = dsum($norm_w_random*$norm_w_random->transpose*$VarCov);
	#	my $V_random = dsum($w_star*$w_star->transpose*$VarCov);
	my $random_chi_square_df1 = ($B_random**2)/$V_random;
	
	my $stouffer_w_random = $w_star/dsum($w_star**2);
	my $B_stouffer_random = dsum($B*$stouffer_w_random)/sqrt(dsum($stouffer_w_random*$stouffer_w_random->transpose*$VarCov));
	
	my $back =  {
		'B_stouffer_fix' => $B_stouffer_fix,
		'B_stouffer_random' => $B_stouffer_random,
		'B_fix' => $B_fix,
		'B_random' => $B_random,
		'V_fix' => $V_fix,
		'V_random' => $V_random,
		'Chi_fix' => $fix_chi_square_df1,
		'Chi_random' => $random_chi_square_df1,
		'Q' => $Q,
		'I2' => $I_squared,
		'tau_squared' => $tau_squared,
		'N' => $N,
		'Z_P_fix' => gsl_cdf_gaussian_P( -1 * $B_fix ,sqrt($V_fix )),
		'Z_P_random' => gsl_cdf_gaussian_P( -1 * $B_random,sqrt($V_random) ),
		'Q_P' => 1- gsl_cdf_chisq_P($Q ,$N -1),
	};
	return $back;
}

sub print_OUT {
	my $string = shift;
	print scalar localtime(), "\t$string\n";
	#	print LOG scalar localtime(), "\t$string\n";
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
 
 
 
 =head1 OPTIONS
 
 =over 8
 
 =item B<-help>
 
 Print help message
 
 =item B<-man>
 
 print complete documentation
 
 
 =back
 
 =head1 DESCRIPTION
 
 TODO
 
 
 =cut
