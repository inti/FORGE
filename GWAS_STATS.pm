package GWAS_STATS;
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
	use PDL::GSL::CDF;
	use PDL::Primitive;
	use PDL::NiceSlice;
	use PDL::Stats::Basic;
	#use PDL::LinearAlgebra; # commented until re-implement the use of simulation to calculate p-values
	use PDL::Bad;
};
if ($@) { 
	print "Some libraries does not seem to be in you system. quitting\n";
	exit(1);
}

our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);

@EXPORT = qw( get_makambi_chi_square_and_df calculate_LD_stats get_lambda_genomic_control number_effective_tests z_based_gene_pvalues);	# symbols to export by default
@EXPORT_OK = qw(  get_makambi_chi_square_and_df calculate_LD_stats get_lambda_genomic_control number_effective_tests z_based_gene_pvalues); # symbols to export on request


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
	my $chi_stat = dsum(-2 * $pvalues->log * $w); # and the chi-square for the combine probability
	$chi_stat = ( $chi_stat/2 ) * $df;
	my $sum = dsum(-2 * $pvalues->log * $w);
	if (defined $v){
		print_OUT("df: $df; var: $varMf_m; w: $w; chi-stat:$chi_stat; SUM: $sum");
		print_OUT("P-values: " . $pvalues . "");
	}   
	return ($chi_stat,$df);
}


sub z_based_gene_pvalues {
	my $gene = shift;
	if ( ref($gene->{'pvalues'}) eq 'ARRAY'){ $gene->{'pvalues'} = pdl @{ $gene->{'pvalues'} }; }
	if ( ref($gene->{'effect_size'}) eq 'ARRAY'){ $gene->{'effect_size'} = pdl @{ $gene->{'effect_size'} }; }
	if ( ref($gene->{'effect_size_se'}) eq 'ARRAY'){ $gene->{'effect_size_se'} = pdl @{ $gene->{'effect_size_se'} }; }
	
	if (scalar @{$gene->{'geno_mat_rows'}} < 2){
		return(-9);
	}
	my $cov = $gene->{cor_ld_r}**2;
	
	my $pvals = $gene->{'pvalues'};
	$pvals->index( which($pvals == 1) ) .= double 1-2.2e-16;
	$pvals->index( which($pvals == 0) ) .= double 2.2e-16;
	my $B = -1*gsl_cdf_ugaussian_Pinv($pvals);
	my $observed_stat = undef;
	if (not defined $gene->{'effect_size_se'}){
		$observed_stat = get_fix_and_radom_meta_analysis($B,$gene->{'effect_size_se'},undef,$cov);
	} else {
		my $se = 1/$gene->{'weights'};
		$observed_stat = get_fix_and_radom_meta_analysis($B,$se,undef,$cov);
	}
	return($observed_stat);	
}


sub get_lambda_genomic_control {
	my $p = shift;
	$p = double 1 - pdl $p;
	my $chi = gsl_cdf_chisq_Pinv($p,1);
	return $chi->median/0.456;
}



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
	my $B_fix = dsum($B*$W)/$W->dsum;
	my $norm_w_fix = $W/$W->dsum;
	my $V_fix = dsum($norm_w_fix*$norm_w_fix->transpose*$VarCov);
	#my $V_fix = dsum($W*$W->transpose*$VarCov);
	my $fix_chi_square_df1 = ($B_fix**2)/$V_fix;
	my $fix_Z = $B_fix/sqrt($V_fix);
	
	my $B_stouffer = dsum($B*$norm_w_fix)/sqrt($V_fix);
	
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
	my $B_random = dsum( $B*$w_star)/$w_star->dsum;
	my $norm_w_random = $w_star/$w_star->dsum;
	my $V_random = dsum($norm_w_random*$norm_w_random->transpose*$VarCov);
	#	my $V_random = dsum($w_star*$w_star->transpose*$VarCov);
	my $random_chi_square_df1 = ($B_random**2)/$V_random;
	my $random_Z = $B_random/sqrt($V_random);
	
	my $back =  {
		'B_stouffer' => $B_stouffer,
		'B_fix' => $B_fix,
		'B_random' => $B_random,
		'V_fix' => $V_fix,
		'V_random' => $V_random,
		'Chi_fix' => $fix_chi_square_df1,
		'Chi_random' => $random_chi_square_df1,
		'Z_fix' => $fix_Z,
		'Z_random' => $random_Z,
		'Q' => $Q,
		'I2' => $I_squared,
		'tau_squared' => $tau_squared,
		'N' => $N,
	};
	return $back;
}

sub calculate_LD_stats {
=head1 Docs
	 
	 From Lon Cardon
	 
	 The EM algorithm is used to estimate the recombination fraction (parameter theta)
	 between two genetic markers under the assumption of Hardy-Weinberg equilibrium.
	 
	 For SNPs, consider the 3 x 3 table of marker 1 with alleles 'A' and 'a' 
	 and marker 2 with alleles 'B' and 'b'.  Number the cells as follows:
	 
	 AA  |  Aa  |  aa
	 -------------------
	 BB |  0  |   1  |  2
	 Bb |  3  |   4  |  5
	 bb |  6  |   7  |  8
	 
	 Note that all cells arise from unique haplotypes (e.g., cell 0 can only comprise
	 two 'AB' haplotypes; cell 1 has one 'AB' and one 'aB', etc), _except_ the double
	 heterozygote cell 4, which can have either AB/ab or Ab/aB.  We need to estimate the
	 probability of these two events.  We assume that pair AB/ab arises when no
	 recombination occurs (1 - theta) whereas pair Ab/aB arises in the presence of recombination
	 (probability theta).  Thus, for the four haplotype possibilities at two markers 
	 (AB,ab,Ab,aB), we have:
	 
	 prAB=prAb=praB=prab=.25;  
	 
	 nAB=(float)(2*cells[0]+cells[1]+cells[3]); 
	 nab=(float)(2*cells[8]+cells[7]+cells[5]);
	 nAb=(float)(2*cells[6]+cells[7]+cells[3]);
	 naB=(float)(2*cells[2]+cells[1]+cells[5]);
	 
	 N = nAB + nab + nAb + naB + 2*cell4;
	 
	 while(theta-thetaprev > CONVERGENCE_CRITERIA) 
	 {
	 thetaprev=theta;
	 prAB=(nAB + (1-theta)*cells[4])/N;
	 prab=(nab + (1-theta)*cells[4])/N;
	 prAb=(nAb + theta*cells[4])/N;
	 praB=(naB + theta*cells[4])/N;
	 theta=(prAb*praB)/(prAB*prab + prAb*praB); 
	 }
	 
	 
	 i.e., D = prAB - frq(A)*frq(B), r^2 = D^2/(frq(A)*frq(B)*frq(a)*frq(b)),
	 Dprime = D/Dmax, etc.
	 
=cut
    my $first  = shift;
    my $second = shift;
	
    my %genotypes = (	'AABB' => 0, # cell 0
	'AaBB' => 0, # cell 1
	'aaBB' => 0, # cell 2
	
	'AABb' => 0, # cell 3
	'AaBb' => 0, # cell 4
	'aaBb' => 0, # cell 5
	
	'AAbb' => 0, # cell 6
	'Aabb' => 0, # cell 7
	'aabb' => 0  # cell 8
	);
	
	for (my $i = 0; $i < scalar @$first; $i++){
		$genotypes{'AABB'}++ if (( $first->[$i] == 1) and ($second->[$i] == 1));	
		$genotypes{'AaBB'}++ if (( $first->[$i] == 2) and ($second->[$i] == 1));	
		$genotypes{'aaBB'}++ if (( $first->[$i] == 3) and ($second->[$i] == 1));	
		
		$genotypes{'AABb'}++ if (( $first->[$i] == 1) and ($second->[$i] == 2));	
		$genotypes{'AaBb'}++ if (( $first->[$i] == 2) and ($second->[$i] == 2));	
		$genotypes{'aaBb'}++ if (( $first->[$i] == 3) and ($second->[$i] == 2));	
		
		$genotypes{'AAbb'}++ if (( $first->[$i] == 1) and ($second->[$i] == 3));	
		$genotypes{'Aabb'}++ if (( $first->[$i] == 2) and ($second->[$i] == 3));	
		$genotypes{'aabb'}++ if (( $first->[$i] == 3) and ($second->[$i] == 3));	
	}
	
	my ($nAB,$nAb,$naB,$nab) = double 0.0;
    $nAB = 2*$genotypes{'AABB'} + $genotypes{'AaBB'} + $genotypes{'AABb'};
    $nab = 2*$genotypes{'aabb'} + $genotypes{'Aabb'} + $genotypes{'aaBb'};
    $nAb = 2*$genotypes{'AAbb'} + $genotypes{'Aabb'} + $genotypes{'AABb'};
    $naB = 2*$genotypes{'aaBB'} + $genotypes{'AaBB'} + $genotypes{'aaBb'};
    my $AaBb = double $genotypes{'AaBb'};
	
 	
	if (defined $v){
		print "Haplotype freqs\n";
		print  "Table: $genotypes{'AABB'}  $genotypes{'AaBB'} $genotypes{'aaBB'}\n";
	    print  "Table: $genotypes{'AABb'}  $genotypes{'AaBb'} $genotypes{'aaBb'}\n";
	    print  "Table: $genotypes{'AAbb'}  $genotypes{'Aabb'} $genotypes{'aabb'}\n";
	}
	
    my $N = double $nAB + $nab + $nAb + $naB + 2*$AaBb;
    my $theta = 0.5;
    my $thetaprev = 2;
   	my $iterations = 0;
    while( abs($theta-$thetaprev) > 0.0001 ) {	
		$thetaprev = $theta;
		eval{
			$theta = (($nAb + $theta*$AaBb)*($naB + $theta*$AaBb))/
		    (($nAB + (1-$theta)*$AaBb)*($nab + (1-$theta)*$AaBb) + 
			($nAb + $theta*$AaBb)*($naB + $theta*$AaBb));
		};
		if ($@){ $theta = 0.5; } #included to avoid division by 0
		$iterations++;
    }
	
    # now calculate stats
    #my ( $f_A, $f_B ) = major_freqs( \%people );
	my ( $f_A, $f_B ) = 0;
	my $n = 0;
	map {  
		if ($_ > 0){
			$f_A += 2 if ($_ == 1);
			$f_A += 1 if ($_ == 2);
			$n++;
		}
	} @$first;
    $f_A /= 2*$n;
	$n = 0;
	map {  
		if ($_ > 0){
			$f_B += 2 if ($_ == 1);
			$f_B += 1 if ($_ == 2);
			$n++;
		}
	} @$second;
    $f_B /= 2*$n;
	
	#	print "$f_A $f_B\n";
	#getc;
    my $D;
    my $r2;
	my $r;
    eval{
		$D  = ($nAB+(1-$theta)*$AaBb)/$N - ($f_A*$f_B);
		$r2 = $D*$D/($f_A*$f_B*(1-$f_A)*(1-$f_B)); 
		$r = $D/sqrt($f_A*$f_B*(1-$f_A)*(1-$f_B));
    };
	
    if ($@){
		$D = 0;
		$r2 = 0; #for some cases is not possible to calculate the r2 due to a 0 in the divisor
		$r = 0;
    }
    
    my $Dmax = 0;
    my $d_prime;
    
    if ($D < 0){
		$Dmax = $f_A*$f_B if ($f_A*$f_B < ((1-$f_A)*(1-$f_B)));
		$Dmax = (1-$f_A)*(1-$f_B) if ($f_A*$f_B >= ((1-$f_A)*(1-$f_B)));	
    }
    if ($D > 0){
		$Dmax = $f_A*(1-$f_B) if ($f_A*(1-$f_B) < (1-$f_A)*$f_B);
		$Dmax = (1-$f_A)*$f_B if ($f_A*(1-$f_B) >= (1-$f_A)*$f_B);
    }
    eval{
		$d_prime = $D/$Dmax;
    };
    if ($@){
		$d_prime = 0;
    }
	
	$D = sclr $D;
	$r2 = sclr $r2;
	$r = sclr $r;
	$N = sclr $N;
	$d_prime = sclr $d_prime;
	$theta = sclr $theta;
    my $o = { 'D'=> $D,
		'r2' => $r2,
		'r' => $r,
		'theta' => $theta,
		'N' =>  $N,
		'd_prime' =>  $d_prime,
		#'people' => scalar(keys %people)
	};
    return $o;
	
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
		#print_OUT("  AA Aa aa");
		#print_OUT("BB $cells[0]  $cells[1]  $cells[2]");
		#print_OUT("Bb $cells[3]  $cells[4]  $cells[5]");
		#print_OUT("bb $cells[6]  $cells[7]  $cells[8]");
	}
	my $h1 = sqrt($cells[0]);
	my $h2 = sqrt($cells[2]);
	my $h3 = sqrt($cells[6]);
	my $h4 = sqrt($cells[8]);
	my $corr = 2*($h1 * $h4 - $h2 * $h3)/sqrt(4*($h1 + $h2)*($h3 + $h4)*($h1 + $h3)*($h2 + $h4));
	return($corr);
}

	
# this subroutine calculate the number of effective test by the Galwey and Gao method.
sub number_effective_tests {
	my $mat = shift;
	# calculate the eigen value of the correlation matrix
	my $eigens = eigens ${$mat};
	# normalize the values
	my $eigens_norm =  pdl sort { $b <=> $a } list ($eigens/(dsum $eigens) );
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
	if (defined $v){ 
		#print_OUT(" simpleM_Gao = $simpleM; Meff_galwey=$Meff_galwey $numerator $denominator"); 
	}
	return($simpleM,$Meff_galwey);
}

1;