package Pareto_Distr_Fit;
use strict;
use diagnostics;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
use PDL::GSL::CDF;
use PDL::Stats::Basic;
use PDL::Stats::Distr;
use PDL::LinearAlgebra qw(mchol);
use PDL::GSL::INTERP;

our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);

@EXPORT = qw( Pgpd );				# symbols to export by default
@EXPORT_OK = qw( Pgpd );			# symbols to export on request


sub Pgpd {
=h
	 Computing permutation test P-value of the GPD approximation
		
		 Input:    x0         original statistic
				   y          permutation values
				   N          number of permutation values
				   Nexc       number of permutations used to approximate the tail
				   alpha      confidence level
		 Output:   Phat       estimated P-value
				   Phatci     confidence intervals of the estimated P-value
	 Original code in Matlab by Theo Knijnenburg, Institute for Systems Biology (Dec 11 2008)
=cut

	my $x0 = shift;
	my $y =  shift;
	my $Nexc = shift;
	my $alpha = shift;
	my $N = $y->nelem;
	my $Padth = 0.05;

	if ($Nexc >= $N-1) {
		$Nexc = int(0.1*$N);
	}
	
	#print "$N $Nexc $alpha\n";
	
	my ($PW2_less_005,$PA2_less_005,$W2,$A2,$cov,$a,$k,$t,$frac) = undef;
	my $fitted = 0;
	do {
		$y = pdl reverse list qsort $y;
		# Defining the tail	
		my $z = $y->(1:$Nexc);
		$t = davg($y($Nexc:$Nexc+1));
		$z = $z-$t;
		
		$frac = $Nexc/$N;
		# Fitting the tail and computing the Pvalue
		($a,$k,$cov) = Jin_ZHANG_pareto_fit($z);
		#print "$Nexc $a $k $cov\n";
		my $p = 1- gsl_cdf_pareto_P($z,$a,$k);
		($PW2_less_005,$PA2_less_005,$W2,$A2) = gpdgoft( $p ,$k); # goodness-of-fit test
		#print "$PW2_less_005,$PA2_less_005,$W2,$A2\n";
		
		if ( $A2 eq 'nan'){
			$fitted = $PA2_less_005;
		} else {
			$fitted = $PW2_less_005;
		}
		$Nexc -= 10;
	} until ( ( $fitted == 1 ) or ($Nexc < 11) );	
	if ($fitted == 1){
		my ($Phat,$Phatci_low,$Phatci_up) = gpdPval($x0-$t,$a,$k,$cov,$alpha);
		$Phat = $frac*$Phat;
		$Phatci_up = $frac*$Phatci_up;
		$Phatci_low =$frac*$Phatci_low; 
		return($Phat,$Phatci_low,$Phatci_up);
	} else {
		return(-1,-1,-1);
	}
}

sub gpdPval {
	my $x0 = shift;
	my $a = shift;
	my $k = shift;
	my $cov = shift;
	my $alpha = shift;


	my $Phat = 1 - gsl_cdf_pareto_P($x0,$a,$k);
	
	my $gaussian_random_mat = grandom 2, 10000; 
	my $reapmat= ones 2, 10000;
	my $par = pdl $a,$k;
	$reapmat = $reapmat * $par; 

	my $chol = mchol($cov);
	my $Q = ($gaussian_random_mat  x $chol->transpose ) + $reapmat;
	my $sample =  1 - gsl_cdf_pareto_P($x0,$Q->(0,),$Q->(1,));
	$sample->inplace->setnantobad;
	$sample =  $sample->where(  $sample->isgood() );
	
	my $Phatci_low = pct($sample,$alpha/2);
	my $Phatci_up =  pct($sample,1-$alpha/2);	
	
	return($Phat,$Phatci_low,$Phatci_up);	
}



sub Jin_ZHANG_pareto_fit {
	my $x = shift; # x is the sample data from the GPD
	my $n = $x->nelem;
	$x = qsort $x;
	my $p = pdl (3 .. 9);
	$p /= 10;
	my $xp = $x->(rint($n*(1-$p)-0.5));
	my $m = 20 + rint( sqrt($n) );
	my $xq = $x->(rint($n*(1-$p*$p)+0.5));
	my $k = log($xq/$xp-1)/log($p);
	$k->inplace->setnantobad;
	$k->inplace->setbadtoval(0);
	my $a = $k*$xp/(1-$p**$k);
	$a->inplace->setnantobad;
	$a->inplace->setbadtoval(0);	
	$a->where( $k == 0 ) .= (-$xp->where( $k == 0 )/log($p->where( $k == 0 )) ) ;
	$k = -1;
	my $tmp1 = pdl 1 .. $m;
	my $j_minus05_overm = ($tmp1 -0.5)/$m;
	my $leading = (($n-1)/($n+1))*(1/$x->($n-1));
	my $b =  $leading - ((1- $j_minus05_overm**$k)/$k) * (1/ median($a))/2;
	my @tmp = ();
	for (my $i = 0; $i < $m; $i++) { 
		push @tmp, $n * lx($b->($i),$x);
	}
	my $L = flat pdl @tmp;
	@tmp = ();
	for (my $i = 0; $i < $m; $i++) {
		push @tmp, 1/dsum( exp( $L - $L->( $i) ) );
	}
	my $w = flat pdl @tmp;
	$b = dsum($b*$w);
	$k = -davg(log(1-$b*$x));
	my $sigma = $k/$b;
	
	my $cov = mones 2,2;
	$cov(0,0) .= 2*($sigma**2)*(1-$k);
	$cov(1,0) .= $sigma*(1 - $k);
	$cov(0,1) .= $sigma*(1 - $k);
	$cov(1,1) .= (1 - $k)**2;
	return($sigma, $k,$cov);
}

sub lx {
	my $b= shift;
	my $x = shift;
	
	my $k = - davg( log( 1-$b*$x ) );
	if ($b == 0) {
		return( $k-1-log(davg($x)) );
	} else {
		return( $k-1+log($b/$k) );
	} 
}	




sub gpdgoft{
	my $p = shift;
	my $k = shift;
	
	$p = qsort $p;
	$p->where( (1 - $p) < 2.2e-16 ) .= 2.2e-16;

	my ($W2,$A2) =  fit_stats($p);
		
	my $W2_less_005 = 0; # 0 no and 1 yes
	my $A2_less_005 = 0; # 0 no and 1 yes
	
	$k = 0.5 if ($k > 0.5);
	
	my $k_vals = pdl (-0.9,-0.5,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5);
	my $W2_vals  = pdl (0.115,0.124,0.137,0.114,0.153,0.16,0.171,0.184,0.201,0.222);
	my $A2_vals = pdl (0.771,0.83,0.903,0.935,0.974,1.02,1.074,1.14,1.221,1.321);
	
	my $trend_A2 =  PDL::GSL::INTERP->init('linear',$k_vals,$A2_vals);
	my $trend_W2 =  PDL::GSL::INTERP->init('linear',$k_vals,$W2_vals);
	
	
	my $critic_w2 = undef;
	my $critic_a2 = undef;
	
	if ($k < -0.9 or $k > 0.5) {
		$critic_w2 = $trend_W2->eval($k,{ Extrapolate => 1 } );
		$critic_a2 = $trend_A2->eval($k,{ Extrapolate => 1 } );
	} else {
		$critic_w2 = $trend_W2->eval($k,{ Extrapolate => 0 } );
		$critic_a2 = $trend_A2->eval($k,{ Extrapolate => 0 } );
	}

	$A2_less_005 = 1 if ( $A2 >= $critic_a2 );
	$W2_less_005 = 1 if ( $W2 >= $critic_w2 );
	
	return($A2_less_005,$W2_less_005,$W2,$A2);
}

sub fit_stats { 
	my $p = shift;
	my $n = $p->nelem;
	my $i = 1 + sequence $n;
	
	# Cramer-von Mises statistic
	my $W2 = dsum( ($p - ( (2*$i-1)/(2*$n) ))^2 ) + 1/(12*$n);
	# Anderson Darling statistic
	my $A2 = -$n -  (1/$n) * dsum( ((2*$i-1) * (log($p)) + log(double 1 - double $p->($n+1-$i - 1)) ));
	return($W2,$A2);
}

1;