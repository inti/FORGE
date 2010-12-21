#!/usr/bin/perl -w
use strict;
use PDL;
use PDL::Matrix;
use PDL::Primitive;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::LinearAlgebra::Real;

sub wt_moments{
	my ($x, $w) = @_; 
	$w = pvt_check_w($w, $x->getdim(1));
	my $h1 = 1/(1 - dsum($w * $w));
	my $m = dsumover($w * $x->transpose);
	my $v = $h1 * (dsumover($w * transpose($x**2)) - dsumover($w * $x->transpose)**2);
	return({ 'mean' => $m, 'var' => $v});
}
sub wt_scale{
	my $x = shift;
	my $w = shift;
	my $center = shift;
	my $scale = shift;
	$w = pvt_check_w($w,$x->getdim(1));
    my $wm = { 'mean' => $x->xchg(0,1)->daverage, 'var' => $x->xchg(0,1)->stdv**2 }; 
    if ($center eq 'TRUE'){
        $x = $x - $wm->{'mean'};
    }
    if ($scale eq 'TRUE') {
        my $sc = sqrt($wm->{'var'});
        $x = $x/$sc;
    }
    return($x);
}
sub pvt_check_w {
	my ($w, $n) = @_; 
    if (not defined $w) {
        $w = ones $n;
		$w /= $n;
		$w /= $w->dsum;
    } else {
        if ($w->getdim(0) != $n) {
            warn("Weight vector has incompatible length. Is [ " . $w->getdim(0) . " ] should be [ $n ]\n");
        }
        $w = ones $n;
		$w /= $n;
		$w /= $w->dsum;
        my $sw = sum($w);
        if ($sw != 1){ $w = $w/$sw; }
	}
    return($w);
}

sub make_positive_definite {

	my ($m, $tol) = @_; 
	my $d = $m->getdim(1);
	if ($m->getdim(1) != $m->getdim(0)) { die("Input matrix is not square!"); }
	my ($es,$esv) = eigens $m;
	if (not defined $tol){ $tol = $d * max(abs($esv)) * 2e-16; }
	my $delta = 2 * $tol;
	my $tau = $delta - $esv;
	$tau->( which($tau < 0) ) .= 0;
	my $dm = $es x stretcher($tau) x transpose($es);
	return($m + $dm)
}

=head
sub is_positive_definite
function (m, tol, method = c("eigen", "chol")) 
{
    method = match.arg(method)
    if (!is.matrix(m)) 
        m = as.matrix(m)
		if (method == "eigen") {
			eval = eigen(m, only.values = TRUE)$values
			if (is.complex(eval)) {
				warning("Input matrix has complex eigenvalues!")
				return(FALSE)
			}
			if (missing(tol)) 
				tol = max(dim(m)) * max(abs(eval)) * .Machine$double.eps
				if (sum(eval > tol) == length(eval)) 
					return(TRUE)
					else return(FALSE)
						}
    if (method == "chol") {
        val = try(chol(m), silent = TRUE)
        if (class(val) == "try-error") 
            return(FALSE)
			else return(TRUE)
				}
}
=cut

sub cov_shrink {
	my ($x, $lambda, $lambda_var, $w, $collapse, $verbose) = @_;
    my $n = $x->getdim(1);
	if (not defined $lambda) { $lambda = -1; }
	if (not defined $lambda_var) { $lambda_var = -1; } 
	if (not defined $w) { 
		$w = ones $n;
		$w /= $n;
		$w /= $w->sum;
	}
	my $sc_data = pvt_svar($x, $lambda_var, $w, $verbose);
	my $sc = $sc_data->{'vs'}->sqrt;
	my $c = pvt_powscor($x, 1, $lambda, $w, $collapse, $verbose);
	my $cor = $c->{powr};
	# I have commented this if because I do not know when is useful. 
	# it seems to be some kind of error checking that i do not understand.
	#if ( $c->isempty){ 
	#	$c->{powr} = $c->{powr} * $sc * $sc;
	#} else{
		$c->{powr} = transpose($c->{powr}*$sc)* $sc;
	#}
	my $back = { 
		'lambda_var' => $sc_data->{"lambda_var"},
		'lambda_var_estimated' => $sc_data->{"lambda_var_estimated"},
		'lambda_cor' => $c->{"lambda"},
		'lambda_cor_estimated' => $c->{"lambda_estimated"},
		'cov' => $c->{powr},
		'cor' => $cor,
	};
	return($back);
}



sub pvt_get_lambda {
	my ($x, $lambda, $w, $verbose, $type,$target) = @_;
    my ($kind,$func) = undef;
	my $lambda_estimated = undef;
    
    if ($lambda < 0) {
		if ($type eq "correlation") {
			$kind = "lambda (correlation matrix):";
			$lambda = pvt_corlambda($x, $w, $target);
		}
		if ($type eq "variance") {
			$kind = "lambda.var (variance vector):";
			$lambda = pvt_varlambda($x, $w, $target);
		}
        if ($verbose) {
            print  ("Estimating optimal shrinkage intensity [ $kind ]");
        }
		$lambda_estimated = "TRUE";
        if ($verbose) {
            print  ("lambda [ $lambda ]");
        }
    } elsif ($lambda > 1){  
		$lambda = 1;
		$lambda_estimated = "FALSE";
		if ($verbose) {
			print  ("Specified shrinkage intensity $kind [ $lambda ]");
		}
    }
    if ($type eq "correlation") {
        return( { 'lambda' => $lambda, 'lambda_estimated' => $lambda_estimated } );
    }
    if ($type eq "variance") {
        return( { 'lambda_var' => $lambda, 'lambda_var_estimated' => $lambda_estimated } );
    }
}

sub pvt_powscor {
	my ($x, $alpha, $lambda, $w, $collapse, $verbose) = @_;
    my $xs = wt_scale($x, $w, 'TRUE', 'TRUE');
	my $h1 = 1/(1 - dsum($w * $w));
    my $z = pvt_get_lambda($xs, $lambda, $w, $verbose, "correlation", 0);
    my $p = $xs->getdim(0);
	my $powr = null;
	# for debug
	#	$alpha = 0.5;
    if ($z->{'lambda'} == 1 or $alpha == 0) {
        if ($collapse) {
            $powr = ones $p;
        } else {
            $powr = zeroes $p,$p;
			$powr->diagonal(0,1) .= ones $p;
        }
    } elsif ($alpha == 1) {
		my $cross = $xs*$w->transpose->sqrt;
        my $r0 = $h1 * mcrossprod( $xs*$w->transpose->sqrt);
        $powr = (1 - $z->{'lambda'}) * $r0;
        $powr->diagonal(0,1) .= ones $powr->getdim(0);
    } else {
		# These methods have not been implemented
		# program will throw an error.
        my $zeros = $xs;
        my $svdxs = fast_svd($xs);
        my $m = length($svdxs->{'d'});
        my $UTWU = transpose($svdxs->{'u'}) x ( $svdxs->{'u'}*$w);
        my $C = $UTWU*$svdxs->{'d'}*$svdxs->{'u'}->transpose;
        $C = (1 - $z->{'lambda'}) * $h1 * $C;
        $C = ($C + transpose($C))/2;
        if ($lambda == 0) {
            if ($m < $p - dsum($zeros)){ 
                warning("Estimated correlation matrix doesn't have full rank - pseudoinverse used for inversion.\n");
				$powr = $svdxs->{'v'} x ((mpower($C, $alpha) x transpose($svdxs->{'v'})) );
			}
		} else {
            my $F = $m->diagonal(0,1) - mpower($C/$z->{'lambda'} + $m->diagonal(0,1), $alpha);
            $powr = ( $p->diagonal(0,1) - $svdxs->{'v'} x ($F x transpose($svdxs->{'v'}) )) * ($z->{'lambda'})**$alpha;
        }
        $powr->diagonal(0,1) = 1;
        #rownames(powr) = colnames(xs)
        #colnames(powr) = colnames(xs)
    }
	$powr = {
		'powr' => $powr,
		'lambda' => $z->{lambda},
		'lambda_estimated' => $z->{lambda_estimated},
		'class' => "shrinkage",
	};
    return($powr);
}

sub fast_svd {
	my ($m, $tol) = @_; 
    my $n = $m->getdim(1);
    my $p = $m->getdim(0);
    my $EDGE_RATIO = 2;
    if ($n > $EDGE_RATIO * $p) {
        return(psmall_svd($m, $tol));
    } elsif ($EDGE_RATIO * $n < $p) {
        return(nsmall_svd($m, $tol));
    } else {
        return(positive_svd($m, $tol));
    }
}
=head
sub psmall.svd {
	my ($m, $tol) = @_; 
    my $B = mcrossprod($m);
    my $s = svd($B, 0);
    if (not defined ($tol)) { 
        my $tol = $B->getdim(1) * maximum($s->{'d'}) * 2.220446e-16;
		my $Positive = which($s->{'d'} > $tol);
		my $d = sqrt($s->{'d'}->($Positive));
		my $v = $s->{'v'}->(, $Positive);
		my $u = m %*% v %*% diag(1/d, nrow = length(d))
		return(list(d = d, u = u, v = v))
}
=cut
sub pvt_svar {
	my ($x, $lambda_var, $w, $verbose) = @_;
	$w = pvt_check_w($w,$x->getdim(1));
    my $xs = wt_scale( $x,$w,'TRUE', 'FALSE');
	my $wt_mom = wt_moments($xs,$w);
	my $v = $wt_mom->{var};
	my $tgt = median($v);
    my $z = pvt_get_lambda($xs, $lambda_var, $w, 0, "variance", $tgt);
    my $vs = $z->{'lambda_var'} * $tgt + (1 - $z->{'lambda_var'}) * $v;
	my $back = {
		'vs' => $vs,
		"lambda_var" => $z->{'lambda_var'},
		"lambda_var_estimated" => $z->{'lambda_var_estimated'},
		"class" => "shrinkage",
	};
    return($back);
}

sub pvt_varlambda {
	my ($xc,$w, $target) = @_;
    my $w2 = dsum($w * $w);
    my $h1 = 1/(1 - $w2);
    my $h1w2 = $w2/(1 - $w2);
    my $zz = $xc**2;
    my $q1 = dsumover(  $zz->transpose*$w  );
    my $q2 = dsumover( ($zz->transpose**2)*$w ) - $q1**2;
    my $numerator = dsum($q2);
	my $lambda = undef;
    my $denominator = dsum(($q1 - $target/$h1)**2);
    if ($denominator == 0) {
		$lambda = 1;
	} else {
		$lambda = ($numerator/$denominator) * $h1w2; 
		if ($lambda > 1) { $lambda =1; }
	}
	return($lambda);
}

sub pvt_corlambda {
	my ($xs, $w, $target) = @_;
    my $w2 = dsum($w * $w);
    my $h1w2 = $w2/(1 - $w2);
    my $sw = sqrt($w);
	my $Q1_squared = (mcrossprod($xs*$sw->transpose))**2;
    my $Q2 = mcrossprod(($xs**2)*$sw->transpose)  - $Q1_squared;
    my $denominator = dsum($Q1_squared) - dsum( $Q1_squared->diagonal(0,1));
    my $numerator = dsum($Q2) - dsum( $Q2->diagonal(0,1));
	my $lambda = undef;
	if ($denominator == 0) {
        $lambda = 1;
	} else { 
		$lambda = ($numerator/$denominator) * $h1w2;
		if ($lambda > 1) { $lambda = 1; }
	}
	return($lambda);
}

1;