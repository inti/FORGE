#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use PDL;
use PDL::GSL::CDF;
use PDL::Primitive;
use IO::File;
use Data::Dumper;
use Statistics::RankCorrelation;
use PDL::Stats::Basic; 

our ( $help, $man, $out, $snpmap, $bfile, $assoc, $gene_list, @genes,
      $all_genes, $analysis_chr, $report, $spearman, $affy_to_rsid,
      $v,  $lambda,$print_cor, $pearson_genotypes,$distance);

GetOptions(
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
);

# define distance threshold,
defined $distance or $distance = 20;
print scalar localtime(), "\t", "Max SNP-to-gene distance allowed [ $distance ] kb\n";
# $report is use to specify how often report the advance when reading input files
defined $report or $report = 50_000;
# tell if analysis is restricted to a specific chromosome
defined $analysis_chr and print scalar localtime(), "\t", "Restricting analysis to chromosome [ $analysis_chr ]\n";
#  defined lambda value for Genomic Control correction
defined $lambda or $lambda = 1;
# tell if user wants to print the *.correlation file 
defined $print_cor and print scalar localtime(), "\t", "Defined -print_cor: I will print the *.correlation file (it is bulky)\n";
#set output file if not set already.
defined $out or $out = "gene_based_fisher_v041.OUT";
# tell if user wants to correct p-values by genomic control 
if ($lambda != 1){ print scalar localtime(), "\t", "SNP p-value will be corrected with lambda = [ $lambda ]\n";}

# i will read the gene_list and i will load data for just this genes to speed up.
if ( not defined $all_genes ) { # in case user want to analyse all genes
  if ( not defined @genes ) { # in case user gave a list of genes in the command line
    print scalar localtime(), "\t", "Reading Gene List from [ $gene_list ]\n";
    # read file with gene list and store gene names.
    open( GL, $gene_list ) or die $!;
    @genes = <GL>;
    chomp(@genes);
    close(GL);
  } else { 
	print scalar localtime(), "\t", "Read Gene List command line [ @genes  ]\n"; 
  }
} else {
	print scalar localtime(), "\t", "Going to analyze all genes on [ $snpmap ] file.\n"; 
}

# Now lets going to read the affy id to rsid mapping. This is used to keep all ids in the
# same nomenclature
my %affy_id = ();
if ( defined $affy_to_rsid ) { # if conversion file is defined
   print scalar localtime(), "\t", "Reading AFFY to rsID mapping from [ $affy_to_rsid ]\n";
   open( AFFY, $affy_to_rsid ) or die $!;
   while (my $affy = <AFFY>){
      chomp($affy);
      my @b = split(/\t+/,$affy);
      $affy_id{$b[0]} = $b[1];
   }
   close(AFFY);
}

# Read file with genetic association results.
print scalar localtime(), "\t","Reading association file: [ $assoc ]\n";
# creat hash to store SNP information
my %assoc_data = ();
open( ASSOC, $assoc ) or die ("I can not find [ $assoc ]\n");
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
   exists $header{"SNP"} or die("Not p-value columns available, these are the headers i found [ " . (keys %header) . " ]\n");
   exists $header{"P"} or die("Not p-value columns available, these are the headers i found [ " . (keys %header) . " ]\n");
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
   if ($lambda == 1) {
	$assoc_data{ $data[$header{"SNP"}] } = {
                'pvalue'       => $data[$header{"P"}] ,
        };
   } elsif ($lambda > 1) { 
   	$assoc_data{ $data[$header{"SNP"}] } = {
			# transform the p-value on a chi-square, correct it by the infaltion factor and tranfsormit again on a p-value
      		'pvalue'       => 1-  gsl_cdf_chisq_P( gsl_cdf_chisq_Pinv( $data[$header{"P"}], 1 )/$lambda, 1 ),
   	};
   } else { 
	die("\nPlease check the lambda value is correct\n"); 
   }
}
close(ASSOC);

if (scalar keys %assoc_data == 0){ 
	print "\nNo SNPs with genetic assocition to used in the analysis\n\n";
	exit(0);
}
print scalar localtime(), "\t", "[ ", scalar keys %assoc_data ," ] SNPs with association data\n";

#read snp-to-gene mapping and store in a hash with key equal gene name and value
# an array with the snps in the gene.
my @bim = ();
my @fam = ();
my %bim_ids = ();
if (defined $bfile) {
  # read the bim file with snp information and fam file with sample information
  @bim = @{ read_bim("$bfile.bim") };
  map { $bim_ids{$_->{snp_id}} = ""; } @bim;
  @fam = @{ read_fam("$bfile.fam") };
}

print scalar localtime(), "\t", "Reading gene to SNP mapping file from [ $snpmap ]\n";

open( MAP, $snpmap ) or die $!;
my %gene = ();
my %snp_to_gene = ();
while ( my $read = <MAP> ) {
   chomp($read);
   # the line is separate in gene info and snps. the section are separated by a tab.
   my ($chr,$start,$end,$ensembl,$hugo,$gene_status,$gene_type,$description,@m) = split(/\t+/,$read);
   if ($m[0] !~ m/:{3}/){ $description .= splice(@m,0,1); }

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

   # go over mapped snps and change converte affy ids to rsid.
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
	  if (defined $v ){ print "Mapping [ $s ] to [ $ensembl ]\n";}
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
	if (defined $v ){ print "Gene $ensembl $hugo included in the analysis with [ ", scalar @{ $gene{$ensembl}->{snps} }, " ] mapped SNPs.\n"; }
   }
   
}
close(MAP);

print scalar localtime(), "\t", "[ ",scalar keys %gene ," ] Genes will be analyzed \n";
print scalar localtime(), "\t", "[ ",scalar keys %snp_to_gene ," ] SNPs mapped to Genes and with association results will be analyzed \n";


if (scalar keys %gene == 0){
        print "\nNo genes mapped. Quitting here with.\n\n";
        exit(0);
}

# start a hash to store the SNP-to-SNP correlation values
my %correlation = ();

# if provided get the SNP-to-SNP correlation values from a tab separated file with 3 cols:snp1 snp2 correlatio_value
if (defined $spearman){
   print scalar localtime(), "\t", "Reading SNP correlation from [ $spearman ]\n";
   open( SPRMN, $spearman ) or die $!;
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
print scalar localtime(), "\t", "Output file will be written to [ $out ]\n";
my $out_fh = new IO::File;
$out_fh->open(">$out");
print $out_fh  "Ensembl_ID\tHugo_id\tgene_type\tchromosome\tstart(ensemblv51+20kb)\tend(ensemblv51+20kb)\tmin_p\tmin_p_sidak_corrected_by_number_of_tests\tfisher_combined_gene-pvalue\tchi-square\tdegrees_of_freedom\tGalwey_pvalue\tGalwey_chi_square\tGalwey_degrees_freedom\tnumber_of_snps\tnumber_effective_tests(Gaoetal;PDMID:19434714)\n";
#start count to report advance 
my $count = 0;
# if user defined a genotypes file. read genotypes and store a genotyop matrix (rows: samples, cols: genotypes)for each gene 
if (defined $bfile) {
	  # open binary genotype file
	  my $bed = new IO::File;
	  $bed->open("<$bfile.bed");
	  binmode($bed);
	  # start a line counter to evaluate when we are in the first line. the first line has information
	  # on the structure of the file
	  $line         = 0;
	  # start a snp and sample index. the genotype file does not have information about samples or snps. So in order
	  #  to identify them we need to match them with information from the bin and fam files.
	  my $index_sample = 0;
	  my $index_snp    = 0;
	  my @tmp_genos = ();
	  # put the whole binary file in memory to speed up things
	  print scalar localtime(), "\t","Reading genotypes from [ $bfile.bed ]\n";
	  # go over each binary information line
	  foreach my $d (<$bed>) {
		 # every line is made up of 8 blocks of information. so specify the length of the whole set of blocks  	
		 my $l = ( length $d ) * 8;
		 # now extract the amount of information corresponding to the 8 blocks from $d
		 foreach  my $data ( unpack( "B$l", $d ) ) {
			# now put the 8 blocks or bytes on an array to process them one by one.
			my @bytes = ( $data =~ m/\d{8}/g );
			# if it is the first line we need to get the information regarding the structure of the file
			# this information is on the first 3 bytes. the first two say that this is a PLINK PED file in binary format
			# the 3rd says wheter it is on SNP or individual major format. see PLINK web page for more details http://pngu.mgh.harvard.edu/purcell/plink/ 
			if ( $line == 0 ) {
			   if ( ( $bytes[0] eq '01101100' ) and ( $bytes[1] eq '00011011' ) ) { # check wheter it is a PLINK binary file
				  if ( $bytes[2] eq '00000001' ) { # check wether it is on SNP major format
					 print scalar localtime(), "\t","Binary file is on SNP-major format\n";
				  }
				  elsif ( $bytes[2] eq '00000000' ) { # check wether it is on individual major format
					 print scalar localtime(), "\t","Binary file is on Individual-major format\n";
					 die("This format is not suppported at the moment\n");
				  } else { die("I can not recognize if the file has SNP-major or Individual-major format\n"); }
				  # increase line count and delete the first 3 bytes, next ones will have genotype information
				  $line++;
				  splice( @bytes, 0, 3 );
			   }
			   else { die("This does not seem to be a PLINK binary file\nBye for now\n\n"); }
			}
			# Now read the actual genotypes and generate the genotype matrix for each gene.
			# in the SNP major format (the only supported at the moment) the information is
			# all genotypes for the SNP in the order of the samples extracted from the fam file.
			# SNPs will come in the order of the bim file
			foreach  my $set (@bytes) {
			   # each byte will have 8 numbers which correspond to 4 genotypes  
			   my @g = ( $set =~ m/\d\d/g );
			   # the genotypes are reversed so I need to reverse them to read the wrigth genotype
			   while ( my $geno = reverse pop(@g) ) {
				  # in case this SNP is mapped to an SNP i will process it. if not just skip it
				  if (exists $snp_to_gene{${ $bim[$index_snp] }{snp_id} }) {
					 if    ( $geno eq '00' ) {  # homozygote
						push @tmp_genos, 1;
					 } elsif ( $geno eq '11' ) { # -- other homozygote
						push @tmp_genos, 3;
					 } elsif ( $geno eq '01' ) { # -- heterozygote
						push @tmp_genos, 2;
					 } elsif ( $geno eq '10' ) { # -- missing genotype
						push @tmp_genos, 0;
					 } else { print scalar localtime(), "\t","This genotype is not recognize [ $geno ]\n"; }    # genotype not recogniize
					 # increase the sample count after reading the genotype
					 $index_sample++;
					 # if reached the last sample it is moment to store all genotypes for this SNP.
					 if ( $index_sample == scalar @fam ) {
						# store SNP genotypes only if it has association data
						if (exists $assoc_data{ ${ $bim[$index_snp] }{snp_id} } ){
						   # for every gene mapped to this SNP, push inside the gentype matrix this SNP genotypes.
						   foreach my $gn (@{ $snp_to_gene{${ $bim[$index_snp] }{snp_id} } }){
							  if (defined $v){ print "Adding SNP [  ${ $bim[$index_snp] }{snp_id}  ] to genotypes of $gn\n"; }
							  next if ( grep $_ eq ${ $bim[$index_snp] }{snp_id} , @{ $gene{$gn}->{geno_mat_rows} } );
							  # make a piddle on PDL with the genotypes
							  my $tmp_pdl = pdl @tmp_genos;
							  # add the genotypes to the genotype matrix
							  $gene{$gn}->{genotypes} = $gene{$gn}->{genotypes}->glue(1,$tmp_pdl);
							  push @{ $gene{$gn}->{geno_mat_rows} }, ${ $bim[$index_snp] }{snp_id};
							  push @{ $gene{$gn}->{pvalues} }, $assoc_data{ ${ $bim[$index_snp] }{snp_id} }->{pvalue};
							  if (scalar @{ $gene{$gn}->{snps} } == scalar @{ $gene{$gn}->{geno_mat_rows} }){
								&gene_pvalue($out_fh,$gn);
								delete($gene{$gn});
								$count++;# if there are more than 100 genes change the $report variable in order to report every ~ 10 % of genes.
								unless (scalar keys %gene < 100){
								  my $report = int((scalar keys %gene)/100 + 0.5)*10;
								  &report_advance($count,$report,"GENES");
								}
							  }
						   }
						}
						# increase SNP count and set sample id to 0
						$index_snp++;
						$index_sample = 0;
						# empty the genotype information
						@tmp_genos = ();
						# report the advance on reading the SNPs
						#&report_advance($index_snp,$report,"SNPs");
						# regarding whether the SNP genotype finish at the middle of the byte
						# the next SNP starts at the new byte, so jump to the next byte
						last;
					 }
				  } else { # if the SNP is not map to gene, increase sample count
					 $index_sample++;
					 if ( $index_sample == scalar @fam ) {
						# increase SNP count and set sample id to 0
						$index_snp++;
						$index_sample = 0;
						# empty the genotype information
						@tmp_genos = ();
						#&report_advance($index_snp,$report,"SNPs");
						# regarding whether the SNP genotype finish at the middle of the byte
						# the next SNP starts at the new byte, so jump to the next byte
						last;
					 }
				  }
			   }
			}
		 }
	  }
} else {
	  print scalar localtime(), "\t", "WARNING: Gene p-values will be calculated with the precomputed correlation only. If correlation for some SNPs pairs are missing you may get wrong results, please check your inputs for completeness.\n";
}


# loop over all genes and calculate the gene p-values
print scalar localtime(), "\t","Starting to Calculate gene p-values\n";
# if there are more than 100 genes change the $report variable in order to report every ~ 10 % of genes.
unless (scalar keys %gene < 100){
	  $report = int((scalar keys %gene)/100 + 0.5)*10;
}



foreach my $gn (keys %gene){
  &gene_pvalue($out_fh,$gn);
  delete($gene{$gn});
  $count++;
  &report_advance($count,$report,"genes");	
}

$out_fh->close();

# if the user want to get the correlation values print the *.correlation file
if (defined $print_cor){
	  open (COR,">$out.correlation") or die$!;
	  foreach my $snp1 (keys %correlation) {
		 foreach my $snp2 (keys %{$correlation{$snp1}}  ) {
			next if ($snp1 eq $snp2);
			printf COR ("$snp1 $snp2 %.3f\n",abs($correlation{ $snp1 }{ $snp2 }));
			delete($correlation{ $snp1 }{ $snp2 });
		 }
	  }
	  close(COR);
}
print scalar localtime(), "\t","Well Done!!\n";
exit;
sub gene_pvalue {
	my $out_fh = shift;
	my $gn = shift;
   if (defined $v){ print "____ $gn ____\n"; }
   if (not defined $bfile){ @{ $gene{$gn}->{geno_mat_rows}} = @{ $gene{$gn}->{snps}}; }
   my $n_snps = scalar @{ $gene{$gn}->{geno_mat_rows}};
   next if ($n_snps == 0);
   # if the gene has just 1 SNP we make that SNP's pvalue the gene p-value under all methods
	  if ($n_snps == 1){
			if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\tNA\tNA\t%0.3e\tNA\tNA\t1\t1\n",$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp}); }
			printf $out_fh ("$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\tNA\tNA\t%0.3e\tNA\tNA\t1\t1\n",$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp},$gene{$gn}->{minp});
			delete($gene{$gn});
			return();
	  }
   #if the user defined a genotype file then we need to get the number of SNPs in the gene.
   # if the user did not defined neither a genotype nor a file with the SNP-SNP correlations exit the program
   if (not defined $bfile and not defined $spearman) { die("You MUST specify file with the SNP-SNP correlations or a genotype file to calculate them by myself\n\n"); } 
   
   if (defined $v){ print scalar localtime(), "\t", "Calculating correlation matrix for $n_snps SNPs with association data\n"; }
   # initialize correlation matrix with zeroes
   my $cor = zeroes $n_snps,$n_snps;
   # set equal weight to all SNPs
   my $weight = 1/($n_snps);
   if (defined $v){ print "Weigth = [ $weight ]\n"; }
   # initialize the second term for the variance of the statistics  for the Makambi method
   my $second_term = 0;
   # loop over all pairs of SNPs and get their correlation values
   for ( my $i = 0; $i < $n_snps; $i++){
      for ( my $p = $i; $p < $n_snps; $p++){
	 # if the correlation is known
	 if (defined $v){ print "working in $i:$gene{$gn}->{geno_mat_rows}->[$i]; $p:$gene{$gn}->{geno_mat_rows}->[$p]\n"; }
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
		$second_term += $weight*$weight*(3.25*abs($c) + 0.75*abs($c)**2);
      } else { # if correlation is unknown calculate the LD (r) between the snps
		# if the user did not provided a genotype file and tell that this pair of SNPs does not have a correlaiton value
		if (not defined $bfile){print " WARNING: this pair of snps does not have a correlation value: [$gn => $gene{$gn}->{snps}->[$p] ] and [ $gene{$gn}->{snps}->[$i] ]\n"; }
		# compute the correlation
		my $c = "";
		if (defined $pearson_genotypes) { # if requested used the pearson correlation for genotypes
		  $c = pearson_corr_genotypes([list $gene{$gn}->{genotypes}->slice(",($i)")],[list $gene{$gn}->{genotypes}->slice(",($p)")]);
		} else { # by default use the spearman rank correlation
		  my $d = Statistics::RankCorrelation->new([list $gene{$gn}->{genotypes}->slice(",($i)")],[list $gene{$gn}->{genotypes}->slice(",($p)")]);
		  $c = $d->spearman;
		}
		set $cor, $i,$p, $c;
		set $cor, $p,$i, $c;
		# increase the second term of the Makambi method
		$second_term += ($weight**2)*(3.25*(abs($c)) + 0.75*(abs($c)**2));
		# store the correlation value in case is needed later.
		$correlation{ ${ $gene{$gn}->{geno_mat_rows} }[$i] }{${ $gene{$gn}->{geno_mat_rows} }[$p]} = $c;
		$correlation{ ${ $gene{$gn}->{geno_mat_rows} }[$p] }{${ $gene{$gn}->{geno_mat_rows} }[$i]} = $c;
      }
    }
  }

   if (defined $v){ print $cor,"\n";}
   if (defined $v){ print scalar localtime(), "\t","Calculating effetive number of tests: "; }
   # calculate number of effective tests by the Gao ($k) and Galwey ($Meff_galwey) method.	  
   my ($k,$Meff_galwey) = number_effective_tests(\$cor);
   # calculate first term of the variance, the variance and the degrees of freedom of the statatistics for the Makambi method
   my $first_term = 4*($n_snps*$weight**2);
   my $variance = $first_term + $second_term;
   my $degrees_freedom = 8/$variance;
   
   # get the log of the SNP p-value
   $gene{$gn}->{pvalues} = pdl @{ $gene{$gn}->{pvalues} };
   # calculate the chi-square statistic and the p-value for the Galwey method
   my $chi_stat_galwey = -2 * ($Meff_galwey/$n_snps) * $gene{$gn}->{pvalues}->log->sum;
   my $fisher_p_value_galwey = 1 - gsl_cdf_chisq_P( $chi_stat_galwey, 2*$Meff_galwey );
   
   # calculate the chi-square statistics for the Makambi method and its p-value
   my $chi_stat = sumover(-2 * $gene{$gn}->{pvalues}->log * $weight);
   $chi_stat = ( $chi_stat/2 ) * $degrees_freedom;
   my $fisher_p_value =   1 - gsl_cdf_chisq_P($chi_stat, $degrees_freedom );

   # print out the results
   if (defined $v){ printf (scalar localtime() . "\t$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t%0.3e\t%0.5f\t%0.5f\t%2d\t%3d || $weight\t$gene{$gn}->{pvalues}->log\t@{ $gene{$gn}->{geno_mat_rows} }\n",$gene{$gn}->{minp},1-(1-$gene{$gn}->{minp})**$k,$fisher_p_value,$chi_stat,$degrees_freedom,$fisher_p_value_galwey,$chi_stat_galwey,$Meff_galwey,$n_snps,$k); }
   printf $out_fh ("$gn\t$gene{$gn}->{hugo}\t$gene{$gn}->{gene_type}\t$gene{$gn}->{chr}\t$gene{$gn}->{start}\t$gene{$gn}->{end}\t%0.3e\t%0.3e\t%0.3e\t%0.5f\t%0.5f\t%0.3e\t%0.5f\t%0.5f\t%2d\t%3d\n",$gene{$gn}->{minp},1-(1-$gene{$gn}->{minp})**$k,$fisher_p_value,$chi_stat,$degrees_freedom,$fisher_p_value_galwey,$chi_stat_galwey,$Meff_galwey,$n_snps,$k);
}

# this subrutine take an array and return a hash were every element of the line is
# a key and the value is the index in the array
sub get_header {
   my $in = shift;
   my %back = ();
   for (my $i = 0;$i< scalar @$in; $i++){
      $back{$$in[$i]} = $i;
   }
   return(\%back);
}

# this subrutine calculate the number of effective test by the Gawey and Gao method.
sub number_effective_tests {
   my $mat = shift;
   # calculate the eigen value of the correlation matrix
   my $eigens = eigens ${$mat};
   # normalize teh values
   my $eigens_norm =  pdl sort { $b <=> $a } list ($eigens/(sumover $eigens) );
   # calculate number of offective test per Gao et al Gao, X., Becker, L. C., Becker, D. M., Starmer, J. D. & Province, M. A. Avoiding the high Bonferroni penalty in genome-wide association studies. Genet. Epidemiol. (2009).
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
   if (defined $v){ print " simpleM_Gao = $simpleM; Meff_galwey=$Meff_galwey $numerator $denominator\n"; }
   return($simpleM,$Meff_galwey);
}

# this subrutine recives a list of SNPs and return the log of the their p-values
sub get_logs {
   my $snps = shift;
   my @stats = ();
   if (defined $v) { print "____ ASSOCIATIONS ____ \n  [ ",scalar @$snps," ] SNPs \n"; }
   foreach my $s (@{$snps}){
      my $p = $assoc_data{$s}->{pvalue};
      next if ($p eq "NA");
      if ($p == 0){ push @stats ,1;}
      else {
         push @stats ,log($p);
      }
      if (defined $v) {  print "\t$s\t$p\n"; }
   }
   if (defined $v) { print "______________________\n"; }
   return(\@stats);
}

# this subrutine print the adnvance of of something
# it receives 3 parameter: index, rep and tag. index is the adnvance so far, report how often to report and tag is name for the printing
# it will print a report when index if divisible by rep
sub report_advance {
	  my ($index,$rep,$tag) = @_;
   if (( $index/$rep - int($index/$rep)) == 0) {print scalar localtime(), "\t", " '->Done with [ $index ] $tag\n"}
}

# this subrutine read the fam file and stores the informaiton in an array
# the elements of the array are pseudohashes with all the sample's information
sub read_fam {
   my $fam = shift;
   print scalar localtime(), "\t", "Reading samples info from [ $fam ]\n";
   open( FAM, $fam ) or die $!;
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
   print scalar localtime(), "\t", "[ ", scalar @back, " ] samples read\n";
   return ( \@back );
}

# this subrutine read the bim file and store information about the SNPs
# each element of the array returned is a pseudohash with the SNP information
sub read_bim {
   my $bim = shift;
   print scalar localtime(), "\t","Reading SNPs info from [ $bim ]\n";
   open( BIM, $bim ) or die $!;
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
   print scalar localtime(), "\t", "[ ", scalar @back, " ] SNPs on BED file\n";
   return ( \@back );
}


sub pearson_corr_genotypes {
  # implemented as in S. Wellek, A. Ziegler, Hum Hered 67, 128 (2009).
  # genotypes must be coded as 1,2 and 3,any other coding will be use. missing genotypes can be set to anything different of 1, 2 or 3.
  # this method will fail is more than 2 out of the four homozygoues haplotypes have counts 0. In such case i recomend to use the default option of the spearman correlation. 
  my $snp1 = shift;
  my $snp2 = shift;
  # initialize to 0 all values for the hapltype cells
  my @cells = (0,0,0,0,0,0,0,0,0);
  # count how many oberservation of every haplotype are
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
	print "  AA Aa aa\n";
	print  "BB $cells[0]  $cells[1]  $cells[2]\n";
	print  "Bb $cells[3]  $cells[4]  $cells[5]\n";
	print  "bb $cells[6]  $cells[7]  $cells[8]\n";
  }
  my $h1 = sqrt($cells[0]);
  my $h2 = sqrt($cells[2]);
  my $h3 = sqrt($cells[6]);
  my $h4 = sqrt($cells[8]);
  my $corr = 2*($h1 * $h4 - $h2 * $h3)/sqrt(4*($h1 + $h2)*($h3 + $h4)*($h1 + $h3)*($h2 + $h4));
  return($corr);
}