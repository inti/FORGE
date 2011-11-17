#!/usr/bin/perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Getopt::Long;
use IO::File;

our ( $Usage, $help, $out, $chr, $gene_list,$distance,$distance_five_prime,$distance_three_prime);

GetOptions(
   'help|h'   => \$help,
   'out|o=s'   => \$out, #name of the output file
   'chr=i@'  => \$chr,
   'distance|d=i' => \$distance,
   'distance_five_prime|five=i' => \$distance_five_prime,
   'distance_three_prime|three=i' => \$distance_three_prime,
   'gene_list=s' => \$gene_list,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $help);


open (OUT,">$out") or die $!;
my $LOG = new IO::File; 
$LOG->open(">$out.log") or print_OUT("I can not open [ $out.log ] to write to",$LOG) and exit(1);

print_OUT("Opened output and log file [ $out ] and [ $out.log ]",$LOG);

defined $chr or @{ $chr } = (1..22,'X','Y','MT');
print_OUT("Will analyze chromosomes [ @$chr ]");
# start and get registry
print_OUT("Contacting Ensembl DB",$LOG);
my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous', -verbose => '0', );

#start Core adaptors
print_OUT("Generating Data adaptors",$LOG);
my $genes = $reg->get_adaptor("human","core","gene");
my $sliceAdap = $reg->get_adaptor("human","core","slice");

#start Variation adaptors
my $dbVar = $reg->get_DBAdaptor("human","variation");
my $varAdap = $dbVar->get_VariationAdaptor();
my $varFeatAdap = $dbVar->get_VariationFeatureAdaptor();

print_OUT("Fetching Information from Ensembl database version " . $reg->software_version() . "",$LOG);

if (defined $distance) {
	$distance_three_prime = $distance;
	$distance_five_prime = $distance;
    print_OUT("Max SNP-to-gene distance allowed [ $distance ] kb",$LOG);
} elsif ( (defined $distance_three_prime) or (defined $distance_five_prime)){
    defined $distance_three_prime or $distance_three_prime = 50;
    defined $distance_five_prime or $distance_five_prime = 50;
	print_OUT("Max SNP-to-gene distance allowed for 3-prime [ $distance_three_prime ] and 5-prime [ $distance_five_prime ] kb",$LOG);
} else {
	$distance = 50;
	$distance_three_prime = $distance;
	$distance_five_prime = $distance;
	print_OUT("Max SNP-to-gene distance allowed [ $distance ] kb",$LOG);
}

my %selected_genes = ();
if (defined $gene_list){
    open(GL,$gene_list)or die print_OUT("Cannot open [ $gene_list ]");
    my @tmp = <GL>;
    chomp(@tmp);
    map { $selected_genes{$_} = $_; } @tmp;
}

my @annot;
my $size = undef;

foreach my $c (@{ $chr }) {
   if ($c == 23){ $c = 'X'}
   elsif ($c == 24){ $c = 'Y'}
   elsif ($c == 26){ $c = 'MT'}
   elsif ($c == 25){ die("I do know chromosome 25!! I know form 1 to 22 and X (23), Y (24) and MT (26)\nPlease choose one or find another programe :)\n");}
   my $chr_slice = $sliceAdap->fetch_by_region('chromosome',$c);   

   print_OUT("Fetching genes for chromosome $c",$LOG);
    
   my @genes = @ { $genes->fetch_all_by_Slice($chr_slice) };
   print_OUT("   '-> [" . scalar @genes . " ] genes",$LOG);
    my $gene_counter = 0; 
    my $total_genes = scalar @genes;
    if (defined $gene_list){ $total_genes = scalar keys %selected_genes; }
    while (my $g = shift @genes){
        if (defined $gene_list){ 
            next if ((not defined $selected_genes{$g->display_id}) and ( not defined $selected_genes{$g->external_name()}));
            delete $selected_genes{$g->display_id};
        }
        my $right_distance = 0;
        my $left_distance = 0;
        if ($g->strand() < 0){ 
                # gene is on reverse strand: <-----    start <- end, 3' <- 5'
                $right_distance = $distance_five_prime;
                $left_distance = $distance_three_prime;
        } elsif ($g->strand() > 0){
                # gene is on reverse strand: ----->    start -> end, 5' -> 3'
                $right_distance = $distance_three_prime;
                $left_distance = $distance_five_prime;
        }
        my $gene_slice_extended = $sliceAdap->fetch_by_region('chromosome',$c, $g->start - $left_distance,$g->end + $right_distance);
        my $vfs = $varFeatAdap->fetch_all_by_Slice($gene_slice_extended); 
        next if (scalar @{$vfs} == 0);
        print OUT join "\t", (  $g->seq_region_name, # chromosome
                            $g->start, # start
                            $g->end, # end
                            $g->strand, # strand
                            $g->display_id, # ensembl id 
                            $g->external_name(), # hugo name 
                            $g->status, #status  
                            $g->biotype(),# biotype 
        #$g->description() # description
                          );
        foreach my $vf (@{$vfs}){
            print  OUT "\t",$vf->variation_name,":",$vf->start + $g->start - $left_distance - 1,":",$vf->allele_string,":",$vf->strand; 
        }
        print OUT "\n";
        print scalar localtime," \t", progress_bar(++$gene_counter,$total_genes);
        if (defined $gene_list){ 
            if (scalar keys %selected_genes == 0) { 
                print "\n";
                print_OUT("Finished",$LOG);
                exit;
            }
        }
    }
    print "\n";
}

print_OUT("Finished",$LOG);

exit;

# wget-style. routine by tachyon
# at http://tachyon.perlmonk.org/
sub progress_bar {
    my ( $got, $total, $width, $char ) = @_;
    $width ||= 25; $char ||= '=';
    my $num_width = length $total;
    sprintf "|%-${width}s| Done with [ %${num_width}s ] genes of [ %s (%.2f%%) ]\r", 
    $char x (($width-1)*$got/$total). '>', 
    $got, $total, 100*$got/+$total;
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

sub coord_system {
   my $feat = shift;
   if ( $feat->seq_region_name() =~ /[A-Za-z]/ ){ return('undef'); }
   else {return(1);}
}

sub feat_2_chrCoord {
   my $feat =shift;
   if ($feat->coord_system_name ne 'chromosome') { $feat->transform('chromosome'); }
   return($feat);
}

sub dist_2_objects {
   my $gn = shift;
   my $snp = shift;
   if ( $gn->seq_region_end < $snp->seq_region_start){
      return( int(($gn->seq_region_end - $snp->seq_region_start)/1e3) );
   } elsif ( $gn->seq_region_start > $snp->seq_region_start ){
      return( int(($snp->seq_region_start - $gn->seq_region_start)/1e3) );
   } else { return(0); } 
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
 -distance		Distance to 5-prime of the gene (in kb)
 -distance_three_prime	Distance to 3-prime of the gene (in kb)
 -distance_five_prime	Distance to 5-prime of the gene (in kb)
 -set_stat		statistics to calculate over the sub-networks
 -gc_correction		Determine the lamda for the genomic control correction from the input p-values
 -z_score		input values are z_scores (the absolute values will be used)
 -gs_coverage		number [0,1]. Fraction of the gene-set that must be covered by the
 experiment for the gene set to be considered in the analysis
 
 Multivariate Normal Distribution sampling
 -mnd			Estimate significance by sampling from a multivatiate normal distribution
 -mnd_n         Number of MND simulations to calculate gene p-value (default=1000000)
 -mnd_gene_corr Calculate the correlation between gene-statistics from by simulattions 
 
 Output modifiers:
 -append			Append results to output file rather than overwrite it
 -add_file_name		add the input file name to the result
 
 Permutations:
 -perm			number of permutations
 
=head1 OPTIONS
 
=over 8
 
=item B<-help>
 
 Print help message
 
=item B<-man>
 
 print complete documentation
 
=item B<-report>
 
 how often to report advance. Provide an integer X and the program will report adnvance after X networks are analyzed.
  
=back
 
 
 
=cut