#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our ( $help, $man, $sif, $min_size,
    $extract,$gene_list, $max_size,
    $outfile, $reduce, %net,$max_path,
    $overlap,
);


GetOptions(
	'help' => \$help,
	'man' => \$man,
    'sif=s' => \$sif,
    'extract_sif' => \$extract,
    'gene_list=s' => \$gene_list,
    'out=s' => \$outfile,
    'reduce_network' => \$reduce,
    'min_size=i' => \$min_size,
    'max_size=i' => \$max_size,
    'max_path' => \$max_path,
    'overlap' => \$overlap,
    'gmt=i' => \$gmt,
) or pod2usage(0);

my $VERSION = "0.2";
print scalar localtime(), "\t", "Running version  [ $VERSION ]\n";	

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);
if (not defined $outfile){
    print scalar localtime(), "\t",  "Please define the output file\n";
    pod2usage(0); 
}
if (not defined $sif){
    print scalar localtime(), "\t", "Please define the network file (in SIF format)\n";
    pod2usage(0); 
}
defined $max_size or $max_size = 1000;
defined $min_size or $min_size = 2;
defined $max_path or $max_path = 2;


if (defined $reduce){
    # read netowrk data and fill the %net hash with it
    &read_sif($sif);
    open (OUT,">$outfile") or die $!;
    print scalar localtime(), "\t", "Printing out to [ $outfile ]\n";	
    my @out_lines = ();
    my %out_nets = ();
    foreach my $node (keys %net) {
	# skip this node if it has more partners than the max subnetwork size permited
	next if (scalar (keys %{ $net{$node}->{'partners'} }) > $max_size);
	my %world = ();
	$world{$node} = 0;
	# make sub-network with all nodes around the seed node up to a max distance equal to $max_path
	# add each of these nodes to %world, so we know which ones we have seen already.
	for (my $i = 0; $i < $max_path; $i++){
            # loop over nodes so far.
            foreach my $node_so_far (keys %world){
                    foreach my $partner (keys %{ $net{$node_so_far}->{'partners'} } ){
                            next if (exists $world{$partner});
                            $world{$partner} = $i;      
                    }
            }
            # discard this sub-network if it is larger than the defined limits.
            last if ( scalar keys %world > $max_size);
	}
        # discard this sub-network if it is larger or smaller than the defined limits.
        next if ( scalar keys %world > $max_size);
        next if ( scalar keys %world < $min_size);
        my $key = join "", keys %world;
        $out_nets{$key} = [keys %world]; 
    }
    my $counter = 1;
    foreach my $net (keys %out_nets){
        print OUT "subnet_$counter\tsubnet_$counter\t";
        print OUT join "\t", sort {$a cmp $b } @{ $out_nets{$net} };
        print OUT "\n";
        $counter++;
    }
}

if (defined $extract) {
    open (GENES,$gene_list) or die $!;
    my %genes = ();
    while (<GENES>){
            chomp($_);
            $genes{lc($_)}= "";
    }
    close(GENES);

    open (OUT,">$outfile") or die $!;
    print scalar localtime(), "\t", "Printing out to [ $outfile ]\n";	
    open (SIF,$sif) or die $!;
    while (<SIF>){
            chomp($_);
            my @int = split(/\t/,$_);
            next unless ((exists $genes{lc($int[0])}) and (exists $genes{lc($int[2])}));
            print OUT "@int\n";
    }
}

if (defined $overlap){
    my %genes = ();
    my @pathways = ();
    my @p_genes = ();
    open( GMT, $gmt ) or print "Cannot open [ $gmt ]" and die $!;
	while ( my $path = <GMT> ) {
		my ( $p_name, $p_desc, @p_genes ) = split( /\t/, $path );
		foreach my $gn ( @p_genes ) {
			$gn = lc($gn);
			if ( $gn =~ m/\// ) {
				$gn =~ s/\s+//g;
				my @genes = split( /\/{1,}/, $gn );
				map {
					$total++;
					push @p_genes, $_;
                    $genes{ $_ } = $p_name;
				} @genes;
			} else {
				$total++;
				push push @p_genes, $_;
                $genes{ $_ } = $p_name;
			}
		}
		my $p = { 'name' => $p_name,
            'desc' => $p_desc ,
            'N_all' => $total,
            'genes' => \@p_genes,
        };
		
		# this section reduces the gene sets to sets of genes of non-overlapping recombination intervals
		# it first check in which recombination intervals the genes are
		# then it will keep those that do not share recomb inter with other genes.
		# the cases where more than one gene belong to the same recombination interval
		# are solve by choosing the best p-value per recombination interval.
		push @pathways, $p;
	}
	close(GMT);
    

}
print scalar localtime(), "\t", "Done!!\n";

exit;

sub read_sif {
	my $file = shift;
	print scalar localtime(), "\t", "Reading network from [ $file ]\n";
	open (IN,$file) or die $!;
	while (my $l = <IN>){
		chomp($l);
		my @data = split(/\t/,$l);

		foreach my $n (@data[(0,2)]){ new_node($n);  }; 
		fill_edge($data[0],$data[2],$data[1]);		
	}
	close(IN);
	
}

sub fill_edge {
        my ($id1,$id2,$edge) = @_;
        if (not exists $net{$id1}->{'partners'}->{$id2} ){
                $net{$id1}->{'partners'}->{$id2} = {
                                                     'count' => 0,
                                                     'edge_type' => [],
                                                    };
        }
        if (not exists $net{$id2}->{'partners'}->{$id1} ){
                $net{$id2}->{'partners'}->{$id1} = {
                                                     'count' => 0,
                                                     'edge_type' => [],
                                                    };
        }
        if ($id1 eq $id2){
                $net{$id1}->{'partners'}->{$id2}->{'count'}++;
                push @{  $net{$id1}->{'partners'}->{$id2}->{'edge_type'}  }, $edge;
                $net{$id1}->{'self'}++;
        } else {
                $net{$id1}->{'partners'}->{$id2}->{'count'}++;
                push @{  $net{$id1}->{'partners'}->{$id2}->{'edge_type'}  }, $edge;
		$net{$id1}->{'degree'} = scalar keys %{ $net{$id1}->{'partners'} };
                
		$net{$id2}->{'partners'}->{$id1}->{'count'}++;
                push @{  $net{$id2}->{'partners'}->{$id1}->{'edge_type'}  }, $edge;
        	$net{$id2}->{'degree'} = scalar keys % { $net{$id2}->{'partners'} };
	}
}

sub new_node {
	my $id = shift;
	if (not exists $net{$id}){
        	$net{$id} = {   
                		'partners' => {},
                                'degree' => 0,
                                'degree_degree_corr' => 0,
				'self' => 0,
                                };
        }
}

__END__

=head1 NAME

 Perl implementation of the PAGE: parametric analysis of gene set enrichment. As bonus includes correction for gene clusters.

=head1 SYNOPSIS

script [options]

	General options
 	-h, --help		print help message
 	-m, --man		print complete documentation
        -sif                    Network file in SIF format
        -out                    output file
        
        Extract subnetwork mode
        -extract_sif            set the program in extract mode
        -gene_list              file with gene ids. All interaction between genes in the list will be extracted
        
        Reduce Network mode     produce a set of gene-sets based on the network
        -reduce_network         set program in reduce mode
        -min_size               min size for the gene-set extracted  
        -max_size               max size for the gene-set extracted
        -max_path               max path length from the seed node to define a gene-set

=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-sif>

Network file in SIF format

=item B<-out>

Output file
        
=item B<-extract_sif>

Set the program in extract mode
        
=item B<-gene_list>

File with gene ids. All interaction between genes in the list will be extracted
        
=item B<-reduce_network>

Set program in reduce mode

=item B<-min_size>

Min size for the gene-set extracted  

=item B<-max_size>

Max size for the gene-set extracted

=item B<-max_path>

Max path length from the seed node to define a gene-set

=back

=head1 DESCRIPTION

B<This program> will read a file with a network in SIF format and perform some useful things, like extracting a subnetwork for a set of genes or generating gene sets based on the network tructure. Contacme at intipedroso@gmail.com for more details.


=cut
