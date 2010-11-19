#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $VERSION = "0.1";


our ( $help, $man, $out, $fam,$size,$N, $qt );

GetOptions(
    'help|h' => \$help,
    'man' => \$man,
    'out|o=s' => \$out,
    'fam=s' => \$fam,
    'ped=s' => \$fam,    
    'n_subsets|n=i' => \$N, 
    'subset_size|size=f' => \$size,
    'qt' => \$qt,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(0) if (not defined $out);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);

defined $size or $size = 0.1;
defined $N or $N = 10;

my @cases = ();
my @controls = ();

open (FAM,$fam) or die $!;
while (my $line = <FAM>){
    chomp($line);
    my @data = split(/\s+/,$line);
    next if ($data[5] eq 'NA');
    next if ($data[5] == -9);
    if ( $data[5] == 1){ push @controls, "$data[0] $data[1]"; }
    elsif ( $data[5] == 2){ push @cases, "$data[0] $data[1]"; }
    else { push @controls, $data[0]; }
}
close(FAM);


print scalar localtime, "\tThere are [ " . scalar @cases . " ] and [ " . scalar @controls . "] controls\n";
print scalar localtime, "\tWill make [ $N ] subsets of [ " . int(scalar @cases*$size) . " ] cases and [ " . int(scalar @controls*$size) . " ] controls\n";
print scalar localtime, "\tResults will be written to [ $out.subsetNumber ]\n";

my @new_sets = ();
for my $i (1 .. $N){
    my @case_index  = @{get_rand_index(scalar @cases,int(scalar @cases*$size))}; 
    my @control_index  = @{get_rand_index(scalar @controls,int(scalar @controls*$size))}; 
    push @{$new_sets[$i-1]}, (@cases[@case_index], @controls[@control_index]);
    open (OUT,">$out.$i") or die $!;
    print OUT join "\n",(@cases[@case_index],@controls[@control_index]);
    close (OUT);
}
print scalar localtime, "\tWell Done\n";

exit;
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




__END__

=head1 NAME

 Running network analysis by greedy search

=head1 SYNOPSIS

script [options]

 	-h, --help		print help message
 	-m, --man		print complete documentation
        -fam                    PLINK FAM file 
        -ped                    PLINK PED file
        -n_subsets, -n          number of sample subsets to make
        -subset_size, -size     size of each subset
 	-out, -o		Name of the output file
        -qt                     Defines that phenotype is a quantitative trait


=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-fam>

PLINK FAM file 

=item B<-ped>

PLINK PED file

=item B<-n_subsets, -n>

number of sample subsets to make

=item B<-subset_size, -size>

size of each subset

=item B<-out, -o>

Name of the output file

=item B<-qt>

Defines that phenotype is a quantitative trait

=back

=head1 DESCRIPTION

TODO


=cut
