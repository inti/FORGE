#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;

my $VERSION = "0.9.5.6";


our ( $help, $man, $out, $files, $cols,$outfile,$keep,$space);

GetOptions(
   'help|h' => \$help,
   'man' => \$man,
   'file|f=s@' => \$files,
   'column|c=i@' => \$cols,
   'out|o=s' => \$outfile,
   'keep' => \$keep,
   'space' => \$space,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);
pod2usage(0) if (not defined $files);
pod2usage(0) if (not defined $outfile);

my $key_number = scalar @$cols/scalar @$files;

# make col number run from 0
for (my $i = 0; $i < scalar @$cols; $i++){ $cols->[$i]--; }

print_OUT("You provided [ $key_number ] keys per file");
print_OUT("Reading input files");
my %key_data = ();
foreach my $f (@$files){
    print_OUT("   '-> Reading [ $f ]");
    my @cols_to_extract = splice(@$cols,0,$key_number);
    print_OUT("   '-> merging by columns [ @cols_to_extract ]");
    open (IN,$f) or print_OUT("Cannot open [$f] ") and exit(1);    
    while (my $line = <IN>) {
        next if ($line =~ m/#/);
        chomp($line);
        my @data= split(/[\s+\t+]/,$line);
        my $key = join "", @data[@cols_to_extract];
        unless ($f eq $files->[0]){
            if (not defined $keep){
                if (scalar @cols_to_extract == 1){
                    splice(@data,@cols_to_extract,1);
                } else {
                    for (my $i = 0; $i < @cols_to_extract; $i++ ) {
                        splice(@data,$cols_to_extract[$i] - $i,1);
                    }
                }
            }
        }
        push @{ $key_data{$key}->{line} }, [@data];
        $key_data{$key}->{count}++;
        if ($f eq $files->[0]) { $key_data{$key}->{index} = $.; }
    }
    close(IN);
}

print_OUT("Writting to [ $outfile ]");
open (OUT,">$outfile") or print_OUT("Cannot open [ $outfile ] ") and exit(1);
my @lines_out = ();
foreach my $k (sort { $key_data{$a}->{index} <=> $key_data{$b}->{index} } keys %key_data){
    next if ($key_data{$k}->{count} < scalar @$files);
    my $delim = "\t";
    if (defined $space){ $delim = " "; }
    my $ol = join $delim, map {  @{$_}  } @{ $key_data{$k}->{line} }[0..(scalar @{ $key_data{$k}->{line} }) - 1];
    $ol .= "\n";
    push @lines_out,$ol;
    
    if (scalar @lines_out > 500){
        print OUT @lines_out;
        @lines_out = ();
    }
    
}
print OUT @lines_out;

print_OUT("Files Merged"); 
exit;

sub print_OUT {
  my $string = shift;
  print scalar localtime(), "\t$string\n";
}



__END__

=head1 NAME

to use do something like
$ perl row_merge.pl -f file1 -c 1 -c 2 -c 6 -f file2 -c 2 -c 4 -c 6 -o file1_file2_merged
where the columns specified after a file tell the program which fields to use to merger the files. In this way you can specify diff columns for different files


=head1 SYNOPSIS

script [options]

 	-h, --help		print help message
 	-m, --man		print complete documentation
        -file, -f               Input files to be merge
        -column, -c             columns by which to merge the file
        -out, -o                output file
        -keep                   keep the columns by which the files were merge. Default is to keep them only in the first file
        -space                  Output file separated by spaces. Default is tab.

=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-file, -f>

Input files to be merge

=item B<-column, -c>

Columns by which to merge the file

=item B<-out, -o>

Output file

=item B<-keep>

Keep the columns by which the files were merge. Default is to keep them only in the first file

=item B<-space>

Output file separated by spaces. Default is tab.



=back

=head1 DESCRIPTION

intipedroso@gmail.com

=cut

