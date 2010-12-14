#!/usr/bin/perl -w
use strict;
use IO::File;
use Pod::Usage;

pod2usage(0) if ($ARGV[0] =~ m/h/);
pod2usage(0) if ($ARGV[0] =~ m/help/);

my $file = IO::File->new();
my $index = IO::File->new();
$file->open("<$ARGV[0]") or print("I can not open [ $ARGV[0] ]") and exit(1);
$index->open("+>$ARGV[0].idx") or print("I can not open [ $ARGV[0] ]") and exit(1);
build_index(*$file, *$index);
$file->close;
$index->close;

exit;


sub build_index {
    my $data_file  = shift;
    my $index_file = shift;
    my $offset     = 0;
	
    while (<$data_file>) {
        print $index_file pack("N", $offset);
        $offset = tell($data_file);
    }
}

__END__

=head1 NAME
 
 Generate binary index file.
 
=head1 SYNOPSIS
 
 perl make_binary_index.pl filename
 
 Output file will be filename.idx
 
 
=cut

