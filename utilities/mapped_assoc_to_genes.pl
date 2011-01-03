#!/usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

our ( $help,$man,$list, $dist, $assoc,$keep, $affy_to_rsid, $out);

GetOptions(
   'gene_list=s'    => \$list, 
   'distance=i'   => \$dist,
   'assoc=s@' => \$assoc,
   'affy_to_rsid=s' => \$affy_to_rsid,
   'keep_annot' => \$keep,
   'out|o=s' => \$out,
)or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 1) if (defined $man);
pod2usage(0) if (not defined $assoc);
pod2usage(0) if (scalar @ARGV == 0);

defined $dist or $dist = 20;
print scalar localtime(), "\t", "SNP to gene mapping set to [ $dist ]\n";
my @genes = ();
if (defined $list){
	print scalar localtime(), "\t", "Reading gene list from [ $list ]\n";
	open (LIST,$list) or die $!;
	my @genes = <LIST>;
	chomp(@genes);
	close(LIST);
}
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

print scalar localtime(), "\t", "Reading association files\n";
my %assoc_data = ();
my $info_counter = 0;
foreach my $file (@$assoc){
	chomp($file);
	print scalar localtime(), "\t", "   '-> [ $file ]\n";
	open( ASSOC, $file ) or die ("I can not find [ $file ]\n");
	my $line = 0;
	my %header = ();
	while (  my $a = <ASSOC> ) {
		chomp($a);
		$a =~ s/^\s+//;
		$a =~ s/^\t+//;
		$a =~ s/\s+/ /g;
		$a =~ s/\t+/ /g;
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
		my $info = "";
		foreach my $h (sort keys %header){
			next if ($h eq "SNP");
			next if ($h eq "P");
			$info .= "\t$h:$data[$header{$h}]";
			$info_counter++;
		}
		push @{ $assoc_data{ $data[$header{"SNP"}] }{$file} }, { "p" => $data[$header{"P"}], "info" => $info };
	}
	close(ASSOC);
}
print scalar localtime(), "\t"," [ ", scalar (keys %assoc_data)," ] SNPs with association data\n";
print scalar localtime(), "\t", "Writting output to [ $out ]\n";
open (OUT,">$out") or die $!;
print scalar localtime(), "\t", "Reading SNP to gene mapping files\n";

foreach my $file ( @ARGV ){
	chomp($file);
	print scalar localtime(), "\t", "   '-> [ $file ]\n";
	open( MAP, $file ) or die $!;
	while ( my $read = <MAP> ) {
		chomp($read);
		# the line is separate in gene info and snps. the section are separated by a tab.
		my ($chr,$start,$end,$ensembl,$hugo,$gene_status,$gene_type,$description,@m) = split(/\t+/,$read);
		if ($m[0] !~ m/:{3}/){ $description .= splice(@m,0,1); }
		
		#if (defined $list){
		#	next unless (grep $_ eq $hugo, @genes);
		#	next unless (grep $_ eq $ensembl, @genes);
		#}
		
		# get all mapped snps within the distance threshold,
		my @mapped_snps = ();
		foreach my $s (@m) {
			my ($id,$pos,$allele,$strand) = split(/\:/,$s);
			next if (not defined $id);
			if (( $pos >= $start) and ($pos <= $end)){ 
				push @mapped_snps, { 'id' => $id, 'pos' => $pos}; 
			} elsif ( ( abs ($pos - $start) <= $dist*1_000 ) or ( abs ($pos - $end) <= $dist*1_000 )) { 
				push @mapped_snps, { 'id' => $id, 'pos' => $pos}; 
			}
		}
		
		next if (scalar @mapped_snps == 0);
		# get the gene position info
		#check if gene was in the list of genes i want to analyze
# 		next unless ( ( grep $_ eq $hugo, @genes ) or ( grep $_ eq $ensembl, @genes ) );
		next if (scalar @mapped_snps == 0);
		foreach my $snp (@mapped_snps){ 
			next unless (exists $assoc_data{$snp->{id}});
			#print OUT "$ensembl\t$hugo\t$chr\t$start\t$end\t$snp->{id}\t$snp->{pos}";
			foreach my $stdy (@$assoc){
				if (not exists $assoc_data{$snp->{id}}{$stdy}) { 
					print OUT "$ensembl\t$hugo\t$gene_status\t$gene_type\t$chr\t$start\t$end\t$snp->{id}\t$snp->{pos}\t$stdy\tNA";
					for (my $i = 0; $i < $info_counter; $i++){ print OUT "\tNA";}
					print OUT "\n";
				}
				foreach my $stat ( @{ $assoc_data{$snp->{id}}{$stdy} } ){ 
					print OUT "$ensembl\t$hugo\t$gene_status\t$gene_type\t$chr\t$start\t$end\t$snp->{id}\t$snp->{pos}\t$stdy\t$stat->{p}\t$stat->{info}\n";
				}
			}
		}
	}
	close(MAP);
}
print scalar localtime(), "\t", "Mapeados!!\n";
exit;

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



__END__

=head1 NAME

 Running network analysis by greedy search

=head1 SYNOPSIS

script [options] -- snp-2-gene-mapping-files

	-h, --help		print help message
	-m, --man		print complete documentation
	-gene_list		gene of list to analyze 
   	-distance		max distance between SNPs and genes
   	-assoc			file with SNP p-values. It need a header with at least SNP and P columns
				Multiple file can be provided by giving the option more than once \$assoc,
   	-affy_to_rsid		mapping between affy and rsids file
   	-keep_annot		print in the output file all other columns from the input files
	-out, -o		output file name
        
=head1 OPTIONS

=over 8

=item B<-help>

Print help message
  
=item B<-man>

print complete documentation

=item B<-gene_list>

gene of list to analyze 

=item B<-distance>

max distance between SNPs and genes

=item B<-assoc>

file with SNP p-values. It need a header with at least SNP and P columns. Multiple file can be provided by giving the option more than once.

=item B<-affy_to_rsid>

mapping between affy and rsids file

=item B<-keep_annot>

print in the output file all other columns from the input files

=item B<-out, -o>

output file name

=back

=head1 DESCRIPTION

TODO


=cut