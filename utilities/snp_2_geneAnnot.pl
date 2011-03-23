#!/usr/bin/perl -w
use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Getopt::Long;

our ( $Usage, $help, $out, $chr, $shift);

GetOptions(
   'help|h'    => \$help,
   'out|o=s'   => \$out, #name of the output file
   'chr=i@'  => \$chr,
   'window_size|w=i' => \$shift,
);

my $VERSION = "0.1";
print "SNP-to-gene annotation script version [ $VERSION ]\n";
defined $chr or @{ $chr } = (1..22,'X','Y','MT');
defined $out or die("please give the name of the output file (you forgot option -o !!!!)\n");
open (OUT,">$out") or die $!;

# start and get registry
my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous', -verbose => '0', );

#start Core adaptors
my $genes = $reg->get_adaptor("human","core","gene");
my $sliceAdap = $reg->get_adaptor("human","core","slice");

#start Variation adaptors
my $dbVar = $reg->get_DBAdaptor("human","variation");
my $varAdap = $dbVar->get_VariationAdaptor();
my $varFeatAdap = $dbVar->get_VariationFeatureAdaptor();

# start population adaptor
my $pop_name = 'CSHL-HAPMAP:HapMap-CEU'; #we only want LD in this population
my $pop_adaptor = $dbVar->get_PopulationAdaptor; #get adaptor for Population object
my $pop = $pop_adaptor->fetch_by_name($pop_name); #get population object from database
my $ldFeatContAdap = $dbVar->get_LDFeatureContainerAdaptor; #get adaptor for LDFeatureContainer object




warn "Fetching Information from Ensembl Version ", $reg->software_version(),"\n";

my @annot;
my $size = undef;

foreach my $c (@{ $chr }) {
   if ($c == 23){ $c = 'X'}
   elsif ($c == 24){ $c = 'Y'}
   elsif ($c == 26){ $c = 'MT'}
   elsif ($c == 25){ die("I do know chromosome 25!! I know form 1 to 22 and X (23), Y (24) and MT (26)\nPlease choose one or find another programe :)\n");}
   my $gene_slice = $sliceAdap->fetch_by_region('chromosome',$c);   
   print"Fetching genes for chromosome $c\n";
   my @genes = @ { $genes->fetch_all_by_Slice($gene_slice) };
   print "\t", scalar @genes,"\n";
   !defined $size and $size = 999999999999;
   for (my $pos = 0; $pos <= $size + $shift; $pos += $shift){
      
      my $slice = $sliceAdap->fetch_by_region('chromosome',$c,$pos,$pos + $shift);
      if (!defined $size) {
         print "Defining chr$c size to ", $slice->seq_region_length(),"\n"; 
         $size = $slice->seq_region_length();
      }
      print "Fetching variations for chromosome $c from ",$pos," ",$pos + $shift,"\n";
      my $ldFeatCont = $ldFeatContAdap->fetch_by_Slice($slice,$pop,$pos + $shift);
      my @vars;
      map { push @vars, @ { $varFeatAdap->fetch_all_by_Variation($_) } } @{ $ldFeatCont->get_variations() } ;
      print "\t", scalar @vars,"\n";
      next if (scalar @vars == 0);
      my $var_index = 0;
      my $anchor = 0;
      foreach my $v (@vars){
         next unless ($v->var_class() eq 'snp');
         $v = feat_2_chrCoord($v);
         my $state = coord_system($v);
         next if ($state eq 'undef');
         #print $v->display_id," ", $v->seq_region_name," ",$v->seq_region_start," ",$v->seq_region_end," ",$v->var_class(),"\n";
         unless ($var_index == 0) {
            my $distance_vars = dist_2_objects($vars[$anchor],$v);
            if ($distance_vars == 0) {
               $annot[$var_index + 1] = $annot[$var_index];
               $var_index++;
               next;
            } elsif (abs($distance_vars) > 2e6) { last; }
         }
         
         for (my $i = 0; $i < scalar @genes;$i++){
            my $distance = dist_2_objects($genes[$i],$v);
            if (abs($distance) <= 1e3) {
               push @{ $annot[$var_index] } , {'variation'=>$v, 'gene'=> $genes[$i], 'distance'=>$distance, 'index' => $i};
            }
         }
         $var_index++;
         $anchor = $var_index;
      }
      foreach my $index (@annot) {
         map {
            print OUT $_->{'variation'}->display_id," ", $_->{'variation'}->seq_region_name," ",$_->{'variation'}->var_class()," ",$_->{'variation'}->seq_region_start," ";
            print OUT $_->{'distance'}," ",$_->{'gene'}->stable_id,":",$_->{'gene'}->external_name(),":",$_->{'gene'}->biotype,"\n";
         } @$index;
      }
      @annot = ();
   }
}


exit;


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