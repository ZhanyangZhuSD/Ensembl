#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use File::Basename;
use Carp;
use Getopt::Long;
my $help = 0;
my $EnsemblFile = "Homo_sapiens.GRCh37.87.chr.gtf"; #  GRCh37 gene info
# ftp site: http://ftp.ensembl.org/pub/grch37/release-106/gtf/homo_sapiens/
# select file - Homo_sapiens.GRCh37.87.chr.gtf.gz

#https://uswest.ensembl.org/info/website/upload/gff.html
#Fields
#Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
#
#seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#source - name of the program that generated this feature, or the data source (database or project name)
#feature - feature type name, e.g. Gene, Variation, Similarity
#start - Start position* of the feature, with sequence numbering starting at 1.
#end - End position* of the feature, with sequence numbering starting at 1.
#score - A floating point value.
#strand - defined as + (forward) or - (reverse).
#frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
#*- Both, the start and end position are included. For example, setting start-end to 1-2 describes two bases, the first and second base in the sequence.
#
#Note that where the attributes contain identifiers that link the features together into a larger structure, these will be used by Ensembl to display the features as joined blocks.
# An example of gene entry:
# 11      ensembl_havana  gene    2465914 2870339 .       +       .       gene_id "ENSG00000053918"; gene_version "11"; gene_name "KCNQ1"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
my $debug=0;
#
my @GeneralInfoKeys = qw{seqname source feature start end score strand frame attribute};
my %AllowedGeneralInfoKeys = map { $_ => 1 } @GeneralInfoKeys;
#
my @GeneInfoKeys = qw{pos_start pos_end gene_version chrom strand gene_id gene_name gene_source gene_biotype};
my %AllowedGeneInfoKeys = map { $_ => 1 } @GeneInfoKeys;
#
BEGIN { # add the pm path to the library search path: 
      #my $module_path = abs_path(dirname($0)) . "/pm";
      #my $module_path = abs_path(dirname($0));
      my $module_path = abs_path(dirname($0)) . "/../pm";
      unshift @INC, $module_path;
}

use ZZ::EnsemblGene;

MAIN:
{
   GetOptions ('f|EnsemblFile=s'=> \$EnsemblFile,
               'h|help'		=> \$help);
   if ($help) { 
     print "$0 -f EnsemblFile\n";
     exit;
   }
   if (! -e $EnsemblFile) { 
     print "Error - EnsemblFile file does not exist <$EnsemblFile>\n";
     print "$0 -f EnsemblFile\n";
     exit 1;
   }
   my $ESBFH;
   open ($ESBFH, "<$EnsemblFile");
  # read in all lines:
   my $ESBGenes = {}; # key is gene_id (ENSG00000223972)
   my $nline = 0;
   my $ngene = 0;
   while (my $line=<$ESBFH>) { 
     chomp($line);
     if ($line =~ m/^#!/) {
      next;    # skip comment line
     }
     my (@data) = split('\t', $line);
     my $nd = @data; 
     for (my $i=0; $i<$nd; $i++){ 
        $data[$i] = trim($data[$i]);
     }
#feature - feature type name, e.g. Gene, Variation, Similarity
# identify the gene entry from other entries in the same file:
     if ($data[2] eq 'gene') { 
       $nline++;
       my $gene = ZZ::EnsemblGene->new($line); 
       #print "$line\n";
       #$gene->printGeneInfo();
       if ($gene->{valid}) { # only keep valid entry:
         $ngene++;
         my $gene_id = $gene->getItem('gene_id');
         if (exists($ESBGenes->{$gene_id})) { 
            print "Error - duplicated gene_id <$gene_id>\n";
         } else {   
            $ESBGenes->{$gene_id} = $gene;
         }
       }
     }
   }   
   close($ESBFH);
   print "Number of lines = $nline; Number of Genes = $ngene\n";
   # find out all gene_source gene_biotype values:
   my $gene_source = {};
   my $gene_biotype = {};
   foreach my $gid (keys %$ESBGenes) {
     my $gs = $ESBGenes->{$gid}->getItem('gene_source');
     if ( !exists($gene_source->{$gs}) ) { 
       $gene_source->{$gs} = 1;
     }  
     my $gbt = $ESBGenes->{$gid}->getItem('gene_biotype');
     if ( !exists($gene_biotype->{$gbt}) ) {
       $gene_biotype->{$gbt} = 1;
     }
   }
   my @gss = keys (%$gene_source);
   my @gbts = keys (%$gene_biotype);
   
   print "All gene_source values: @gss\n";
   print "All gene_biotype values @gbts\n";

   # ensembl_havana ensembl insdc havana
   # IG_V_pseudogene TR_V_pseudogene lincRNA pseudogene snRNA IG_C_pseudogene rRNA TR_J_pseudogene TR_D_gene IG_D_gene TR_V_gene IG_J_pseudogene polymorphic_pseudogene TR_C_gene IG_J_gene sense_intronic IG_C_gene 3prime_overlapping_ncrna Mt_rRNA sense_overlapping IG_V_gene processed_transcript snoRNA TR_J_gene antisense Mt_tRNA miRNA protein_coding misc_RNA
}

#
# remove \' for database loading:
#
sub removeSingleQuota { 
  my $text = shift;
  $text =~ s/\'//g;
  return $text;
}

#
#formt single quota into double single quota for database loading use
#
sub formatSingleQuota {
  my $text = shift;
  $text =~ s/\'/\'\'/g;
  return $text;
}

#
# remove all white space from a string:
#
sub removeAllSpaces { 
  my $text = shift;
  $text =~ s/\s//g;
  return $text;
}
#
# format date of #2015-08-07T03:39:28Z
#
sub formatDate {
  my $date = shift; #2015-08-07T03:39:28Z
  $date =~ s/\"//g;
  $date = trim($date);
  my ($date, $time) = split(/[TZ]/, $date);
  return ($date);
}

sub ltrim { my $s = shift; chomp($s); if (!empty($s)) {$s =~ s/^\s+//;     };  return $s };
sub rtrim { my $s = shift; chomp($s); if (!empty($s)) {$s =~ s/\s+$//;     };  return $s };
sub  trim { my $s = shift; chomp($s); if (!empty($s)) {$s =~ s/^\s+|\s+$//g}; return $s };
sub empty {
     my $s = shift;
     if ($s =~ /^(\s|,)*$/) {
        return 1;
     } else {
        return 0;
     }
}

