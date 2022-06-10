package ZZ::EnsemblGene;
#
# generate ZZ::EnsemblGene  - generate a gene object from an entry of Homo_sapiens.GRCh37.87.chr.gtf
#
#  Zhanyang Zhu
#   
use strict;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use Carp;
use File::Basename;
# 
# The whole gene annotation data file can be downloaded from Ensembl: http://ftp.ensembl.org/pub/grch37/release-106/gtf/homo_sapiens/
# select file - Homo_sapiens.GRCh37.87.chr.gtf.gz
#
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
#An example of gene entry:
#11      ensembl_havana  gene    2465914 2870339 .       +       .       gene_id "ENSG00000053918"; gene_version "11"; gene_name "KCNQ1"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
#m
#21      havana  gene    48073470        48073934        .       -       .       gene_id "ENSG00000230982"; gene_version "1"; gene_name "DSTNP1"; gene_source "havana"; gene_biotyp                     e "pseudogene";
#21      ensembl_havana  gene    48110676        48111157        .       +       .       gene_id "ENSG00000212932"; gene_version "2"; gene_name "RPL23AP4"; gene_source "ensembl_ha                     vana"; gene_biotype "pseudogene";
#MT      insdc   gene    577     647     .       +       .       gene_id "ENSG00000210049"; gene_version "1"; gene_name "MT-TF"; gene_source "insdc"; gene_biotype "Mt_tRNA";
#MT      insdc   gene    648     1601    .       +       .       gene_id "ENSG00000211459"; gene_version "2"; gene_name "MT-RNR1"; gene_source "insdc"; gene_biotype "Mt_rRNA";
#
#
my @GeneralInfoKeys = qw{seqname source feature start end score strand frame attribute};
my %AllowedGeneralInfoKeys = map { $_ => 1 } @GeneralInfoKeys;
#
my @GeneInfoKeys = qw{chrom pos_start pos_end strand gene_id gene_version gene_name gene_source gene_biotype};
my %AllowedGeneInfoKeys = map { $_ => 1 } @GeneInfoKeys;
#
# Gene object constructor:
# The gene object is invalid if any paring error occurs
#
sub new {
  my $class = shift;
  my $self = {};
  bless($self, $class);
  my $line = shift;
  $self = $self->parseGeneEntry($line);
  $self->{valid} = 1;
  return $self;
}

#
# Create the gene object by parsing the gene line entry:
#
sub parseGeneEntry { 
  my $self = shift;
  my $line = shift;
  chomp($line);
  my (@data) = split('\t', $line);
  my $nd = @data;
  for (my $i=0; $i<$nd; $i++){
    $data[$i] = trim($data[$i]);
  }
  if ($data[2] ne 'gene') { 
     $self->{valid} = 0;
     print "Error - not gene entry (it was <$data[2]>\n";
     print "Erro  - line=$line\n";
     return $self;
  }
  $self->{data}->{$GeneInfoKeys[0]} = $data[0];   # chrom
  # skip - source feature
  $self->{data}->{$GeneInfoKeys[1]} = $data[3];   # pos_start
  $self->{data}->{$GeneInfoKeys[2]} = $data[4];   # pos_end
  # skip - score
  $self->{data}->{$GeneInfoKeys[3]} = $data[6];   # strand
  # skip - frame
  # attribute:
  my (@att) = split(/;/, $data[8]);  # gene_id "ENSG00000053918"; gene_version "11"; ...
  my $natt = @att; 
  for (my $i=0; $i<$natt; $i++) { 
    my $tmpS = trim($att[$i]);
    my ($k, $v) = split(/\s/, $tmpS);
    $k = trim($k);
    ($v, undef) = extractQ($v); # remove quotation marks
    $v = trim($v);
    if (exists($AllowedGeneInfoKeys{$k})) {
      $self->{data}->{$k} = $v;
    } else { 
      print "Error - unrecognized gene attribute <$k> with value <$v>\n";
      print "Error - line=$line\n";
      $self->{valid} = 0;
    }
  }   
  # make sure all require gene info is here:
  my $nk = @GeneInfoKeys;
  for (my $i=0; $i<$nk; $i++) { 
    if (!(exists($self->{data}->{$GeneInfoKeys[$i]}))) { 
       print "Error - missing gene info $GeneInfoKeys[$i]\n";
       print "Error - line=$line\n"; 
       $self->{valid} = 0;
    }
  }
  return $self;
}

#
# Given a gene info key, return its value:
# Make sure the key is allowed
#
sub getItem { 
  my $self = shift;
  my $item = shift;
  if (exists($AllowedGeneInfoKeys{$item})) { 
     my $v = $self->{data}->{$item};
     return $v;
  } else { 
     return undef;
  }
} 

# Assign a value for gene info key:
# only certain keys are allowed
sub setItem { 
  my $self = shift;
  my $item = shift;
  my $value = shift;
  if (exists($AllowedGeneInfoKeys{$item})) { 
     $self->{data}->{$item} = $value;
  } else {
    print "Error: not allowed key: <$item>\n";
    exit 1;
  }
  return $self;
}

#
#print all the gene information with key/value pair:
#
sub printGeneInfo { 
  my $self = shift;
  my $nk = @GeneInfoKeys;
  for (my $i=0; $i<$nk; $i++) {
    my $k = $GeneInfoKeys[$i];
    my $v = $self->{data}->{$GeneInfoKeys[$i]};
    print "$k=$v; "
  }
  print "\n";
} 

#
# Parse and set an item as key and value pair:
# If key is provided, check if the string contains the key
# 
sub parseSetItem{
  my $self = shift;
  my $line = shift;
  my $selecteditemKey = shift;
  my ($itemKey, $v) = $self->parseItem($line);
  #print "$itemKey, $v \n";
  #print "$line\n";
  $self->setItem($itemKey, $v);
  if (defined($selecteditemKey)) { 
    if ($selecteditemKey ne $itemKey) { 
       print "Error: line keys do not match (<$selecteditemKey> vs <$itemKey>)\n";
       print $line . "\n";
       exit 1;
    }
  }
  return $self;
}

#
# Parse out the fisrt key sperated by ":", return the rest as value:
# The value may contain ":"
#
sub parseItem { 
  my $self = shift;
  my $line = shift;;
  chomp($line);
  my ($itemKey, @tmp) = split(/:/, $line); # split only by first ":"
  my $v = join(":", @tmp);  # rejoin just in case name contains ":"
  $v=trim($v);
  return ($itemKey, $v);
}

#
# extract string between the first " and the last ":
#
sub extractQ {
  my $oneString = shift;
  my $i1 = index $oneString, qq(");
  my $i2 = rindex $oneString, qq(");
  if ($i1 != - 1 && $i2-$i1>1) {
     my $Qs = substr $oneString, $i1+1, $i2-$i1-1;
     my $RestS = substr $oneString, $i2+1;
     $Qs = trim($Qs);
     $RestS = trim($RestS);
     return ($Qs, $RestS);
  } else {
     return (undef, undef);
  }
}

# trimming functions:
sub ltrim { my $s = shift; chomp($s); if (!empty($s)) {$s =~ s/^\s+//;     };  return $s };
sub rtrim { my $s = shift; chomp($s); if (!empty($s)) {$s =~ s/\s+$//;     };  return $s };
sub  trim { my $s = shift; chomp($s); if (!empty($s)) {$s =~ s/^\s+|\s+$//g}; return $s };

sub isOnlyDot {
     my $s = shift;
     $s = trim($s);
     if ($s =~ /^\.$/) {
        return 1;
     } else {
        return 0;
     }
}

# check if string is empty string
sub empty {
     my $s = shift;
     if ($s =~ /^(\s)*$/) {
        return 1;
     } else {
        return 0;
     }
}

1;
