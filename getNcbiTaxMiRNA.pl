#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
getNcbiTaxMiRNA.pl

PURPOSE:

INPUT:
--mirna_table <string>      : blast output table
--accession_to_tax <string> : Link between access id and tax id.
--names <string>            : names (names.dmp) file

OUTPUT:
STDOUT

NOTES:
#~/scripts/getNcbiTax.pl --mirna_table ./myblast.table --accession_to_tax ./gi_taxid_nucl.dmp --names ../ncbi_tax2/taxa.sqlite.dmp > ./stdout.txt 2> stderr.txt

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
National Research Council - Biomonitoring
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $mirna_table, $accession_to_tax, $names);
my $verbose = 0;

GetOptions(
   'mirna_table=s'         => \$mirna_table,
   'accession_to_tax=s'    => \$accession_to_tax,
   'names=s'               => \$names,
   'verbose'               => \$verbose,
   'help'                  => \$help
);
if ($help) { print $usage; exit; }

die "--mirna_table missing\n" unless($mirna_table);
die "--accession_to_tax missing\n" unless($accession_to_tax);
die "--names missing\n" unless($names);


## MAIN
my %hash_accession;
my %hash_accession_to_tax;

open(IN, "<".$mirna_table) or die "Can't open $mirna_table\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $subject = $row[1];
   $subject =~ s/ .*$//;
   $hash_accession{$subject} = "";
}
close(IN);

print STDERR "Done reading mirna table\n";
print STDERR Dumper(\%hash_accession);

my %hash_tax;

open(IN, "<".$accession_to_tax) or die "Can't open $accession_to_tax\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $accession = $row[1];
   my $tax_id = $row[2];
   if(exists $hash_accession{$accession}){
      $hash_tax{$tax_id}{ACCESSION} = $accession;
      $hash_accession_to_tax{$accession}{TAXID} = $tax_id;
   }
}
close(IN);
print STDERR "Done reading accession_to_tax table\n";
print STDERR Dumper(\%hash_tax);

## Example of file format
#INSERT INTO "species" VALUES(1580,1578,'Lactobacillus brevis','species','1580,1578,33958,186826,91061,1239,1783272,2,131567,1');
#
#1580,1578;Lactobacillus brevis;species;1580,1578,33958,186826,91061,1239,1783272,2,131567,1
#
#1580,1578;
#Lactobacillus brevis;
#species;
#1580,1578,33958,186826,91061,1239,1783272,2,131567,1


my %hash_names;
open(IN, "<".$names) or die "Can't open $names\n";
while(<IN>){
   chomp;
   next if $_ =~ m/synonym/;
   my $string = $_;
   
   my $alt_rank = "null";
   if($string =~ m/INSERT INTO "(\S+)"/){
     $alt_rank = $1;
   }
   #print STDERR "Alt rank:".$alt_rank."\n";
   
   if($string =~ m/VALUES\((.*)\);/){
      my $match = $1; 
      #print STDERR "match before: ".$match."\n";
      $match =~ s/',/';/g;
      $match =~ s/,'/;'/g;
      $match =~ s/'//g;
      #print STDERR $match."\n";
      my @row = split(/;/, $match);
   
      my $row0 = $row[0];
      my @row0 = split(/,/, $row0); #  1122198,48076
      my $tax_id = $row0[0];        #  OK
      my $name = $row[1];           #  Marinobacterium georgiense DSM 11526
      my $rank = $row[2];           #  no rank 
      my $lineage = $row[3];        #  1122198,48076,48075,72275,135622,1236,1224,2,131567,1
   
      $hash_names{$tax_id}{NAME} = $name;
      #if($rank eq "no rank" && $alt_rank eq "null"){
      #  $hash_names{$tax_id}{RANK} = $rank;
      #}elsif($rank eq "no rank" && $alt_rank ne "null"){ # alt_rank should be "species"
      #  $hash_names{$tax_id}{RANK} = $alt_rank;
      #}elsif($rank ne "no rank" && $alt_rank ne "null"){
      #  $hash_names{$tax_id}{RANK} = $rank;
      #} 
      $hash_names{$tax_id}{RANK} = $rank;
      $hash_names{$tax_id}{LINEAGE} = $lineage;
    
   }
   #exit if($. > 10);
}
close(IN);
print STDERR "Done reading names sql dmp file\n";
print STDERR Dumper(\%hash_names);

my %final_hash;

for my $tax_id (keys %hash_tax){
   
   if(exists $hash_names{$tax_id}){
  
      #print STDOUT $hash_names{$tax_id}{LINEAGE}."\n"; 
      my @lineage = split(/,/, $hash_names{$tax_id}{LINEAGE});
      foreach (@lineage){
         next if($_ eq "131567"); # this is the tax id for cellular organisms and should be excluded, because all species will be assigned 'cellular organisms'
         if(exists $hash_names{$_}){
             #print STDOUT $hash_names{$_}{RANK}."=>".$hash_names{$_}{NAME}.";";

            if($hash_names{$_}{RANK} eq "superkingdom"){
               $final_hash{$tax_id}{kingdom} = $hash_names{$_}{NAME};

            }elsif($hash_names{$_}{RANK} eq "phylum"){
               $final_hash{$tax_id}{phylum} = $hash_names{$_}{NAME};

            }elsif($hash_names{$_}{RANK} eq "class"){
               $final_hash{$tax_id}{class} = $hash_names{$_}{NAME};

            }elsif($hash_names{$_}{RANK} eq "order"){
               $final_hash{$tax_id}{order} = $hash_names{$_}{NAME};

            }elsif($hash_names{$_}{RANK} eq "family"){
               $final_hash{$tax_id}{family} = $hash_names{$_}{NAME};

            }elsif($hash_names{$_}{RANK} eq "genus"){
               $final_hash{$tax_id}{genus} = $hash_names{$_}{NAME};

            }elsif($hash_names{$_}{RANK} eq "species"){
               $final_hash{$tax_id}{species} = $hash_names{$_}{NAME};
            }
         }
      }
   }
}
print STDERR "Done generating final hash\n";
#print STDERR Dumper(\%final_hash);

#print STDOUT "gene_id\tAccession\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";

open(IN, "<".$mirna_table) or die "Can't open $mirna_table\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $subject = $row[1];
   $subject =~ s/ .*$//;
   #$subject =~ m/gi\|(\d+)\|/;
   my $accession = $subject;
   my $geneid = $row[0];
     
   next if(!defined $hash_accession_to_tax{$accession}{TAXID}); 
   my $tax_id = $hash_accession_to_tax{$accession}{TAXID};
   next if(!defined $tax_id);
   next if($tax_id eq "");

   print STDOUT $_."\t";

   #for my $tax_id (keys %final_hash){
   print STDOUT $hash_tax{$tax_id}{ACCESSION}."\t";
   
   if(defined $final_hash{$tax_id}{kingdom}) { print STDOUT $final_hash{$tax_id}{kingdom}."\t";}  else { print STDOUT "NULL\t"; }
   if(defined $final_hash{$tax_id}{phylum})  { print STDOUT $final_hash{$tax_id}{phylum}."\t";}   else { print STDOUT "NULL\t"; }
   if(defined $final_hash{$tax_id}{class})   { print STDOUT $final_hash{$tax_id}{class}."\t";}    else { print STDOUT "NULL\t"; }
   if(defined $final_hash{$tax_id}{order})   { print STDOUT $final_hash{$tax_id}{order}."\t";}    else { print STDOUT "NULL\t"; }
   if(defined $final_hash{$tax_id}{family})  { print STDOUT $final_hash{$tax_id}{family}."\t";}   else { print STDOUT "NULL\t"; }
   if(defined $final_hash{$tax_id}{genus})   { print STDOUT $final_hash{$tax_id}{genus}."\t";}    else { print STDOUT "NULL\t"; }
   if(defined $final_hash{$tax_id}{species}) { print STDOUT $final_hash{$tax_id}{species};}       else { print STDOUT "NULL"; }
   print STDOUT "\n";
}
close(IN);
