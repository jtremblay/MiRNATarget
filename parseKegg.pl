#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
parseKegg.pl

PURPOSE:

INPUT:
--infile <string>         : output of blastp or diamond-blastp against kegg.
--ko <string>             : ko main file
--genetoko <string>       : file to do linking between gened id and ko id
--genetofunction <string> : Represents the genes.pep (faa) file, but with link to gene function.
--genes_desc <string>     : Tab file in which we have gene_id\tgene_description
                            The objective of this option is to populate genes that do not match
                            any KO in the ko file.

OUTPUT:
STDOUT


NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
National Research Council - Genomics and Microbiomes
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile, $ko, $genetoko, $genes_desc);
my $verbose = 0;

GetOptions(
   'infile=s'     => \$infile,
   'ko=s'         => \$ko,
   'genetoko=s'   => \$genetoko,
   'genes_desc=s' => \$genes_desc,
   'verbose'      => \$verbose,
   'help'         => \$help
);
if ($help) { print $usage; exit; }

## MAIN

my %hash;

# Parse blastp table
open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $geneId = $row[0];
   $hash{$geneId} = $row[1];
}
close(IN);

# Parse genetoko link table/file
my %hash_link;
my %hash_link_noprefix;
open(IN, "<".$genetoko) or die "Can't open $genetoko\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $geneId = $row[0];
   my $ko = $row[1];
   $ko =~ s/KO://;
   $hash_link{$geneId} = $ko;

   $geneId =~ s/\S\S\S\://;
   $hash_link_noprefix{$geneId} = $ko;

}
close(IN);

# Go through the genes_desc table
my %hash_genes_desc;
open(IN, "<".$genes_desc) or die "Can't open $genes_desc\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $geneId = $row[0];
   $hash_genes_desc{$geneId} = $row[1];
}
close(IN);

my %hash_ko;
my $counter = 0;
open(IN, "<".$ko) or die "Can't open $ko\n";
my $currKO;
my $currCat;
while(<IN>){
   chomp;

   if($_ =~ m/^ENTRY/){
      $currCat = "ENTRY";
   }elsif($_ =~ m/^NAME/){
      $currCat = "NAME";
   }elsif($_ =~ m/^DEFINITION/){
      $currCat = "DEFINITION";    
   }elsif($_ =~ m/^PATHWAY/){
      $currCat = "PATHWAY";
   }elsif($_ =~ m/^MODULE/){
      $currCat = "MODULE";
   }elsif($_ =~ m/DBLINKS/){
      $currCat = "DBLINKS";
   }elsif($_ =~ m/GENES/){
      $currCat = "GENES";
   }elsif($_ =~ m/BRITE/){
      $currCat = "BRITE";
   }
   #print STDERR $currCat."\t".$_."\n";

   if($currCat eq "ENTRY"){
      if($_ =~ m/ENTRY\s+(K\d{5})/){
         $currKO = $1;
         $hash_ko{$currKO}{ENTRY} = $1;
      }
   }
   
   if($currCat eq "NAME"){
      if($_ =~ m/NAME\s+(.*)$/){
         $hash_ko{$currKO}{NAME} = $1;
      }
   }
   
   if($currCat eq "DEFINITION"){
      if($_ =~ m/DEFINITION\s+(.*)\s\[.*\]$/){
          $hash_ko{$currKO}{DEFINITION} = $1;
      }elsif($_ =~ m/DEFINITION\s+(.*)/){
          $hash_ko{$currKO}{DEFINITION} = $1;
      }
   }

   if($currCat eq "PATHWAY"){
      if($_ =~ m/^.*(ko\d{5})\s+(.*)$/){
          if(exists $hash_ko{$currKO}{PATHWAY_ID}){
            $hash_ko{$currKO}{PATHWAY_ID} .= "==".$1;
          }else{
            $hash_ko{$currKO}{PATHWAY_ID} = $1;
          }
          if(exists $hash_ko{$currKO}{PATHWAY_DESC}){
            $hash_ko{$currKO}{PATHWAY_DESC} .= "==".$2;
          }else{
            $hash_ko{$currKO}{PATHWAY_DESC} = $2;
          }
      }
   }

   if($currCat eq "MODULE"){
      if($_ =~ m/^.*(M\d{5})\s+(.*)$/){
          if(exists $hash_ko{$currKO}{MODULE_ID}){
            $hash_ko{$currKO}{MODULE_ID} .= "==".$1;
          }else{
            $hash_ko{$currKO}{MODULE_ID} = $1;
          }
          if(exists $hash_ko{$currKO}{MODULE_DESC}){
            $hash_ko{$currKO}{MODULE_DESC} .= "==".$2;
          }else{
            $hash_ko{$currKO}{MODULE_DESC} = $2;
          }
      }
   }

   if($_ =~ m/^\\\\\\$/){
      $currCat = "";
      $counter++;
   }
   
   #last if($counter == 20);

}
close(IN);

#print STDERR Dumper(\%hash_ko);

# Print final table
my %hash_second_pass;

print STDOUT "#query\tkegg_gene_id\tKO_id\tNAME\tENTRY\tDEFINITION\tPATHWAY_ID\tPATHWAY_DESC\tMODULE_ID\tMODULE_DESC\n";
for my $geneId (keys %hash){
   my $kegg_gene_id = $hash{$geneId};

   if(exists $hash_link{$kegg_gene_id}){
       #print STDOUT $hash{$geneId}."\t".$geneId;
      print STDOUT $geneId."\t".$kegg_gene_id;
      my $ko = $hash_link{$kegg_gene_id};
      print STDOUT "\t".$ko;
      if(exists $hash_ko{$ko}{NAME}){         print STDOUT "\t".$hash_ko{$ko}{NAME};             }else{ print STDOUT "\tNULL";   }
      if(exists $hash_ko{$ko}{ENTRY}){        print STDOUT "\t".$hash_ko{$ko}{ENTRY};            }else{ print STDOUT "\tNULL";   }
      if(exists $hash_ko{$ko}{DEFINITION}){   print STDOUT "\t".$hash_ko{$ko}{DEFINITION};       }else{ print STDOUT "\tNULL";   }
      if(exists $hash_ko{$ko}{PATHWAY_ID}){   print STDOUT "\t".$hash_ko{$ko}{PATHWAY_ID};       }else{ print STDOUT "\tNULL";   }
      if(exists $hash_ko{$ko}{PATHWAY_DESC}){ print STDOUT "\t".$hash_ko{$ko}{PATHWAY_DESC};     }else{ print STDOUT "\tNULL";   }
      if(exists $hash_ko{$ko}{MODULE_ID}){    print STDOUT "\t".$hash_ko{$ko}{MODULE_ID};        }else{ print STDOUT "\tNULL";   }
      if(exists $hash_ko{$ko}{MODULE_DESC}){  print STDOUT "\t".$hash_ko{$ko}{MODULE_DESC}."\n"; }else{ print STDOUT "\tNULL\n"; }
   }else{
       print STDERR "ko with no id: ".$geneId." => ".$kegg_gene_id."\n";
       $hash_second_pass{$geneId} = $kegg_gene_id;
   }
}
#print STDERR Dumper(\%hash_link_noprefix);
#exit;

# Then get kegg gene ids that somehow does not contain abc: prefix...

print STDERR "Then get kegg gene ids that somehow does not contain abc: prefix...\n";

for my $geneId (keys %hash_second_pass){
    my $kegg_gene_id = $hash_second_pass{$geneId};

    if(exists $hash_link_noprefix{$kegg_gene_id}){
      print STDOUT $geneId."\t".$kegg_gene_id;
      my $ko = $hash_link_noprefix{$kegg_gene_id};
      print STDOUT "\t".$ko;
      if(exists $hash_ko{$ko}{NAME}){         print STDOUT "\t".$hash_ko{$ko}{NAME};             }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{ENTRY}){        print STDOUT "\t".$hash_ko{$ko}{ENTRY};            }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{DEFINITION}){   print STDOUT "\t".$hash_ko{$ko}{DEFINITION};       }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{PATHWAY_ID}){   print STDOUT "\t".$hash_ko{$ko}{PATHWAY_ID};       }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{PATHWAY_DESC}){ print STDOUT "\t".$hash_ko{$ko}{PATHWAY_DESC};     }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{MODULE_ID}){    print STDOUT "\t".$hash_ko{$ko}{MODULE_ID};        }else{ print STDOUT "\tNULL"; }
      if(exists $hash_ko{$ko}{MODULE_DESC}){  print STDOUT "\t".$hash_ko{$ko}{MODULE_DESC}."\n"; }else{ print STDOUT "\tNULL\n"; }
      delete($hash_second_pass{$geneId});
   }else{
       print STDERR "ko with still no id: ".$geneId." => ".$kegg_gene_id."\n";
   }
}        

# Then get all the kegg genes that do not have a KO in the ko file.
print STDERR "Then get all the kegg genes that do not have a KO in the ko file...\n";
for my $geneId (keys %hash_second_pass){
    my $kegg_gene_id = $hash_second_pass{$geneId};

    if(exists $hash_genes_desc{$kegg_gene_id}){
      print STDOUT $geneId."\t".$kegg_gene_id;
      #my $ko = $hash_link_noprefix{$kegg_gene_id};
      my $ko = "NULL";
      print STDOUT "\t".$ko;
      print STDOUT "\tNULL"; # name
      print STDOUT "\tNULL"; # entry
      print STDOUT "\t".$hash_genes_desc{$kegg_gene_id}; # definition
      print STDOUT "\tNULL";
      print STDOUT "\tNULL";
      print STDOUT "\tNULL";
      print STDOUT "\tNULL\n";
      delete($hash_second_pass{$geneId});
   }else{
       print STDERR "ko with still no id: ".$geneId." => ".$kegg_gene_id."\n";
   }
}        

exit;

