#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Iterator::FastaDb;
use POSIX qw(mkfifo);
use Env qw(TMPDIR);
use File::Temp;

my $usage=<<'ENDHERE';
NAME:
extractMiRNATargetsGeneSequence.pl

PURPOSE:

INPUT:
--infile <string>    : tsv file with accession number in 2nd column.
--link <string>      : tsv with 2 columns. 1st col=genbank (GCA_...).
                       2nd col=fasta header.
--indir_fna <string> : Directory where are located fasta files.       

OUTPUT:
STDOUT <string>      : Same as infile. But with function added in last column.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
National Research Council Canada - Genomics and Microbiomes
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile, $link, $indir_fna);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'link=s'      => \$link,
   'indir_fna=s' => \$indir_fna,
   'verbose'     => \$verbose,
   'help'        => \$help
);
if ($help) { print $usage; exit; }

my $tmpdir = File::Temp->newdir(
    "tmpDir-SplitPairsGz-XXXXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 1
);

## MAIN
my %hash;
open(IN, "<".$link) or die "Can't open $link\n";
while(<IN>){
   chomp;
   my @row = split(/\t/, $_);
   my $header = $row[1];
   my $genbank_id = $row[0];
   
   $hash{$header}{name} = $genbank_id;
}
close(IN);

print STDERR Dumper(\%hash);

open(IN, "<".$infile) or die "Can't open $infile\n";
while(<IN>){
    chomp;
    my @row = split(/\t/, $_);
    my $header = $row[1];
    my $start = $row[-2];
    my $end = $row[-1];
    print STDERR $start." ".$end."\n";
  
    if(exists $hash{$header}){
        #print STDERR $header." exists\n";
        
        # Then do fna to get coordinates. 
        my $fasta_file = $indir_fna."/".$hash{$header}{name}.".fna.gz";
        if(!-e $fasta_file){
            print STDERR $fasta_file." missing!\n";
        }else{
            print STDERR "Processing $fasta_file\n";
        }
        #open(IN_FNA, "<".$fasta_file) or die "Can't open $fasta_file\n";
        my $pipe = "$tmpdir/reads.pipe";
        system("mkfifo $pipe");
        system("gunzip -c ".$fasta_file." > $pipe &");
        my $ref_fasta_db = Iterator::FastaDb->new($pipe) or die("Unable to open Fasta file, $pipe\n");
        #my $found = 0;
        while( my $curr = $ref_fasta_db->next_seq() ) {
            my ($fasta_header) = $curr->header =~ m/^>(\S+)/;
            print STDERR "fasta_header: ".$fasta_header."\n";
            print STDERR "header: ".$header."\n";
            if($header eq $fasta_header){#&& defined $hash{$header}{start} && defined($hash{$header}{end})){
                my $gene_seq = substr $curr->seq, ($start - 1), (($end - $start) + 1);
                print STDOUT ">".$header."_".$start."_".$end."\n".$gene_seq."\n";
                #$found = 1;
            }
        }
        system("rm $pipe");
    }
}
close(IN);
print STDERR "Completed.\n";
