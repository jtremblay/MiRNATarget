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
addFunctionToMiRNATargets.pl

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
my ($help, $infile, $link, $indir_gff, $indir_fna);
my $verbose = 0;

GetOptions(
   'infile=s'    => \$infile,
   'link=s'      => \$link,
   'indir_gff=s' => \$indir_gff,
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
    my $hit_start = $row[5];
    my $hit_end = $row[6];
  
    if(exists $hash{$header}){
        #print STDERR $header." exists\n";
        # First do gff to get coordinates. 
        my $genbank_file = $indir_gff."/".$hash{$header}{name}.".gff";
        if(!-e $genbank_file){
            print STDERR $genbank_file." missing!\n";
            next;
        }else{
            #next;
            print STDERR "Processing $genbank_file\n";
        }
        open(IN_GFF, "<".$genbank_file) or die "Can't open $genbank_file\n";
        while(my $line = <IN_GFF>){
            next if($line =~ m/^#/);
            chomp;
            my @row = split(/\t/, $line);
            next if($row[2] ne "gene");
            if($row[0] eq $header){
                # find coords.
                my $start = $row[3];
                my $end = $row[4];
                if($hit_start >= $start && $hit_end <= $end && $row[2] eq "gene"){
                    print STDOUT $_."\t".$start."\t".$end."\n";
                    $hash{$header}{start} = $start;
                    $hash{$header}{end} = $end;
                }#elsif($hit_start <= $start && $hit_end <= $end && $row[2] eq "gene"){#else if hit outside boundaries...
                #    print STDOUT $_."\t".$start."\t".$end."\n";
                #    $hash{$header}{start} = $start;
                #    $hash{$header}{end} = $end;
                #}elsif($hit_start >= $start && $hit_end >= $end && $row[2] eq "gene"){#else if hit outside boundaries...
                #    print STDOUT $_."\t".$start."\t".$end."\n";
                #    $hash{$header}{start} = $start;
                #    $hash{$header}{end} = $end;
                #}else{
                #    print STDERR $_."\n".$start."\t".$end."\t".$row[2]."\n".$hit_start."\t".$hit_end."\n";
                #}
            }
        }
        close(IN_GFF);
        
#        # Then do fna to get coordinates. 
#        my $fasta_file = $indir_fna."/".$hash{$header}{name}.".fna.gz";
#        if(!-e $fasta_file){
#            print STDERR $fasta_file." missing!\n";
#        }else{
#            print STDERR "Processing $fasta_file\n";
#        }
#        #open(IN_FNA, "<".$fasta_file) or die "Can't open $fasta_file\n";
#        my $pipe = "$tmpdir/reads.pipe";
#        system("mkfifo $pipe");
#        system("gunzip -c ".$fasta_file." > $pipe &");
#        my $ref_fasta_db = Iterator::FastaDb->new($pipe) or die("Unable to open Fasta file, $pipe\n");
#        my $found = 0;
#        while( my $curr = $ref_fasta_db->next_seq() ) {
#            my ($fasta_header) = $curr->header =~ m/^>(\S+)/;
#            print STDERR "fasta_header: ".$fasta_header."\n";
#            print STDERR "header: ".$header."\n";
#            if($header eq $fasta_header && defined $hash{$header}{start} && defined($hash{$header}{end}) && $found == 0){
#                my $gene_seq = substr $curr->seq, ($hash{$header}{start}-1), (($hash{$header}{end} - $hash{$header}{start})+1);
#                print STDOUT $gene_seq."\n";
#                $found = 1;
#            }
#        }
#        system("rm $pipe");
    }
}
close(IN);

