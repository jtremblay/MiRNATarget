#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
parseMiRNATargets.pl

PURPOSE:
From the output of companion script parseSsearch.pl, this script assigns a penalty score (Expect value) to each miRNA alignment according to the 
scoring scheme found in PMC3499661. Results of this script matches exactly the ones given 
by psRNATarget.

INPUT:
--infile <string>                : Alignments .tsv file format (output of companion script parseSsearch.pl). If no --infile <string> arg is 
                                   specified, will read from standard input.
--E_cutoff <float>               : default = 5.0  - Expectation value is the penalty for the mismatches 
                                   between mature small RNA and the target sequence. A higher value indicates 
                                   less similarity (and possibility) between small RNA and the target candidate. 
                                   The default penalty rule is set up by the scoring schema. Maximum expectation is the cutoff; 
                                   any small RNA-target pair with an expectation less than the cutoff will be discarded in the final result. 
                                   The recommended values are 3.0-5.0 depending on the scoring schema. 
--penalty_multiplier <string>    : default = 1.5. In the seed region (by default, 2-13 nt from the 5' of the miRNA strand), multiply mismatches by 
                                   penalty_multiplier. Only mismatches are multiplied, not G:U pairs and/or the '.' alignment
                                   caracters given by SSEARCH.
--hsp_start <int>                : default = 2 - Beginning of HSP region. --hsp_start and --hsp_end positions affect where in the alignment 
                                   the --penalty_multiplier will be applied.
--hsp_end <int>                  : default = 13 - End of HSP region.
--num_mismatch_seed <int>        : default = 2 - Maximum of allowed mismatches in the seed region, excluding G:U pairs.
--hsp_cutoff <int>               : default = 14 - HSPs (i.e. the one computed by SSEARCH) shorter than this value will be discarded.
--gap_cutoff <int>               : default = 1 - alignments having more than <int> gaps willbe discarded.
--total_mismatches_cutoff <int>  : default = 8 - alignments showing more than <int> mismatches will be discarded.
--GUs_cutoff <int>               : default = 7 - alignments having more than <int> mismatches will be discarded.
--keep_target_suffix             : set flag if you wish to keep the temporary _<int> suffix appended at the end
                                   of subjects IDs.
--rev                            : Set flag if alignments of miRNAs were done on revcomp subject sequences.
--verbose                        : Set flag for debugging.
--maximum_alignment_length <int> : default = 22 - alignments longer than this value will be discarded.
--extra_penalty_query_gap <int>  : default = 1. If gap is located on query (miRNA) sequence, add an extra penalty of <int>.
--alignment_length <int>         : default = 19 - penalty score for each alignment will be computed from 1 to <alignment_length>.

OUTPUT:
STDOUT <string>                  : standard output. Alignments that passed filters.
--outfile_failed <string>        : stanbard error. Alignments that failed to pass filters.

NOTES:
It is not entirely clear which ssearch36 parameters are used by psRNATarget, but using the following SSEARCH followed by the execution of this script with default parameters gives identical results to psRNATarget.
ssearch36 -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 -i <input_mirna.fasta> <reference_genome.fasta> > <output_file>
If you provide an already rev-complemented fasta file, you can omit the -i argument.

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

ENDHERE

## OPTIONS
my ($help, $infile, $E_cutoff, $num_mismatch_seed, $hsp_cutoff, $gap_cutoff, $total_mismatches_cutoff, $GUs_cutoff, $keep_target_suffix, $penalty_multiplier, 
    $rev, $maximum_alignment_length, $extra_penalty_query_gap, $outfile_failed, $keep_tmp_file, $hsp_start, $hsp_end, $alignment_length);
my $verbose = 0;

GetOptions(
  'infile=s'                   => \$infile,
  'E_cutoff=f'                 => \$E_cutoff,
  'num_mismatch_seed=i'        => \$num_mismatch_seed,
  'hsp_cutoff=i'               => \$hsp_cutoff,
  'gap_cutoff=i'               => \$gap_cutoff,
  'GUs_cutoff=i'               => \$GUs_cutoff,
  'penalty_multiplier=f'       => \$penalty_multiplier,
  'rev'                        => \$rev,
  'hsp_start=i'                => \$hsp_start,
  'hsp_end=i'                  => \$hsp_end,
  'keep_tmp_file'              => \$keep_tmp_file,
  'maximum_alignment_length=s' => \$maximum_alignment_length,
  'total_mismatches_cutoff=i'  => \$total_mismatches_cutoff,
  'keep_target_suffix'         => \$keep_target_suffix,
  'extra_penalty_query_gap=i'  => \$extra_penalty_query_gap,
  'outfile_failed=s'           => \$outfile_failed,
  'alignment_length=i'         => \$alignment_length,
  'verbose'                    => \$verbose,
  'help'                       => \$help
);
if ($help) { print $usage; exit; }

print STDERR "##########################################\n";
print STDERR "# Running parseMiRNATargets.pl ...       #\n";
print STDERR "##########################################\n";

$E_cutoff = 5.0 unless($E_cutoff);
$num_mismatch_seed = 2 unless($num_mismatch_seed);
$hsp_cutoff = 14 unless($hsp_cutoff);
$gap_cutoff = 1 unless($gap_cutoff);
$total_mismatches_cutoff = 8 unless($total_mismatches_cutoff);
$GUs_cutoff = 7 unless($GUs_cutoff);
$maximum_alignment_length = 22 unless($maximum_alignment_length);
$extra_penalty_query_gap = 1 unless($extra_penalty_query_gap);
$penalty_multiplier = 1.5 unless($penalty_multiplier);
$hsp_start = 2 unless($hsp_start);
$hsp_end = 13 unless($hsp_end);
if($hsp_start >= $hsp_end){
    die("--hsp_start <int> has to be greater than --hsp_end <int>\n");
}
if($hsp_start < 1){
    die("--hsp_start has to be 1 or greater, but smaller than --hsp_end <int>\n");
}
$alignment_length = 19 unless($alignment_length);
if($alignment_length < $hsp_end){
    die("--alignment_length <int> has to be >= than --hsp_end <int>\n");
}

my $OUT_FAILED;
if($outfile_failed){
    open($OUT_FAILED, ">".$outfile_failed) or die "Can't open $outfile_failed\n";
}

## MAIN
my %seen;
my %stats;

# Print header.
print STDOUT "#query_id\tsubject_id\tmatch_aln\tquery_aln\tsubject_aln\tq_start\tq_end\ts_start\ts_end\texpect_value\tstrand\n";

my $IN = new IO::File;
if($infile){
    open($IN, "<", $infile) or die "Can't open $infile\n";
}else{
    $IN = *STDIN;
}    
while(<$IN>){
    chomp;
    if ($_ =~ /^\s*$/) {
       next;#blank line
    }
    my @row = split(/\t/, $_);
    my $query       = $row[0];
    my $subject     = $row[1];
    my $aln_str     = $row[2];
    my $query_str   = $row[3];
    my $subject_str = $row[4];
    my $q_start     = $row[5];
    my $q_end       = $row[6];
    my $start       = $row[7];
    my $end         = $row[8];
    my $strand      = $row[9];
    my $hsp         = $row[10];
            
    #foreach my $subject (keys %{ $hash{$query} }) {
    next if($subject eq "length"); 
    my $subject_str_test = $row[4];
    $subject_str_test =~ s/\s+//g;
    if(length($subject_str_test) < 20 || length($subject_str_test) > 23){
        next;
    }

    # First avoid duplicates
    my $contig_id = ""; 
    if($subject =~ m/^(.*)_\d+$/){
       $contig_id = $1;
    }else{
        print STDERR $subject."\n";
        die("Problem parsing subject field...\n");
    }
    if(!exists $seen{ $strand."_".$query."_".$contig_id."_".$q_end."-".$q_start.":".$start."-".$end }){
        $seen{ $strand."_".$query."_".$contig_id."_".$q_end."-".$q_start.":".$start."-".$end }++;
    }else{
        $seen{ $strand."_".$query."_".$contig_id."_".$q_end."-".$q_start.":".$start."-".$end }++;
        next;
    }

    my $seed_region_start; my $seed_region_end;
    my $aln_str2; my $query_str2; my $subject_str2;
    
    print STDERR "Processing: ".$query."\t".$subject."\n" if($verbose);

    # For debug only. Should always use $hsp_config = 1
    my $hsp_config = 1; my $offset; my $end_substr;
    if($hsp_config == 1){
        $offset = 0;
        $end_substr = $alignment_length;
    }elsif($hsp_config == 2){
        $offset = 1;
        $end_substr = 20;
    }elsif($hsp_config == 3){
        $offset = 2;
        $end_substr = 20;
    }elsif($hsp_config == 4){
        $offset = length($query_str) - $q_start;
        $end_substr = 19;
    }

    # If alignments (5'-miRNA-3') vs 5'-DNA-3' are being processed, no need to reverse strings as they already are in the correct orientation.
    if($rev){
        $aln_str2       = substr($aln_str,     $offset, $end_substr);
        $query_str2     = substr($query_str,   $offset, $end_substr);
        $subject_str2   = substr($subject_str, $offset, $end_substr);
    
    # If alignments (3'-miRNA-5') vs 5'-DNA-3' are being processed, we have to reverse strings to make sure that we are starting from the 5' end.
    }else{
        print STDERR "     orig aln string:    ".$aln_str."\n" if($verbose);
        print STDERR "     query aln string:   ".$query_str."\n" if($verbose);
        print STDERR "     subject aln string: ".$subject_str."\n" if($verbose);
      
        $aln_str2      = reverse($aln_str);
        $query_str2    = reverse($query_str);
        $subject_str2  = reverse($subject_str);
        $aln_str2      = substr($aln_str2,     $offset, $end_substr);
        $query_str2    = substr($query_str2,   $offset, $end_substr);
        $subject_str2  = substr($subject_str2, $offset, $end_substr);
    }
    print STDERR "    length aln_str:      ".length($aln_str)."\n" if($verbose);
    print STDERR "    length query_str:    ".length($query_str)."\n" if($verbose);
    print STDERR "    length subject_str:  ".length($subject_str)."\n" if($verbose);
    print STDERR "    aln_str   :          ".$aln_str."\n" if($verbose);
    print STDERR "    query_str :          ".$query_str."\n" if($verbose);
    print STDERR "    aln_str   :          ".$subject_str."\n" if($verbose);
    print STDERR "    length aln_str2:     ".length($aln_str2)."\n" if($verbose);
    print STDERR "    length query_str2:   ".length($query_str2)."\n" if($verbose);
    print STDERR "    length subject_str2: ".length($subject_str2)."\n" if($verbose);
    print STDERR "    aln_str2   :         ".$aln_str2."\n" if($verbose);
    print STDERR "    query_str2 :         ".$query_str2."\n" if($verbose);
    print STDERR "    aln_str2   :         ".$subject_str2."\n" if($verbose);
    print STDERR "------------------------------------------\n" if($verbose); 

    # When looping through hsp, we have to keep the aln string part that matches to the actual hsp
    my @query_char = split(//, $query_str2);
    my @subject_char = split(//, $subject_str2);
    my @aln_char = split(//, $aln_str2);

    my $length = scalar(@aln_char);

    # Compute score according (more or less) fo Fahlgren 2009 (and psRNATarget).
    # Careful here because in psRNATarget, what they refer to as G:U pairs actually corresponds
    # to . characters in the alignment strings given by SSEARCH36. Furthermore, G:U pairs (or . pairs) 
    # are NOT multiplied by the penalty error multiplier (default 1.5) in the seed region (2-13 nt), only mismatches are multipled.
    # The following routine gives the exact same penalty scoring results
    # to the ones given by psRNATarget (which they call E or Expect value).
    my $z = 1;
    my $score = 0;
    my $gap_status_q = "closed";
    my $gap_status_s = "closed";
    my $gap_status = "closed";
    my $seed_mismatches = 0;
    my $total_mismatches = 0;
    my $gaps = 0;
    my $GUs = 0;
    my $gap_penalty = 0;
    my $total_extra_penalty_gap_query = 0;
    print STDERR "    ".$q_start."\n" if($verbose);
    for my $el (@aln_char){
        my $el2 = shift(@query_char);
        my $el3 = shift(@subject_char);

        my $is_GU = "no";
        my $curr_score = 0;
        my $mismatch_penalty = 0;
        my $is_mismatch = "no";
        
        if($el eq " " && $el2 ne "-" && $el3 ne "-"){
            $curr_score = 1;
            $is_mismatch = "yes";
            $total_mismatches++;
        }elsif($el eq " "  && $el2 eq "-" && $gap_status_q eq "closed"){ #gap opening
            $curr_score = 2;
            $gap_status_q = "open";
            $is_mismatch = "yes";
            $gaps++;
            $total_mismatches++;
            $total_extra_penalty_gap_query = $total_extra_penalty_gap_query + $extra_penalty_query_gap;
        }elsif($el eq " "  && $el3 eq "-" && $gap_status_s eq "closed"){ #gap opening
            $curr_score = 2;
            $gap_status_s = "open";
            $is_mismatch = "yes";
            $gaps++;
            $total_mismatches++;
        }elsif($el eq " "  && $el2 eq "-" && $gap_status_q eq "open"){ #gap extension
            $curr_score = 0.5;
            $is_mismatch = "yes";
            $gaps++;
            $total_mismatches++;
            $total_extra_penalty_gap_query = $total_extra_penalty_gap_query + $extra_penalty_query_gap;
        }elsif($el eq " "  && $el3 eq "-" && $gap_status_s eq "open"){ #gap extension
            $curr_score = 0.5;
            $is_mismatch = "yes";
            $gaps++;
            $total_mismatches++;
        }elsif($el eq "."){
            $curr_score = 0.5;
            $GUs++;
            $is_GU = "yes";
            $is_mismatch = "no";
            $gap_status_q = "closed";
            $gap_status_s = "closed";
        }elsif($el eq ":"){
            $gap_status_q = "closed";
            $gap_status_s = "closed";
        }else{
            print STDERR "    Did not match anything....: ".$el." ".$el2." ".$el3."\n" if($verbose);
            die;
        }
        
        if($z >= $hsp_start && $z <= $hsp_end){
            if($is_mismatch eq "yes"){
                print STDERR ("    is mismatch in the seed region - will apply penalty... curr_score: ".$curr_score." * 1.5 => ".($curr_score*$penalty_multiplier)."\n") if($verbose);
                $curr_score = $curr_score * $penalty_multiplier;
                $seed_mismatches++;
            }else{
                $curr_score = $curr_score * 1.0;
            }
        }
        $score = $score + $curr_score ;
        print STDERR "    $z:$curr_score    $score    $el2 $el $el3\n" if ($verbose);
        $z++;
    }
    $score = $score + $total_extra_penalty_gap_query;
    my $penalty_score = $score;

    if( $penalty_score <= $E_cutoff && 
        $seed_mismatches <= $num_mismatch_seed && 
        $hsp >= $hsp_cutoff &&
        $gaps <= $gap_cutoff &&
        $total_mismatches <= $total_mismatches_cutoff &&
        length($subject_str) <= $maximum_alignment_length){
        
        unless($keep_target_suffix){
            $subject =~ s/_\d+$//;
        }
       
        if($GUs <= $GUs_cutoff){
            print STDOUT $query."\t".$subject."\t".$aln_str."\t".$query_str."\t".$subject_str."\t".$q_start."\t".$q_end."\t".$start."\t".$end."\t".$penalty_score."\t".$strand."\n";
        }else{
            $stats{GU_reject}++;
        }
    }else{
        $stats{other_reject}++;
        if($outfile_failed){
            print $OUT_FAILED $query."\t".$subject."\t".$aln_str."\t".$query_str."\t".$subject_str."\t".$q_start."\t".$q_end."\t".$start."\t".$end."\t".$penalty_score."\t".$strand."\n";
        }
    }
}

close($OUT_FAILED) if($outfile_failed);
close($IN);
print STDERR Dumper(\%stats) if($verbose);
print STDERR  Dumper(\%seen) if($verbose);


sub complement_IUPAC {
   my $dna = shift;

   # complement the reversed DNA sequence
   my $comp = $dna;
   $comp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy-/TVGHCDKNYSAABWXRtvghcdknysaabwxr-/;
   return $comp;
}

exit;


