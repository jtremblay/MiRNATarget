#!/usr/bin/env python

"""Python script that parses miRNA alignment tsv file (generated with companion script parse_ssearch.py) and compute the psRNATarget penalty scoring scheme for each alignment.
Developed and tested with python 3.9.0

Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
"""

import argparse
import os
import sys
import re
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def parse_command_line_arguments():
    parser = argparse.ArgumentParser(description='Parse miRNA targets')
    parser.add_argument('-i', '--infile', required=False, help='Input file (i.e. output of parse_ssearch.py). This argument is optional as the output of parse_ssearch.py can be piped directly to this script as well.', type=argparse.FileType('r'))
    parser.add_argument('--E_cutoff', type=float, default=5.0, help='E-value cutoff')
    parser.add_argument('--num_mismatch_seed', type=int, default=2, help='Number of seed mismatches')
    parser.add_argument('--hsp_cutoff', type=int, default=14, help='HSP cutoff')
    parser.add_argument('--gap_cutoff', type=int, default=1, help='Gap cutoff')
    parser.add_argument('--GUs_cutoff', type=int, default=7, help='GUs cutoff')
    parser.add_argument('--penalty_multiplier', type=float, default=1.5, help='Penalty multiplier')
    parser.add_argument('--seed_start', type=int, default=2, help='Start of the seed region')
    parser.add_argument('--seed_end', type=int, default=13, help='End of the seed region')
    parser.add_argument('--alignment_length', type=int, default=19, help='Length of the alignment to be scores. In other words, how many bp from the 5p end of the miRNA sequence will be considered in the final penalty score')
    parser.add_argument('--rev', default=False, action=argparse.BooleanOptionalAction, help='Reverse')
    parser.add_argument('--maximum_alignment_length', type=int, default='22', help='Maximum alignment length')
    parser.add_argument('--total_mismatches_cutoff', type=int, default=8, help='Total mismatches cutoff')
    parser.add_argument('--keep_target_suffix', default=False, action=argparse.BooleanOptionalAction, help='Keep target suffix')
    parser.add_argument('--extra_penalty_query_gap', type=int, default=1, help='Extra penalty for query gap')
    parser.add_argument('-o', '--outfile_failed', type=argparse.FileType('w'), help='Output file for failed alignments')
    parser.add_argument('--verbose', default=False, action=argparse.BooleanOptionalAction, help='Verbose output')
    args = parser.parse_args()
    return args

def main(arguments):
    args = parse_command_line_arguments()
    seed_start = args.seed_start
    seed_end = args.seed_end
    alignment_length = args.alignment_length
    outfile_failed = args.outfile_failed
    verbose = args.verbose
    rev = args.rev
    penalty_multiplier = args.penalty_multiplier
    hsp_cutoff = args.hsp_cutoff
    gap_cutoff = args.gap_cutoff
    GUs_cutoff = args.GUs_cutoff
    E_cutoff = args.E_cutoff
    num_mismatch_seed = args.num_mismatch_seed
    extra_penalty_query_gap = args.extra_penalty_query_gap
    total_mismatches_cutoff = args.total_mismatches_cutoff
    keep_target_suffix = args.keep_target_suffix
    maximum_alignment_length = args.maximum_alignment_length


    if seed_start >= seed_end:
        raise Exception("--seed_start <int> has to be greater than --seed_end <int>")

    if seed_start < 1 :
        raise Exception("--seed_start has to be 1 or greater, but smaller than --seed_end <int>")
    
    if alignment_length < seed_end:
        raise Exception("--alignment_length <int> has to be >= than --seed_end <int>")

    fhand_out_failed = ""
    if outfile_failed:
        try:
            fhand_out_failed = open(outfile_failed, 'wb')
        except OSError:
            print("Could not open/write file:", outfile_failed)
            sys.exit()

    seen = {}
    stats = {}
    
    direction = "rev" if args.rev else "fwd"

    print("#query_id\tsubject_id\tmatch_aln\tquery_aln\tsubject_aln\tq_start\tq_end\ts_start\ts_end\texpect_value\tstrand", file=sys.stdout)
 
    # Decides if reading from file or standard input (piped)
    if args.infile:
        infile = os.path.abspath(args.infile.name)
    
    fhand = open(infile, 'r') if args.infile else sys.stdin
    
    for line in fhand:
        line = line.strip()
    
        row = line.split("\t")
        query       = row[0]
        subject     = row[1]
        aln_str     = row[2]
        query_str   = row[3]
        subject_str = row[4]
        q_start     = int(row[5])
        q_end       = int(row[6])
        start       = int(row[7])
        end         = int(row[8])
        strand      = row[9]
        hsp         = int(row[10])
        #print(line)
        if line.startswith("#"):
            continue
    
        subject_str_test = row[4]
        subject_str_test = subject_str_test.replace(" ", "")
        if(len(subject_str_test) < 20 or len(subject_str_test) > 23):
            continue

        # First avoid duplicates
        contig_id = ""
        match = re.match(r"^(.*)_\d+$", subject)
        if match:
            contig_id = match.group(1)
        else:
            print(subject, file=sys.stderr)
            raise Exception("Problem parsing subject field...")
       
        if strand+"_"+query+"_"+contig_id+"_"+str(q_end)+"-"+str(q_start)+":"+str(start)+"-"+str(end) not in seen:
            seen[ strand+"_"+query+"_"+contig_id+"_"+str(q_end)+"-"+str(q_start)+":"+str(start)+"-"+str(end) ] = 1
        else:
            seen[ strand+"_"+query+"_"+contig_id+"_"+str(q_end)+"-"+str(q_start)+":"+str(start)+"-"+str(end) ] += 1
            continue
    
        #seed_region_start = None
        #seed_region_end = None
        aln_str2 = str
        query_str2 = str
        subject_str2 = str
        
        if verbose: print("Processing: "+query+"\t"+subject+"\n", file=sys.stderr)

        # For debug only. Should always use hsp_config = 1
        hsp_config = 1
        offset = int
        end_substr = int
        if hsp_config == 1:
            offset = 0
            end_substr = alignment_length
        elif hsp_config == 2:
            offset = 1
            end_substr = 20
        elif hsp_config == 3:
            offset = 2
            end_substr = 20
        elif hsp_config == 4:
            offset = len(query_str) - q_start
            end_substr = 19

        # If alignments (5'-miRNA-3') vs 5'-DNA-3' are being processed, no need to reverse strings as they already are in the correct orientation.
        if rev:
            aln_str2       = aln_str[offset: end_substr]
            query_str2     = query_str[offset: end_substr]
            subject_str2   = subject_str[offset: end_substr]
        # If alignments (3'-miRNA-5') vs 5'-DNA-3' are being processed, we have to reverse strings to make sure that we are starting from the 5' end.
        else:
            if verbose: print("     orig aln string:    "+aln_str, file=sys.stderr)
            if verbose: print("     query aln string:   "+query_str, file=sys.stderr)
            if verbose: print("     subject aln string: "+subject_str, file=sys.stderr)
          
            aln_str2      = aln_str[::-1]
            query_str2    = query_str[::-1]
            subject_str2  = subject_str[::-1]
            aln_str2      = aln_str2[offset: end_substr]
            query_str2    = query_str2[offset: end_substr]
            subject_str2  = subject_str2[offset: end_substr]
        
        if verbose: print("    length aln_str:      "+str(len(aln_str)), file=sys.stderr)
        if verbose: print("    length query_str:    "+str(len(query_str)), file=sys.stderr)
        if verbose: print("    length subject_str:  "+str(len(subject_str)), file=sys.stderr)
        if verbose: print("    aln_str   :          "+aln_str, file=sys.stderr)
        if verbose: print("    query_str :          "+query_str, file=sys.stderr)
        if verbose: print("    subject_str :        "+subject_str, file=sys.stderr)
        if verbose: print("    length aln_str2:     "+str(len(aln_str2)), file=sys.stderr)
        if verbose: print("    length query_str2:   "+str(len(query_str2)), file=sys.stderr)
        if verbose: print("    length subject_str2: "+str(len(subject_str2)), file=sys.stderr)
        if verbose: print("    aln_str2   :         "+aln_str2, file=sys.stderr)
        if verbose: print("    query_str2 :         "+query_str2, file=sys.stderr)
        if verbose: print("    subject_str2 :       "+subject_str2, file=sys.stderr)
        if verbose: print("------------------------------------------", file=sys.stderr)

        # When looping through hsp, we have to keep the aln string part that matches to the actual hsp
        query_char = list(query_str2)
        subject_char = list(subject_str2)
        aln_char = list(aln_str2)
        length = len(aln_char)

        # Compute score according (more or less) fo Fahlgren 2009 (and psRNATarget).
        # Careful here because in psRNATarget, what they refer to as G:U pairs actually corresponds
        # to . characters in the alignment strings given by SSEARCH36. Furthermore, G:U pairs (or . pairs) 
        # are NOT multiplied by the penalty error multiplier (default 1.5) in the seed region (2-13 nt), only mismatches are multipled.
        # The following routine gives the exact same penalty scoring results
        # to the ones given by psRNATarget (which they call E or Expect value).
        z = 1
        score = 0
        gap_status_q = "closed"
        gap_status_s = "closed"
        gap_status = "closed"
        seed_mismatches = 0
        total_mismatches = 0
        gaps = 0
        GUs = 0
        gap_penalty = 0
        total_extra_penalty_gap_query = 0
        if verbose: print("    "+str(q_start), file=sys.stderr)
        for el in aln_char:
            el2 = query_char.pop(0)
            el3 = subject_char.pop(0)

            is_GU = "no"
            curr_score = 0
            mismatch_penalty = 0
            is_mismatch = "no"
            
            if el == " " and el2 != "-" and el3 != "-":
                curr_score = 1
                is_mismatch = "yes"
                total_mismatches+=1
            elif el == " " and el2 == "-" and gap_status_q == "closed": #gap opening
                curr_score = 2
                gap_status_q = "open"
                is_mismatch = "yes"
                gaps+=1
                total_mismatches+=1
                total_extra_penalty_gap_query = total_extra_penalty_gap_query + extra_penalty_query_gap
            elif el == " " and el3 == "-" and gap_status_s == "closed": #gap opening
                curr_score = 2
                gap_status_s = "open"
                is_mismatch = "yes"
                gaps+=1
                total_mismatches+=1
            elif el == " " and el2 == "-" and gap_status_q == "open": #gap extension
                curr_score = 0.5
                is_mismatch = "yes"
                gaps+=1
                total_mismatches+=1
                total_extra_penalty_gap_query = total_extra_penalty_gap_query + extra_penalty_query_gap
            elif el == " " and el3 == "-" and gap_status_s == "open": #gap extension
                curr_score = 0.5
                is_mismatch = "yes"
                gaps+=1
                total_mismatches+=1
            elif el == ".":
                curr_score = 0.5
                GUs+=1
                is_GU = "yes"
                is_mismatch = "no"
                gap_status_q = "closed"
                gap_status_s = "closed"
            elif el == ":":
                gap_status_q = "closed"
                gap_status_s = "closed"
            else:
                print("    Did not match anything: ["+el+"] ["+el2+"] ["+el3+"] "+gap_status_s+" "+gap_status_q, file=sys.stderr)
                sys.exit()
            
            if z >= seed_start and z <= seed_end:
                if is_mismatch == "yes":
                    if verbose: print("    is mismatch in the seed region - will apply penalty... curr_score: "+str(curr_score)+" * 1.5 => "+str(curr_score*penalty_multiplier), file=sys.stderr)
                    curr_score = curr_score * penalty_multiplier
                    seed_mismatches+=1
                else:
                    curr_score = curr_score * 1.0
                
            score = score + curr_score
            if verbose: print("    "+str(z)+":"+str(curr_score)+"    "+str(score)+"    "+el2+" "+el+" "+el3, file=sys.stderr)
            z+=1
        
        score = score + total_extra_penalty_gap_query
        penalty_score = score

        if( penalty_score <= E_cutoff and
            seed_mismatches <= num_mismatch_seed and
            hsp >= hsp_cutoff and
            gaps <= gap_cutoff and
            total_mismatches <= total_mismatches_cutoff and
            len(subject_str) <= maximum_alignment_length):
            
            if not keep_target_suffix:
                subject = subject.replace("_\d+$", "")
           
            if GUs <= GUs_cutoff:
                print(query+"\t"+subject+"\t"+aln_str+"\t"+query_str+"\t"+subject_str+"\t"+str(q_start)+"\t"+str(q_end)+"\t"+str(start)+"\t"+str(end)+"\t"+str(penalty_score)+"\t"+strand, file=sys.stdout)
            else:
                stats.get('GU_reject', 0) + 1
            
        else:
            stats.get('other_reject', 0) + 1
            if outfile_failed:
                print(query+"\t"+subject+"\t"+aln_str+"\t"+query_str+"\t"+subject_str+"\t"+str(q_start)+"\t"+str(q_end)+"\t"+str(start)+"\t"+str(end)+"\t"+str(penalty_score)+"\t"+strand, file=fhand_out_failed)
            
    if outfile_failed: fhand_out_failed.close()
    fhand.close();

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

