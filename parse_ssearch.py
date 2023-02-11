#!/usr/bin/env python

"""parse_ssearch.py. From the default output of ssearch36, this script parses each alignment to a .tsv format that can then be used in the companion 
script parse_mirna_targets.py
python>=3.9.0

NOTES:
It is not entirely clear which ssearch36 parameters are used by psRNATarget, but using the following SSEARCH followed by the execution of this script with default parame
ssearch36 -f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000 -i <input_mirna.fasta> <reference_genome.fasta> > <output_file>
If you provide an already rev-complemented fasta file, you can omit the -i argument.

Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
"""

import argparse
import os
import sys
import re
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def parse_command_line_arguments():
    parser = argparse.ArgumentParser(description='Convert SSEARCH36 default output to a .tsv format output')
    parser.add_argument('-i', '--infile', required=True, help='Input file (i.e. output of ssearch36)', type=argparse.FileType('r'))
    parser.add_argument('--rev', default=False, action=argparse.BooleanOptionalAction, help='Reverse')
    parser.add_argument('--verbose', default=False, action=argparse.BooleanOptionalAction, help='Verbose output')
    args = parser.parse_args()
    return args

def main(arguments):
    args = parse_command_line_arguments()
    infile = os.path.abspath(args.infile.name)
    verbose = args.verbose
    rev = args.rev
    
    hash = {}
    curr_query = "";
    curr_target = "";
    curr_aln = "";
    curr_query_length = int;
    counter_match_nucl_string = 0;
    curr_start = int;
    curr_end = int;
    start = int;
    end = int;
    hsp = int;
    q_end = int;
    q_start = int;
    query_str = str;
    aln_str = str;
    subject_str = str;
    i = 0;
    strand = str;

    fhand = open(infile, 'r') if infile else sys.stdin

    for line in fhand:
        line = line.rstrip()

        if line.startswith("#"):
            continue

        match = re.match(r"\s+\d+>>>(\S+) - (\d+) nt", line)
        if match:
            if verbose: print("\n---------------------------------------------", file=sys.stderr)
            counter_match_nucl_string = 0
            curr_start = 0
            curr_end = 0
            curr_query = match.group(1)
            curr_query_length = match.group(2)
            i = 0

            if verbose: print("curr_query: "+curr_query, file=sys.stderr)
            continue
        
        match = re.match(r"^>>(\S+) .*\((\d+) nt\)$", line)
        if match:
            i+=1
            curr_target = match.group(1)+"_"+str(i)
            continue
        
        #Smith-Waterman score: 279 89.5% identity (100.0% similar) in 19 nt overlap (2-20:52164-52182)
        match = re.match(r"^Smith-Waterman.*\((\d+)-(\d+):(\d+)-(\d+)\)$", line)
        if match:
            diff_start = int(match.group(1)) - 1
            diff_end = int(curr_query_length)
            start = int(match.group(3)) - diff_start
            end = int(match.group(4)) + diff_end
            hsp = int(match.group(4)) - int(match.group(3))
            if rev is True:
                q_start = int(match.group(1))
                q_end = int(match.group(2))
                strand = "-"
            else:
                # If fwd, qstart is in reverse orientation.
                q_start = int(match.group(2))
                q_end = int(match.group(1))
                strand = "+"
            
            continue
        
        match = re.match(r"^>--$", line) #same contig, but in another location
        if match:
            i+=1
            continue
        
        match = re.match(r"^\S+\s+([ACGTU-]*)\s*$", line)
        if match and counter_match_nucl_string == 0:
            # Found a query nucl alignment string.
            curr_start = match.span(1)[0]
            curr_end = match.span(1)[1]
            curr_str = line[curr_start:curr_end]
            if curr_query not in hash:
                hash[curr_query] = {}
                hash[curr_query][curr_target + "_" + str(i)] = 1
            else:
                hash[curr_query][curr_target + "_" + str(i)] = 1

            query_str = curr_str
            
            counter_match_nucl_string = 1
            continue
        
        match = re.match(r"^\s+[\.\:]", line)
        if match:
            if counter_match_nucl_string == 1:
                # Found a match string.
                #extract substring at previously found positions.
                curr_str = line[curr_start:curr_end]
                aln_str = curr_str
                continue
        
        match = re.match(r"^\S+\s+[ACGTBDHUYNVWRMKS-]*\s*$", line)
        if match and counter_match_nucl_string == 1:
            #extract substring at previously found positions.
            curr_str = line[curr_start:curr_end]
            # Found subject string of current match.
            subject_str = curr_str
           
            if curr_target + "_" + str(i) not in hash[curr_query]:
                raise Exception("Problem at " + curr_query + "\n" + curr_target + "_"+str(i)+"\n" + line + "\n" + "line number " +str(sys._getframe().f_lineno) + " in the alignment file")

            counter_match_nucl_string = 0

            print(str(curr_query) + "\t" + str(curr_target) + "\t" + str(aln_str) + "\t" + str(query_str) + "\t" + str(subject_str) + "\t" + str(q_start) + "\t" + str(q_end) + "\t" + str(start) + "\t" + str(end) + "\t" + str(strand) + "\t" + str(hsp), file=sys.stdout)
            continue
    
    fhand.close()
    #print STDERR  Dumper(\%hash) if($verbose)
    #print STDERR "Done parsing ssearch standard format outfmt file.\n"


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

