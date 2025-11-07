#!/usr/bin/env python3

import sys
import os
import re
import json
import pickle
import argparse
from nanoseq_utils import file_chk
from nanoseq_ginterval import GInterval

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq: create genomic intervals')
    
    parser_opt = parser.add_argument_group('optional arguments')
    parser_req = parser.add_argument_group('required arguments')
    
    # Arguments
    parser_req.add_argument('-R', '--ref', action='store',
                            required=True, help="reference sequence")
    
    # Interval filtering arguments
    parser_opt.add_argument('--exclude', action='store', default='MT,GL%%,NC_%,hs37d5',
                            help='List of contigs to exclude. Comma separated, %% acts as a wild card.')
    parser_opt.add_argument('--include', action='store',
                            help='Only include these contigs. Comma separated, %% acts as a wild card.')
    parser_opt.add_argument('--larger', type=int, action='store',
                            default=0, help='Only include contigs larger than this size. (0)')
    
    args = parser.parse_args()
    
    # Build chromosome dictionary
    if args.exclude is None or args.exclude == "":
        excludes = []
    else:
        excludes = [re.compile(istr + "$")
                    for istr in args.exclude.replace("%", ".+").split(',')]
    if args.include is None or args.include == "":
        includes = [re.compile(".+")]
    else:
        includes = [re.compile(istr + "$")
                    for istr in args.include.replace("%", ".+").split(',')]
    
    chrList = []
    rnames = {}
    with open(args.ref + '.fai', 'r') as iofile:
        for iline in iofile:
            ichr = iline.split('\t')[0]
            if any(iregx.match(ichr) for iregx in includes):
                if any(iregx.match(ichr) for iregx in excludes):
                    continue
                ilength = int(iline.split('\t')[1])
                if ilength <= args.larger:
                    continue
                chrList.append(ichr)
                rnames[ichr] = ilength
    
    # Create genomic intervals
    gintervals = []
    for ichr in chrList:
        gintervals.append(GInterval(ichr, 2, rnames[ichr]-1))
    gintervals.sort()
    
    print("\nAnalysing contigs: %s\n" % [g.chr for g in gintervals])
    print(f"Total intervals: {len(gintervals)}")
    print(f"Total bases: {sum(g.l for g in gintervals):,}\n")
    
    # Save intervals as pickle
    with open("gIntervals.dat", 'wb') as iofile:
        pickle.dump(gintervals, iofile)
    
    # Save intervals as BED
    with open("intervals.bed", 'w') as bedfile:
        for g in gintervals:
            bedfile.write(f"{g.chr}\t{g.beg}\t{g.end}\n")

if __name__ == '__main__':
    main()