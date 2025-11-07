#!/usr/bin/env python3

import sys
import os
import re
import json
import pickle
import argparse
from multiprocessing import Pool
from nanoseq_utils import check_dependencies, runCommand
from nanoseq_ginterval import GInterval

def runBamcov(bam, mapQ, window, ichr, out):
    if bam is None:
        return
    out = out+".cov.bed"
    runCommand("bamcov -q %s -w %s -r \"%s\" -o %s %s" %
               (mapQ, window, ichr, out, bam))
    runCommand("bgzip -l 2 -f %s" % (out))
    outdone = re.sub('cov.bed$', 'done', out)
    open(outdone, 'w').close()

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq coverage calculation')
    
    parser_opt = parser.add_argument_group('optional arguments')
    parser_req = parser.add_argument_group('required arguments')
    
    # Arguments
    parser_opt.add_argument('-t', '--threads', type=int,
                            action='store', default=1, help='number of threads (1)')
    parser_req.add_argument('-R', '--ref', action='store',
                            required=True, help="reference sequence")
    parser_req.add_argument('-A', '--normal', action='store',
                            required=True, help="normal BAM / CRAM")
    parser_req.add_argument('-B', '--duplex', action='store',
                            required=True, help="duplex (tumour) BAM / CRAM")
    parser_req.add_argument('--gintervals', action='store',
                            required=True, help="gIntervals.dat file from nanoseq_intervals.py")
    
    # Cov specific arguments
    parser_opt.add_argument('-w', '--win', type=int, action='store',
                            default=100, help='bin size for coverage distribution (100)')
    parser_opt.add_argument('-Q', type=int, action='store',
                            default=0, help="minimum mapQ to include a duplex read (0)")
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies(["bamcov", "bgzip"])
    
    # Check gIntervals file exists
    if not os.path.isfile(args.gintervals):
        sys.exit(f"gIntervals file {args.gintervals} not found!")
    
    # Create output directory
    os.makedirs('cov', exist_ok=True)
    
    with open("cov/args.json", "w") as jsonOut:
        json.dump(args.__dict__, jsonOut)
    
    # Load gIntervals from file
    print(f"Loading intervals from {args.gintervals}\n")
    with open(args.gintervals, 'rb') as iofile:
        gintervals = pickle.load(iofile)
    
    print(f"Loaded {len(gintervals)} intervals\n")
    
    # Extract chromosome list from gIntervals
    chrList = [g.chr for g in gintervals]
    
    print("\nAnalysing contigs: %s\n" % chrList)
    print("Starting cov calculation\n")
    
    with open("cov/nfiles", "w") as iofile:
        iofile.write(str(len(chrList)))
    
    inputs = []
    for ii, ichr in enumerate(chrList):
        if os.path.isfile("cov/%s.done" % (ii + 1)):
            inputs.append((None, None, None, None, None))
        else:
            inputs.append((args.duplex, str(args.Q), str(
                args.win), ichr, "cov/%s" % (ii + 1)))
    
    with Pool(args.threads) as p:
        p.starmap(runBamcov, inputs)
    
    print("Completed cov calculation\n")

if __name__ == '__main__':
    main()