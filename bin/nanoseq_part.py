#!/usr/bin/env python3

import sys
import os
import json
import pickle
import gzip
import argparse
from nanoseq_utils import file_chk
from nanoseq_ginterval import GInterval

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq partition: partition intervals into n jobs using coverage information')
    
    parser_opt = parser.add_argument_group('optional arguments')
    parser_req = parser.add_argument_group('required arguments')
    
    # Arguments
    parser_opt.add_argument('--out', action='store', default='.',
                            help='path of the output files and scratch directory (.)')
    parser_req.add_argument('-R', '--ref', action='store',
                            required=True, help="reference sequence")
    
    # Part specific arguments
    parser_req.add_argument('-n', '--jobs', type=int, action='store',
                            required=True, help='partition dsa,var,indel to this many tasks')
    parser_opt.add_argument('--excludeBED', action='store',
                            help='BED (gz) file with regions to exclude from analysis.')
    parser_opt.add_argument('--excludeCov', type=int, action='store',
                            help='Exclude regions with coverage values higher than this')
    
    args = parser.parse_args()
    
    # Check files
    file_chk(args.ref, ".fai", "Reference", sys.exit)
    
    tmpDir = args.out + "/tmpNanoSeq"
    
    print("Starting partitioning\n")
    
    # Read coverage args
    with open("%s/cov/args.json" % (tmpDir), "r") as jsonIn:
        argscov = json.load(jsonIn)
    
    # Read intervals
    with open("%s/cov/gIntervals.dat" % (tmpDir), 'rb') as iofile:
        gintervals = pickle.load(iofile)
    
    # Read number of coverage files
    with open("%s/cov/nfiles" % (tmpDir), "r") as iofile:
        nfiles = int(iofile.read())
    
    # Read coverage information
    bychr = {}
    bychrAll = {}
    for iindex in range(1, nfiles + 1):
        with gzip.open("%s/cov/%s.cov.bed.gz" % (tmpDir, iindex), 'rt') as iofile:
            for iline in iofile:
                ichr = iline.split('\t')[0]
                istart = int(iline.split('\t')[1])
                iend = int(iline.split('\t')[2])
                icov = int(iline.split('\t')[3])
                if ichr not in bychrAll:
                    bychrAll[ichr] = []
                bychrAll[ichr].append(GInterval(ichr, istart, iend))
                
                if args.excludeCov is not None and icov > args.excludeCov:
                    continue
                
                if ichr not in bychr:
                    bychr[ichr] = []
                bychr[ichr].append(GInterval(ichr, istart, iend))
    
    # Read excluded regions if provided
    if args.excludeBED is not None:
        file_chk(args.excludeBED, ".tbi", "excludeBED", sys.exit)
        with gzip.open(args.excludeBED, 'rt') as iofile:
            for iline in iofile:
                if iline.startswith('#'):
                    continue
                ichr = iline.split('\t')[0]
                istart = int(iline.split('\t')[1])
                iend = int(iline.split('\t')[2])
                if ichr not in bychr:
                    continue
                
                newintervals = []
                for iint in bychr[ichr]:
                    newintervals.extend(iint - GInterval(ichr, istart, iend))
                bychr[ichr] = newintervals
    
    # Combine intervals
    allintervals = []
    for ichr in bychr:
        allintervals.extend(bychr[ichr])
    allintervals.sort()
    
    # Calculate total coverage
    totalcov = sum([iint.l for iint in allintervals])
    targetcov = totalcov / args.jobs
    
    # Partition into jobs
    partitions = []
    currentPartition = []
    currentCov = 0
    
    for iint in allintervals:
        if currentCov + iint.l <= targetcov or len(currentPartition) == 0:
            currentPartition.append(iint)
            currentCov += iint.l
        else:
            partitions.append(currentPartition)
            currentPartition = [iint]
            currentCov = iint.l
    
    if currentPartition:
        partitions.append(currentPartition)
    
    # Write partitions
    for idx, partition in enumerate(partitions, 1):
        with open("%s/part/%s.bed" % (tmpDir, idx), 'w') as iofile:
            for iint in partition:
                iofile.write("%s\t%s\t%s\n" % (iint.chr, iint.beg, iint.end))
    
    # Write number of partitions
    with open("%s/part/nfiles" % (tmpDir), "w") as iofile:
        iofile.write(str(len(partitions)))
    
    # Save args
    with open("%s/part/args.json" % (tmpDir), "w") as jsonOut:
        json.dump(args.__dict__, jsonOut)
    
    print("Completed partitioning into %s jobs\n" % len(partitions))

if __name__ == '__main__':
    main()