#!/usr/bin/env python3

import sys
import os
import json
import pickle
import argparse
import math
import gzip
import copy
from nanoseq_ginterval import GInterval

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq partition: partition intervals into n jobs using coverage information')
    
    parser_opt = parser.add_argument_group('optional arguments')
    parser_req = parser.add_argument_group('required arguments')
    
    # Required arguments
    parser_req.add_argument('-n', '--jobs', type=int, action='store',
                            required=True, help='partition dsa,var,indel to this many tasks')
    parser_req.add_argument('--chrs', type=str, action='store',
                            required=True, help='comma-separated list of contigs to process')
    parser_req.add_argument('--gintervals', type=str, action='store',
                            required=True, help='gIntervals.dat file from nanoseq_intervals.py')
    
    # Optional arguments
    parser_opt.add_argument('--excludeBED', action='store',
                            help='BED (gz) file with regions to exclude from analysis.')
    parser_opt.add_argument('--excludeCov', type=int, action='store',
                            help='Exclude regions with coverage values higher than this')
    parser_opt.add_argument('--win', type=int, action='store',
                            default=100, help='bin size used for coverage calculation (100)')
    
    args = parser.parse_args()
    
    print("Starting partitioning\n")
    
    # Load genomic intervals
    with open(args.gintervals, "rb") as iofile:
        gIntervals = pickle.load(iofile)
    
    # Parse chromosome list
    chr_list = args.chrs.split(',')
    print(f"Processing {len(chr_list)} chromosomes: {chr_list}\n")
    
    # Read coverage information from chromosome-named files
    coverage = []
    cctotal = 0
    chrOffset = {}
    tmpIntervals = []
    
    print("Parsing coverage files\n")
    for chr_name in chr_list:
        cov_file = f"cov_{chr_name}.bed.gz"
        if not os.path.exists(cov_file):
            sys.exit(f"Error: Coverage file {cov_file} not found")
        
        with gzip.open(cov_file, 'rt') as iofile:
            for iline in iofile:
                ichr = str(iline.split('\t')[0])
                ib = int(iline.split('\t')[1])
                ie = int(iline.split('\t')[2])
                cc = int(iline.split('\t')[3])
                cctotal += cc
                
                if args.excludeCov is not None:
                    if cc >= args.excludeCov:
                        tmpIntervals.append(GInterval(ichr, ib+1, ie))
                
                if ib == 0:
                    chrOffset[str(ichr)] = len(coverage)
                coverage.append([ib, cc])
    
    print("Completed parsing coverage files\n")
    
    # Remove regions to exclude from BED file
    if args.excludeBED is not None:
        if not os.path.isfile(args.excludeBED):
            sys.exit(f"excludeBED file {args.excludeBED} not found!")
        if not os.path.isfile(args.excludeBED + ".tbi"):
            sys.exit(f"excludeBED index file {args.excludeBED}.tbi not found!")
        
        with gzip.open(args.excludeBED, 'rt') as iofile:
            for iline in iofile:
                if iline.startswith('#'):
                    continue
                ichr = str(iline.split('\t')[0])
                ib = int(iline.split('\t')[1])
                ie = int(iline.split('\t')[2])
                tmpIntervals.append(GInterval(ichr, ib+1, ie))
        tmpIntervals.sort()
    
    if len(tmpIntervals) > 0:
        # Merge overlapping intervals
        xIntervals = [tmpIntervals.pop(0)]
        while len(tmpIntervals) > 0:
            xIntervals.extend(xIntervals.pop() + tmpIntervals.pop(0))
        
        print(f"Excluding {len(xIntervals)} intervals\n")
        
        # Remove the excluded intervals
        iiresult = []
        for ii in gIntervals:
            ifrag = ii
            for jj in xIntervals:
                if ifrag.chr != jj.chr:
                    continue
                diff = ifrag - jj
                if len(diff) == 2:
                    iiresult.append(diff[0])
                    ifrag = diff[1]
                elif len(diff) == 1:
                    ifrag = diff[0]
                else:
                    break
            else:
                iiresult.append(ifrag)
        gIntervals = iiresult
        
        # Recalculate total coverage after exclusions
        xSumCov = 0
        for iinterval in gIntervals:
            ichar = iinterval.chr
            ibeg = iinterval.beg - 1
            iend = iinterval.end - 1
            for i in range(math.floor(ibeg/args.win), math.floor(iend/args.win) + 1):
                if not ichar in chrOffset:
                    break
                j = i + chrOffset[ichar]
                xSumCov += coverage[j][1]
        cctotal = xSumCov
    
    # Determine genomic intervals for each job
    # Partition so each job has roughly the same amount of coverage
    basesPerCPU = 0
    njobs = args.jobs
    basesPerCPU = cctotal / njobs
    
    print(f"Partitioning {njobs} jobs with {basesPerCPU:.0f} bases/task\n")
    
    sumCov = 0
    oIntervals = []
    intervalsPerCPU = []
    gIntervalsCopy = copy.deepcopy(gIntervals)
    
    while len(gIntervals) > 0:
        iinterval = gIntervals.pop(0)
        ichar = iinterval.chr
        ibeg = iinterval.beg - 1
        iend = iinterval.end - 1
        
        for i in range(math.floor(ibeg/args.win), math.floor(iend/args.win) + 1):
            j = i + chrOffset[ichar]
            sumCov += coverage[j][1]
            
            if sumCov > basesPerCPU:
                jend = min([coverage[j][0] + args.win, iend])
                oIntervals.append(GInterval(ichar, ibeg + 1, jend + 1))
                intervalsPerCPU.append(oIntervals)
                oIntervals = []
                sumCov = 0
                ibeg = jend + 1
        
        if iend >= ibeg:
            oIntervals.append(GInterval(ichar, ibeg + 1, iend + 1))
    
    if len(oIntervals) > 0:
        intervalsPerCPU.append(oIntervals)
    
    while len(intervalsPerCPU) < njobs:
        intervalsPerCPU.append([])
    
    # Verify partitioning is correct
    print("Checking partition of intervals..", end='')
    flatInt = [item for sublist in intervalsPerCPU[0:njobs] for item in sublist]
    nbases1 = sum(i.l for i in flatInt)
    nbases2 = sum(i.l for i in gIntervalsCopy)
    
    if nbases1 != nbases2:
        sys.exit(f"Internal check failed: interval length after partition doesn't match original ({nbases2}, {nbases1})\n")
    
    print("..", end='')
    mIntervals = [flatInt.pop(0)]
    while len(flatInt) > 0:
        mIntervals.extend(mIntervals.pop() + flatInt.pop(0))
    
    for (i, ival) in enumerate(gIntervalsCopy):
        if not (ival == mIntervals[i]):
            sys.exit(f"Internal check failed: mismatch after part (interval {mIntervals[i]} should be {ival})\n")
    
    print(" OK\n")
    
    # Write partitions as separate BED files
    for idx, partition in enumerate(intervalsPerCPU, 1):
        with open(f"part_{idx}.bed", "w") as iofile:
            for iint in partition:
                iofile.write(f"{iint.chr}\t{iint.beg}\t{iint.end}\n")
    
    # Write number of partitions
    with open("nfiles", "w") as iofile:
        iofile.write(str(len(intervalsPerCPU)))
    
    print(f"Completed partitioning into {len(intervalsPerCPU)} jobs\n")
    
    # Print partition summary
    print("Partition summary:")
    for idx, partition in enumerate(intervalsPerCPU, 1):
        total_bases = sum(iint.l for iint in partition)
        n_intervals = len(partition)
        if n_intervals > 0:
            chrs = set(iint.chr for iint in partition)
            print(f"  Partition {idx:3d}: {n_intervals:5d} intervals | {total_bases:12,} bases | chrs: {','.join(sorted(chrs))}")
        else:
            print(f"  Partition {idx:3d}: empty")

if __name__ == '__main__':
    main()