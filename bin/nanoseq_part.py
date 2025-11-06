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
    parser_req.add_argument('-n', '--jobs', type=int, action='store',
                            required=True, help='partition dsa,var,indel to this many tasks')
    parser_req.add_argument('--chrs', type=str, action='store',
                            required=True, help='comma-separated list of contigs to process')
    parser_opt.add_argument('--excludeBED', action='store',
                            help='BED (gz) file with regions to exclude from analysis.')
    parser_opt.add_argument('--excludeCov', type=int, action='store',
                            help='Exclude regions with coverage values higher than this')
    parser_req.add_argument('--gintervals', type=str, action='store',
                            required=True, help='gIntervals.dat file from nanoseq_intervals.py')
    
    args = parser.parse_args()
    
    print("Starting partitioning\n")
    
    # Read intervals
    with open(args.gintervals, "rb") as iofile:
        gintervals = pickle.load(iofile)
    
    # Parse chromosome list
    chr_list = args.chrs.split(',')
    print(f"Processing {len(chr_list)} chromosomes: {chr_list}\n")
    
    # Build expected coverage file names
    cov_files = []
    for chr in chr_list:
        cov_file = f"{chr}.cov.bed.gz"
        if not os.path.exists(cov_file):
            sys.exit(f"Error: Coverage file {cov_file} not found")
        cov_files.append(cov_file)
    
    print(f"Found all {len(cov_files)} coverage files\n")
    
    # Read coverage information
    bychr = {}
    bychrAll = {}
    for cov_file in cov_files:
        with gzip.open(cov_file, 'rt') as iofile:
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
    
    print(f"Total coverage: {totalcov:,} bases")
    print(f"Target coverage per job: {targetcov:,.0f} bases\n")
    
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
    
    # Write partitions to bed file
    with open("partitions.bed", "w") as iofile:
        for idx, partition in enumerate(partitions, 1):
            for iint in partition:
                iofile.write("%s\t%s\t%s\t%s\n" % (iint.chr, iint.beg, iint.end, idx))

    # Write number of partitions
    with open("nfiles", "w") as iofile:
        iofile.write(str(len(partitions)))
    
    # Save partitions as pickle file
    with open("partitions.dat", "wb") as iofile:
        pickle.dump(partitions, iofile)
    
    print("Completed partitioning into %s jobs\n" % len(partitions))

if __name__ == '__main__':
    main()