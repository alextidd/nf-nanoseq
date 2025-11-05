#!/usr/bin/env python3

import sys
import os
import json
import time
import argparse
from multiprocessing import Pool
from nanoseq_utils import file_chk, check_dependencies, runCommand
from nanoseq_ginterval import GInterval

def runDSA(normal, duplex, ref, snp, mask, args_str, interval, out):
    if normal is None:
        return
    
    cmd = "dsa -a %s -b %s -r %s" % (normal, duplex, ref)
    if snp is not None:
        cmd += " -c %s" % snp
    if mask is not None:
        cmd += " -D %s" % mask
    cmd += " %s -i %s -o %s" % (args_str, interval, out)
    
    runCommand(cmd)
    runCommand("bgzip -f %s" % out)
    runCommand("tabix -f -s 1 -b 2 -e 2 %s.gz" % out)

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq DSA: compute tables')
    
    parser_opt = parser.add_argument_group('optional arguments')
    parser_req = parser.add_argument_group('required arguments')
    
    # Arguments
    parser_opt.add_argument('--out', action='store', default='.',
                            help='path of the output files and scratch directory (.)')
    parser_opt.add_argument('-j', '--index', type=int, action='store',
                            help='index of the LSF job array. One based')
    parser_opt.add_argument('-k', '--max_index', type=int,
                            action='store', help='maximum index of the LSF job array')
    parser_opt.add_argument('-t', '--threads', type=int,
                            action='store', default=1, help='number of threads (1)')
    parser_req.add_argument('-R', '--ref', action='store',
                            required=True, help="reference sequence")
    parser_req.add_argument('-A', '--normal', action='store',
                            required=True, help="normal BAM / CRAM")
    parser_req.add_argument('-B', '--duplex', action='store',
                            required=True, help="duplex (tumour) BAM / CRAM")
    
    # DSA specific arguments
    parser_opt.add_argument('-C', '--snp', action='store',
                            help="SNP BED (gz) file")
    parser_opt.add_argument('-D', '--mask', action='store',
                            help="mask BED (gz) file")
    parser_opt.add_argument('-d', type=int, action='store',
                            default=2, help="minimum duplex depth (2)")
    parser_opt.add_argument('-q', type=int, action='store',
                            default=30, help="minimum base quality for normal (30)")
    parser_opt.add_argument('-M', type=int, action='store',
                            help="minimum map quality")
    parser_opt.add_argument('--no_test', action='store_true',
                            help="skip BAM format tests, use with caution")
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies(["dsa", "bgzip", "tabix"])
    
    # Validate arguments
    if ((args.index is None and args.max_index is not None) or
            (args.index is not None and args.max_index is None)):
        sys.exit("Must specify index and max_index for array execution!")
    
    jobArray = args.index is not None and args.max_index is not None
    
    if jobArray and args.threads > 1:
        print("Warning: execution with a job array is single threaded so thread argument is ignored")
    
    # Check files
    ext = os.path.splitext(args.duplex)[1][0:-1] + "i"
    file_chk(args.duplex, ext, "BAM/CRAM", sys.exit)
    ext = os.path.splitext(args.normal)[1][0:-1] + "i"
    file_chk(args.normal, ext, "BAM/CRAM", sys.exit)
    file_chk(args.ref, ".fai", "Reference", sys.exit)
    
    if args.snp is not None:
        file_chk(args.snp, ".tbi", "SNP BED", sys.exit)
    if args.mask is not None:
        file_chk(args.mask, ".tbi", "Mask BED", sys.exit)
    
    if args.index is not None:
        time.sleep(args.index)
    
    tmpDir = args.out + "/tmpNanoSeq"
    
    if args.index is None or args.index == 1:
        with open("%s/dsa/args.json" % (tmpDir), "w") as jsonOut:
            json.dump(args.__dict__, jsonOut)
    
    print("Starting dsa calculation\n")
    
    # Read partition information
    with open("%s/part/nfiles" % (tmpDir), "r") as iofile:
        nfiles = int(iofile.read())
    
    # Build DSA argument string
    args_str = "-d %s -q %s" % (args.d, args.q)
    if hasattr(args, 'M') and args.M is not None:
        args_str += " -M %s" % args.M
    if args.no_test:
        args_str += " --no-test"
    
    # Prepare inputs for each partition
    inputs = []
    for ii in range(1, nfiles + 1):
        if os.path.isfile("%s/dsa/%s.dsa.bed.gz" % (tmpDir, ii)):
            inputs.append((None, None, None, None, None, None, None, None))
        else:
            interval_file = "%s/part/%s.bed" % (tmpDir, ii)
            # Convert BED to interval string
            intervals = []
            with open(interval_file, 'r') as iofile:
                for iline in iofile:
                    parts = iline.strip().split('\t')
                    intervals.append(GInterval(parts[0], int(parts[1]), int(parts[2])).convert2DSAInput())
            interval_str = ' '.join([str(x) for x in intervals])
            
            inputs.append((args.normal, args.duplex, args.ref, args.snp, args.mask,
                          args_str, interval_str, "%s/dsa/%s.dsa.bed" % (tmpDir, ii)))
    
    if args.index is None:
        with Pool(args.threads) as p:
            p.starmap(runDSA, inputs)
    else:
        commands = [[] for i in range(args.max_index)]
        jj = 0
        for ii in range(len(inputs)):
            if ii % args.max_index == 0:
                jj = 0
            commands[jj].append(
                "runDSA(inputs[%s][0], inputs[%s][1], inputs[%s][2], inputs[%s][3], inputs[%s][4], inputs[%s][5], inputs[%s][6], inputs[%s][7])" % 
                (ii, ii, ii, ii, ii, ii, ii, ii))
            jj += 1
        for icmd in commands[args.index - 1]:
            exec(icmd)
    
    print("Completed dsa calculation\n")

if __name__ == '__main__':
    main()