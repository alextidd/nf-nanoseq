#!/usr/bin/env python3

import sys
import os
import json
import time
import argparse
from multiprocessing import Pool
from nanoseq_utils import file_chk, check_dependencies, runCommand
from nanoseq_ginterval import GInterval

def runVariantCaller(normal, duplex, ref, args_str, interval, out):
    if normal is None:
        return
    
    cmd = "variantcaller -a %s -b %s -r %s" % (normal, duplex, ref)
    cmd += " %s -i %s -o %s" % (args_str, interval, out)
    
    runCommand(cmd)
    runCommand("bgzip -f %s" % out)
    runCommand("tabix -f -s 1 -b 2 -e 2 %s.gz" % out)

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq variant caller')
    
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
    
    # Var specific arguments
    parser_opt.add_argument('-a', type=int, action='store',
                            default=50, help="minimum AS-XS (50)")
    parser_opt.add_argument('-b', type=int, action='store', default=0,
                            help="minimum matched normal reads per strand (0)")
    parser_opt.add_argument('-c', type=float, action='store',
                            default=0.02, help="fraction of clips (0.02)")
    parser_opt.add_argument('-d', type=int, action='store',
                            default=2, help="minimum duplex depth (2)")
    parser_opt.add_argument('-f', type=float, action='store', default=0.9,
                            help="minimum fraction of reads for consensus (0.9)")
    parser_opt.add_argument('-i', type=float, action='store', default=1.0,
                            help="maximum fraction of reads with an indel (1.0)")
    parser_opt.add_argument('-m', type=int, action='store',
                            default=8, help="minimum cycle number (8)")
    parser_opt.add_argument('-n', type=int, action='store',
                            default=3, help="maximum number of mismatches (3)")
    parser_opt.add_argument('-p', type=int, action='store', default=0,
                            help="minimum fraction of reads that are proper-pairs (0)")
    parser_opt.add_argument('-q', type=int, action='store',
                            default=60, help="minimum consensus base quality (60)")
    parser_opt.add_argument('-r', type=int, action='store',
                            default=144, help="read length (after 5' trimming) (144)")
    parser_opt.add_argument('-v', type=float, action='store',
                            default=0.01, help="maximum normal VAF (0.01)")
    parser_opt.add_argument('-x', type=int, action='store',
                            default=8, help="maximum cycle number (8)")
    parser_opt.add_argument('-z', type=int, action='store',
                            default=15, help="minimum normal coverage (15)")
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies(["variantcaller", "bgzip", "tabix"])
    
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
    
    if args.index is not None:
        time.sleep(args.index)
    
    tmpDir = args.out + "/tmpNanoSeq"
    
    if args.index is None or args.index == 1:
        with open("%s/var/args.json" % (tmpDir), "w") as jsonOut:
            json.dump(args.__dict__, jsonOut)
    
    print("Starting variant calling\n")
    
    # Read partition information
    with open("%s/part/nfiles" % (tmpDir), "r") as iofile:
        nfiles = int(iofile.read())
    
    # Build variant caller argument string
    args_str = "-a %s -b %s -c %s -d %s -f %s -i %s -m %s -n %s -p %s -q %s -r %s -v %s -x %s -z %s" % (
        args.a, args.b, args.c, args.d, args.f, args.i, args.m, args.n,
        args.p, args.q, args.r, args.v, args.x, args.z)
    
    # Prepare inputs for each partition
    inputs = []
    for ii in range(1, nfiles + 1):
        if os.path.isfile("%s/var/%s.vcf.gz" % (tmpDir, ii)):
            inputs.append((None, None, None, None, None, None))
        else:
            interval_file = "%s/part/%s.bed" % (tmpDir, ii)
            intervals = []
            with open(interval_file, 'r') as iofile:
                for iline in iofile:
                    parts = iline.strip().split('\t')
                    intervals.append(GInterval(parts[0], int(parts[1]), int(parts[2])))
            interval_str = ' '.join([str(x) for x in intervals])
            
            inputs.append((args.normal, args.duplex, args.ref, args_str,
                          interval_str, "%s/var/%s.vcf" % (tmpDir, ii)))
    
    if args.index is None:
        with Pool(args.threads) as p:
            p.starmap(runVariantCaller, inputs)
    else:
        commands = [[] for i in range(args.max_index)]
        jj = 0
        for ii in range(len(inputs)):
            if ii % args.max_index == 0:
                jj = 0
            commands[jj].append(
                "runVariantCaller(inputs[%s][0], inputs[%s][1], inputs[%s][2], inputs[%s][3], inputs[%s][4], inputs[%s][5])" % 
                (ii, ii, ii, ii, ii, ii))
            jj += 1
        for icmd in commands[args.index - 1]:
            exec(icmd)
    
    print("Completed variant calling\n")

if __name__ == '__main__':
    main()