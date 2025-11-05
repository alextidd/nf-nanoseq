#!/usr/bin/env python3

import sys
import os
import json
import time
import argparse
from multiprocessing import Pool
from nanoseq_utils import file_chk, check_dependencies, runCommand, vcfHeader
from nanoseq_ginterval import GInterval

def runIndelCaller(normal, duplex, ref, sample, rb, t3, t5, interval, out, tmpDir, idx):
    if normal is None:
        return
    
    # Step 1: Extract indel candidates
    step1_out = "%s/indel/step1_%s.txt" % (tmpDir, idx)
    cmd1 = "indelCaller_step1.pl -r %s -b %s -rb %s -t5 %s -t3 %s -i '%s' > %s" % (
        ref, duplex, rb, t5, t3, interval, step1_out)
    runCommand(cmd1)
    
    # Step 2: Get normal coverage
    step2_out = "%s/indel/step2_%s.txt" % (tmpDir, idx)
    cmd2 = "indelCaller_step2.pl -r %s -n %s -i %s > %s" % (
        ref, normal, step1_out, step2_out)
    runCommand(cmd2)
    
    # Step 3: Filter and format
    cmd3 = "indelCaller_step3.R %s %s %s" % (step2_out, sample, out)
    runCommand(cmd3)
    
    runCommand("bgzip -f %s" % out)
    runCommand("tabix -f -s 1 -b 2 -e 2 %s.gz" % out)

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq indel caller')
    
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
    
    # Indel specific arguments
    parser_opt.add_argument('-s', '--sample', action='store',
                            default='sample_1', help="sample name in output vcf (sample_1)")
    parser_opt.add_argument('--rb', type=int, action='store',
                            default=2, help="minimum reads in a bundle. (2)")
    parser_opt.add_argument('--t3', type=int, action='store', default=136,
                            help="excess bases above this value are trimmed from 3' (136)")
    parser_opt.add_argument('--t5', type=int, action='store',
                            default=8, help="bases to trim from 5' reads (8)")
    parser_opt.add_argument('-a', type=int, action='store',
                            default=50, help="minimum AS-XS (50)")
    parser_opt.add_argument('-c', type=float, action='store',
                            default=0.02, help="fraction of clips (0.02)")
    parser_opt.add_argument('-z', type=int, action='store', default=15,
                            help="minimum normal coverage (mc) (15)")
    parser_opt.add_argument('-v', type=float, action='store', default=0.01,
                            help="maximum normal VAF (0.01)")
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies(["indelCaller_step1.pl", "indelCaller_step2.pl", 
                       "indelCaller_step3.R", "Rscript", "bgzip", "tabix"])
    
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
        with open("%s/indel/args.json" % (tmpDir), "w") as jsonOut:
            json.dump(args.__dict__, jsonOut)
        
        # Write VCF header
        header = vcfHeader(args)
        header = header.replace('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
                               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % args.sample)
        with open("%s/indel/header.txt" % (tmpDir), "w") as iofile:
            iofile.write(header)
    
    print("Starting indel calling\n")
    
    # Read partition information
    with open("%s/part/nfiles" % (tmpDir), "r") as iofile:
        nfiles = int(iofile.read())
    
    # Prepare inputs for each partition
    inputs = []
    for ii in range(1, nfiles + 1):
        if os.path.isfile("%s/indel/%s.indel.vcf.gz" % (tmpDir, ii)):
            inputs.append((None, None, None, None, None, None, None, None, None, None, None))
        else:
            interval_file = "%s/part/%s.bed" % (tmpDir, ii)
            intervals = []
            with open(interval_file, 'r') as iofile:
                for iline in iofile:
                    parts = iline.strip().split('\t')
                    intervals.append(GInterval(parts[0], int(parts[1]), int(parts[2])))
            interval_str = ' '.join([str(x) for x in intervals])
            
            inputs.append((args.normal, args.duplex, args.ref, args.sample,
                          args.rb, args.t3, args.t5, interval_str,
                          "%s/indel/%s.indel.vcf" % (tmpDir, ii), tmpDir, ii))
    
    if args.index is None:
        with Pool(args.threads) as p:
            p.starmap(runIndelCaller, inputs)
    else:
        commands = [[] for i in range(args.max_index)]
        jj = 0
        for ii in range(len(inputs)):
            if ii % args.max_index == 0:
                jj = 0
            commands[jj].append(
                "runIndelCaller(inputs[%s][0], inputs[%s][1], inputs[%s][2], inputs[%s][3], inputs[%s][4], inputs[%s][5], inputs[%s][6], inputs[%s][7], inputs[%s][8], inputs[%s][9], inputs[%s][10])" % 
                (ii, ii, ii, ii, ii, ii, ii, ii, ii, ii, ii))
            jj += 1
        for icmd in commands[args.index - 1]:
            exec(icmd)
    
    print("Completed indel calling\n")

if __name__ == '__main__':
    main()