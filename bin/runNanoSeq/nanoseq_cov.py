#!/usr/bin/env python3

import sys
import os
import re
import json
import pickle
import time
import argparse
from multiprocessing import Pool
from nanoseq_utils import file_chk, check_dependencies, runCommand, getBAMcontigs
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
    
    # Cov specific arguments
    parser_opt.add_argument('--exclude', action='store', default='MT,GL%%,NC_%,hs37d5',
                            help='List of contigs to exclude. Comma separated, %% acts as a wild card.')
    parser_opt.add_argument('--include', action='store',
                            help='Only include these contigs. Comma separated, %% acts as a wild card.')
    parser_opt.add_argument('--larger', type=int, action='store',
                            default=0, help='Only include contigs larger than this size. (0)')
    parser_opt.add_argument('-w', '--win', type=int, action='store',
                            default=100, help='bin size for coverage distribution (100)')
    parser_opt.add_argument('-Q', type=int, action='store',
                            default=0, help="minimum mapQ to include a duplex read (0)")
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies(["bamcov", "bgzip"])
    
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
    for idir in ('cov', 'part', 'dsa', 'var', 'indel', 'post'):
        if not os.path.isdir(tmpDir+'/'+idir):
            os.makedirs(tmpDir+'/'+idir, exist_ok=True)
    
    if args.index is None or args.index == 1:
        with open("%s/%s/args.json" % (tmpDir, 'cov'), "w") as jsonOut:
            json.dump(args.__dict__, jsonOut)
    
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
    
    # BAM sanity checks
    bamNContigs = {}
    bamNOrder = []
    for iline in getBAMcontigs(args.normal):
        ichr = iline.split('\t')[1].replace('SN:', '')
        ilength = int(iline.split('\t')[2].replace('LN:', ''))
        bamNContigs[ichr] = ilength
        bamNOrder.append(ichr)
    
    bamTContigs = {}
    bamTOrder = []
    for iline in getBAMcontigs(args.duplex):
        ichr = iline.split('\t')[1].replace('SN:', '')
        ilength = int(iline.split('\t')[2].replace('LN:', ''))
        bamTContigs[ichr] = ilength
        bamTOrder.append(ichr)
    
    for icontig in rnames:
        if not (icontig in bamNContigs):
            sys.exit("Reference contig %s was not found in normal BAM" % icontig)
        if bamNContigs[icontig] != rnames[icontig]:
            sys.exit("Length of contig %s in normal BAM doesn't match reference (%s vs %s)" %
                     (icontig, bamNContigs[icontig], rnames[icontig]))
    
    for icontig in rnames:
        if not (icontig in bamTContigs):
            sys.exit("Reference contig %s was not found in duplex BAM" % icontig)
        if bamTContigs[icontig] != rnames[icontig]:
            sys.exit("Length of contig %s in duplex BAM doesn't match reference (%s vs %s)" %
                     (icontig, bamTContigs[icontig], rnames[icontig]))
    
    for icontig in rnames:
        if bamNOrder.index(icontig) != bamTOrder.index(icontig):
            sys.exit("Contigs in BAM files must have the same order ( check order in headers )")
    
    gintervals = []
    for ichr in chrList:
        gintervals.append(GInterval(ichr, 2, rnames[ichr]-1))
    gintervals.sort()
    
    reorderchr = []
    for iint in gintervals:
        reorderchr.append(iint.chr)
    chrList = reorderchr
    print("\nAnalysing contigs: %s\n" % chrList)
    print("Starting cov calculation\n")
    
    if args.index is None or args.index == 1:
        with open("%s/cov/nfiles" % (tmpDir), "w") as iofile:
            iofile.write(str(len(chrList)))
    
    inputs = []
    for ii, ichr in enumerate(chrList):
        if os.path.isfile("%s/cov/%s.done" % (tmpDir, ii + 1)):
            inputs.append((None, None, None, None, None))
        else:
            inputs.append((args.duplex, str(args.Q), str(
                args.win), ichr, "%s/cov/%s" % (tmpDir, ii + 1)))
    
    if args.index is None:
        with open("%s/cov/%s" % (tmpDir, 'gIntervals.dat'), 'wb') as iofile:
            pickle.dump(gintervals, iofile)
        with Pool(args.threads) as p:
            p.starmap(runBamcov, inputs)
    else:
        if args.index == 1:
            with open("%s/cov/%s" % (tmpDir, 'gIntervals.dat'), 'wb') as iofile:
                pickle.dump(gintervals, iofile)
        
        commands = [[] for i in range(args.max_index)]
        jj = 0
        for ii in range(len(inputs)):
            if ii % args.max_index == 0:
                jj = 0
            commands[jj].append(
                "runBamcov(inputs[%s][0], inputs[%s][1], inputs[%s][2], inputs[%s][3], inputs[%s][4])" % (ii, ii, ii, ii, ii))
            jj += 1
        for icmd in commands[args.index - 1]:
            exec(icmd)
    
    print("Completed cov calculation\n")

if __name__ == '__main__':
    main()