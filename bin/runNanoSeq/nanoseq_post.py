#!/usr/bin/env python3

import sys
import os
import json
import glob
import argparse
from nanoseq_utils import file_chk, check_dependencies, runCommand

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq post-processing: gather final files, compute summaries')
    
    parser_opt = parser.add_argument_group('optional arguments')
    parser_req = parser.add_argument_group('required arguments')
    
    # Arguments
    parser_opt.add_argument('--out', action='store', default='.',
                            help='path of the output files and scratch directory (.)')
    parser_req.add_argument('-R', '--ref', action='store',
                            required=True, help="reference sequence")
    parser_req.add_argument('-A', '--normal', action='store',
                            required=True, help="normal BAM / CRAM")
    parser_req.add_argument('-B', '--duplex', action='store',
                            required=True, help="duplex (tumour) BAM / CRAM")
    
    # Post specific arguments
    parser_opt.add_argument('--name', action='store',
                            default='results', help="name for output files (results)")
    parser_opt.add_argument('--triNuc', action='store',
                            help="tri-nucleotide correction file")
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies(["bcftools", "bgzip", "tabix", "variantcaller.R",
                       "efficiency_nanoseq.pl", "efficiency_nanoseq.R"])
    
    # Check files
    ext = os.path.splitext(args.duplex)[1][0:-1] + "i"
    file_chk(args.duplex, ext, "BAM/CRAM", sys.exit)
    ext = os.path.splitext(args.normal)[1][0:-1] + "i"
    file_chk(args.normal, ext, "BAM/CRAM", sys.exit)
    file_chk(args.ref, ".fai", "Reference", sys.exit)
    
    if args.triNuc is not None:
        if not os.path.isfile(args.triNuc):
            sys.exit("triNuc file %s was not found!" % args.triNuc)
    
    tmpDir = args.out + "/tmpNanoSeq"
    
    with open("%s/post/args.json" % (tmpDir), "w") as jsonOut:
        json.dump(args.__dict__, jsonOut)
    
    print("Starting post-processing\n")
    
    # Gather DSA files
    print("Gathering DSA files\n")
    dsa_files = sorted(glob.glob("%s/dsa/*.dsa.bed.gz" % tmpDir))
    if dsa_files:
        dsa_list = ' '.join(dsa_files)
        runCommand("bcftools concat -a %s | bgzip -c > %s/post/%s.dsa.bed.gz" % 
                  (dsa_list, tmpDir, args.name))
        runCommand("tabix -f -s 1 -b 2 -e 2 %s/post/%s.dsa.bed.gz" % (tmpDir, args.name))
    
    # Gather variant files
    print("Gathering variant files\n")
    var_files = sorted(glob.glob("%s/var/*.vcf.gz" % tmpDir))
    if var_files:
        var_list = ' '.join(var_files)
        runCommand("bcftools concat -a %s | bgzip -c > %s/post/%s.var.vcf.gz" % 
                  (var_list, tmpDir, args.name))
        runCommand("tabix -f -p vcf %s/post/%s.var.vcf.gz" % (tmpDir, args.name))
    
    # Gather indel files
    print("Gathering indel files\n")
    indel_files = sorted(glob.glob("%s/indel/*.indel.vcf.gz" % tmpDir))
    if indel_files:
        indel_list = ' '.join(indel_files)
        runCommand("bcftools concat -a %s | bgzip -c > %s/post/%s.indel.vcf.gz" % 
                  (indel_list, tmpDir, args.name))
        runCommand("tabix -f -p vcf %s/post/%s.indel.vcf.gz" % (tmpDir, args.name))
    
    # Run variantcaller.R for post-processing
    if os.path.isfile("%s/post/%s.var.vcf.gz" % (tmpDir, args.name)):
        print("Running variantcaller.R\n")
        cmd = "variantcaller.R %s/post/%s.var.vcf.gz %s/post/%s" % (
            tmpDir, args.name, tmpDir, args.name)
        if args.triNuc is not None:
            cmd += " %s" % args.triNuc
        runCommand(cmd)
    
    # Calculate efficiency metrics
    print("Calculating efficiency metrics\n")
    runCommand("efficiency_nanoseq.pl -duplex %s -dedup %s -ref %s -out %s/post/%s" % (
        args.duplex, args.normal, args.ref, tmpDir, args.name))
    
    print("Completed post-processing\n")
    print("Results are in: %s/post/\n" % tmpDir)

if __name__ == '__main__':
    main()