#!/usr/bin/env python3

import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        description='NanoSeq merge vars: merge var outputs from multiple partitions into CSV files')
    parser.add_argument('--partitions', type=str, required=True,
                        help='Comma-delimited list of partition indices (e.g., "1,2,3,4,5")')
    args = parser.parse_args()
    
    # Parse partition indices
    partition_indices = [int(i.strip()) for i in args.partitions.split(',')]
    
    print(f"\nGenerating CSV files from var output files")
    print(f"\nProcessing {len(partition_indices)} partitions\n")
    
    # Define CSV files to generate
    csvFiles = ['Coverage', 'CallVsQpos', 'PyrVsMask',
                'ReadBundles', 'Burdens', 'Variants', 'DiscardedVariants', 'Mismatches']
    csvIO = {}
    
    # Open output files
    for ifile in csvFiles:
        csvIO[ifile] = open(f'{ifile.lower()}.csv', 'w')
    
    # Write headers
    csvIO['Coverage'].write('count\n')
    csvIO['CallVsQpos'].write('base,qpos,ismasked,count\n')
    csvIO['PyrVsMask'].write('pyrcontext,ismasked,count\n')
    csvIO['ReadBundles'].write('fwd,rev,ismasked,isvariant,count\n')
    csvIO['Burdens'].write('ismasked,isvariant,count\n')
    csvIO['Variants'].write('chrom,chromStart,context,commonSNP,'
                            'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
                            'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
                            'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
                            'bundleType,dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
                            'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
                            'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
                            'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
                            'dplxfwdTotal,dplxrevTotal,left,right,qpos,call,isvariant,pyrcontext,'
                            'stdcontext,pyrsub,stdsub,ismasked,dplxBarcode\n')
    csvIO['DiscardedVariants'].write('chrom,chromStart,context,commonSNP,'
                            'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
                            'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
                            'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
                            'bundleType,dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
                            'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
                            'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
                            'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
                            'dplxfwdTotal,dplxrevTotal,left,right,qpos,call,isvariant,pyrcontext,'
                            'stdcontext,pyrsub,stdsub,ismasked,dplxBarcode,'
                            'dplx_clip_filter,alignment_score_filter,mismatch_filter,matched_normal_filter,'
                            'duplex_filter,consensus_base_quality_filter,indel_filter,five_prime_trim_filter,'
                            'three_prime_trim_filter,proper_pair_filter,vaf_filter\n')
    csvIO['Mismatches'].write('chrom,chromStart,context,commonSNP,'
                              'shearwater,bulkASXS,bulkNM,bulkForwardA,bulkForwardC,bulkForwardG,'
                              'bulkForwardT,bulkForwardIndel,bulkReverseA,bulkReverseC,bulkReverseG,'
                              'bulkReverseT,bulkReverseIndel,dplxBreakpointBeg,dplxBreakpointEnd,'
                              'dplxASXS,dplxCLIP,dplxNM,dplxfwdA,dplxfwdC,dplxfwdG,dplxfwdT,'
                              'dplxfwdIndel,dplxrevA,dplxrevC,dplxrevG,dplxrevT,dplxrevIndel,'
                              'dplxCQfwdA,dplxCQfwdC,dplxCQfwdG,dplxCQfwdT,dplxCQrevA,'
                              'dplxCQrevC,dplxCQrevG,dplxCQrevT,bulkForwardTotal,bulkReverseTotal,'
                              'dplxfwdTotal,dplxrevTotal,left,right,qpos,mismatch,ismasked,dplxBarcode\n')
    
    # Process each var file
    variants_processed = 0
    files_processed = 0
    
    for i in partition_indices:
        ifile = f"var_{i}.tsv"
        
        if not os.path.isfile(ifile):
            sys.exit(f"Error: File {ifile} not found, exiting...")
        
        print(f"Processing {ifile}...")
        files_processed += 1
        
        with open(ifile, 'r') as infile:
            for row in infile:
                if row[0] == '#':
                    continue
                arow = row.strip().split('\t')
                if csvIO.get(arow[0], None):
                    csvIO[arow[0]].write('%s\n' % ','.join(arow[1:]))
                    if arow[0] == 'Variants':
                        variants_processed += 1
    
    # Close all files
    for ifile in csvIO.values():
        ifile.close()
    
    print(f"\nCompleted!")

if __name__ == '__main__':
    main()