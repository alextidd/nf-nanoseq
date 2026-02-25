#!/usr/bin/env python3
# filepath: /nfs/casm/team268im/at31/nextflow/nf-nanoseq/bin/generate_muts_vcf.py

import argparse
import gzip
import os
import sys

def vcfHeader(ref_fai, sample='sample_1'):
    """Generate VCF header with contig information from reference .fai file"""
    header = '##fileformat=VCFv4.2\n'
    header += '##source=NanoSeq pipeline\n'
    header += '##FILTER=<ID=PASS,Description="All filters passed">\n'
    
    # Add contig information from .fai file
    if ref_fai and os.path.isfile(ref_fai):
        contigs = []
        with open(ref_fai, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    chrom = fields[0]
                    length = fields[1]
                    contigs.append((chrom, length))
        
        # Sort contigs (chromosomes first, then others)
        contigs.sort(key=lambda x: chr_sort_key(x[0]))
        
        for chrom, length in contigs:
            header += f"##contig=<ID={chrom},length={length}>\n"
    
    # INFO fields
    header += '##INFO=<ID=BTAG,Number=1,Type=String,Description="Read bundle tag (duplex barcode)">\n'
    header += '##INFO=<ID=BBEG,Number=1,Type=String,Description="Read bundle left breakpoint">\n'
    header += '##INFO=<ID=BEND,Number=1,Type=String,Description="Read bundle right breakpoint">\n'
    header += '##INFO=<ID=TRI,Number=1,Type=String,Description="Pyrimidine context, trinucleotide substitution">\n'
    header += '##INFO=<ID=QPOS,Number=1,Type=Integer,Description="Read position closest to 5-prime end">\n'
    header += '##INFO=<ID=DEPTH_FWD,Number=1,Type=Integer,Description="Read bundle forward reads depth">\n'
    header += '##INFO=<ID=DEPTH_REV,Number=1,Type=Integer,Description="Read bundle reverse reads depth">\n'
    header += '##INFO=<ID=DEPTH_NORM_FWD,Number=1,Type=Integer,Description="Matched normal forward reads depth">\n'
    header += '##INFO=<ID=DEPTH_NORM_REV,Number=1,Type=Integer,Description="Matched normal reverse reads depth">\n'
    header += '##INFO=<ID=DPLX_ASXS,Number=1,Type=Integer,Description="AS-XS for duplex">\n'
    header += '##INFO=<ID=DPLX_CLIP,Number=1,Type=Integer,Description="Clipping for duplex">\n'
    header += '##INFO=<ID=DPLX_NM,Number=1,Type=Integer,Description="Mismatches in duplex">\n'
    header += '##INFO=<ID=BULK_ASXS,Number=1,Type=Integer,Description="AS-XS for bulk">\n'
    header += '##INFO=<ID=BULK_NM,Number=1,Type=Integer,Description="Mismatches in bulk">\n'
    
    # FILTER fields
    header += '##FILTER=<ID=dbsnp,Description="Common SNP site">\n'
    header += '##FILTER=<ID=shearwater,Description="Noisy site">\n'
    
    # Column headers
    header += f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    
    return header

def chr_sort_key(chrom):
    """Generate sort key for chromosome names"""
    import re
    
    # Remove 'chr' prefix for comparison
    chrom_clean = re.sub('^chr', '', chrom)
    
    # Numeric chromosomes come first
    if chrom_clean.isdigit():
        return (0, int(chrom_clean))
    # Then X, Y
    elif chrom_clean == 'X':
        return (1, 0)
    elif chrom_clean == 'Y':
        return (1, 1)
    # Then M/MT
    elif chrom_clean in ('M', 'MT'):
        return (2, 0)
    # Everything else
    else:
        return (3, chrom_clean)

def parse_variants_csv(csv_file):
    """Parse variants.csv file and return list of variant dictionaries

    NOTE: old pipeline (runNanoSeq.py) writes CSVs with header names like
    `chrom,chromStart,context,...,call,...`. To mirror that behaviour we
    provide a simple parser but the main conversion below uses the same
    code-path as runNanoSeq: reading the CSV into a dict-of-lists and
    writing VCF rows by index.
    """
    variants = []
    with open(csv_file, 'r') as f:
        header = f.readline().rstrip('\n').split(',')
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split(',')
            # pad shorter rows so zip won't miss columns
            if len(fields) < len(header):
                fields += [''] * (len(header) - len(fields))
            variant = dict(zip(header, fields))
            variants.append(variant)
    return variants

def variant_to_vcf_line_from_var_dict(var, i):
    """Build VCF line using the same fields/logic as runNanoSeq.py.

    `var` is a dict-of-lists produced by reading `variants.csv` and `i`
    is the index of the variant to write.
    """
    def _get(key):
        return var.get(key, [''])[i] if key in var else ''

    chrom = _get('chrom') or '.'
    chromStart = _get('chromStart') or ''
    try:
        pos = int(chromStart) + 1 if chromStart != '' else '.'
    except ValueError:
        pos = '.'

    # REF comes from middle base of context, ALT comes from call
    context = _get('context')
    ref = context[1] if context and len(context) > 1 else '.'
    alt = _get('call') or '.'

    # determine filter
    ifilter = 'PASS'
    if _get('shearwater') == '1':
        ifilter = 'shearwater'
    if _get('commonSNP') == '1':
        ifilter = 'dbsnp'

    # INFO fields (same order as runNanoSeq)
    info = 'BTAG=%s;BBEG=%s;BEND=%s;TRI=%s;QPOS=%s;DEPTH_FWD=%s;DEPTH_REV=%s;DEPTH_NORM_FWD=%s;DEPTH_NORM_REV=%s;DPLX_ASXS=%s;DPLX_CLIP=%s;DPLX_NM=%s;BULK_ASXS=%s;BULK_NM=%s' % (
        _get('dplxBarcode'), _get('dplxBreakpointBeg'), _get('dplxBreakpointEnd'), _get('pyrsub'), _get('qpos'),
        _get('dplxfwdTotal'), _get('dplxrevTotal'), _get('bulkForwardTotal'), _get('bulkReverseTotal'),
        _get('dplxASXS'), _get('dplxCLIP'), _get('dplxNM'), _get('bulkASXS'), _get('bulkNM'))

    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{ifilter}\t{info}\n"

def main():
    parser = argparse.ArgumentParser(
        description='Generate muts.vcf.gz file from variants.csv')
    
    parser.add_argument('--variants_csv', type=str, required=True,
                        help='Path to variants.csv file')
    parser.add_argument('--ref_fai', type=str,
                        help='Path to reference .fai file (optional, for contig headers)')
    parser.add_argument('--sample', type=str, default='sample_1',
                        help='Sample name for VCF header (default: sample_1)')
    parser.add_argument('--output', type=str, default='results.muts.vcf.gz',
                        help='Output VCF file (default: results.muts.vcf.gz)')
    
    args = parser.parse_args()
    
    print(f"Reading variants from {args.variants_csv}...")

    # Read csv into dict-of-lists (same approach as runNanoSeq.py)
    var = {}
    nVariants = 0
    with open(args.variants_csv, 'r') as iofile:
        iline = iofile.readline().rstrip('\n')
        fields = iline.split(',')
        for ifield in fields:
            var[ifield] = []
        for iline in iofile:
            if not iline.rstrip('\n'):
                continue
            nVariants += 1
            parts = iline.rstrip('\n').split(',')
            # pad to header length
            if len(parts) < len(fields):
                parts += [''] * (len(fields) - len(parts))
            for (i, ival) in enumerate(parts):
                var[fields[i]].append(ival)

    # added so code will not break with older files
    if ('dplxBarcode' not in var) or len(var.get('dplxBarcode', [])) == 0:
        var['dplxBarcode'] = [''] * nVariants

    print(f"Found {nVariants} variants")
    print(f"Writing VCF to {args.output}...")

    # write VCF header and body
    with gzip.open(args.output, 'wt') as vcf:
        vcf.write(vcfHeader(args.ref_fai, args.sample))
        for i in range(nVariants):
            vcf.write(variant_to_vcf_line_from_var_dict(var, i))

    print(f"Successfully wrote {nVariants} variants to {args.output}")

if __name__ == '__main__':
    main()