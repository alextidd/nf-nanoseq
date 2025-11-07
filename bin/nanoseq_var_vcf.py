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
    """Parse variants.csv file and return list of variant dictionaries"""
    variants = []
    
    with open(csv_file, 'r') as f:
        header = f.readline().strip().split(',')
        
        for line in f:
            if not line.strip():
                continue
            
            fields = line.strip().split(',')
            if len(fields) < len(header):
                continue
            
            variant = dict(zip(header, fields))
            variants.append(variant)
    
    return variants

def variant_to_vcf_line(variant):
    """Convert variant dictionary to VCF line"""
    chrom = variant.get('chr', '.')
    pos = variant.get('pos', '.')
    ref = variant.get('ref', '.')
    alt = variant.get('mut', '.')
    
    # Build INFO field
    info_parts = []
    
    if 'btag' in variant and variant['btag']:
        info_parts.append(f"BTAG={variant['btag']}")
    if 'bbeg' in variant and variant['bbeg']:
        info_parts.append(f"BBEG={variant['bbeg']}")
    if 'bend' in variant and variant['bend']:
        info_parts.append(f"BEND={variant['bend']}")
    if 'tri' in variant and variant['tri']:
        info_parts.append(f"TRI={variant['tri']}")
    if 'qpos' in variant and variant['qpos']:
        info_parts.append(f"QPOS={variant['qpos']}")
    if 'dplx_fwd_total' in variant and variant['dplx_fwd_total']:
        info_parts.append(f"DEPTH_FWD={variant['dplx_fwd_total']}")
    if 'dplx_rev_total' in variant and variant['dplx_rev_total']:
        info_parts.append(f"DEPTH_REV={variant['dplx_rev_total']}")
    if 'bulk_fwd_total' in variant and variant['bulk_fwd_total']:
        info_parts.append(f"DEPTH_NORM_FWD={variant['bulk_fwd_total']}")
    if 'bulk_rev_total' in variant and variant['bulk_rev_total']:
        info_parts.append(f"DEPTH_NORM_REV={variant['bulk_rev_total']}")
    if 'dplx_asxs' in variant and variant['dplx_asxs']:
        info_parts.append(f"DPLX_ASXS={variant['dplx_asxs']}")
    if 'dplx_clip' in variant and variant['dplx_clip']:
        info_parts.append(f"DPLX_CLIP={variant['dplx_clip']}")
    if 'dplx_nm' in variant and variant['dplx_nm']:
        info_parts.append(f"DPLX_NM={variant['dplx_nm']}")
    if 'bulk_asxs' in variant and variant['bulk_asxs']:
        info_parts.append(f"BULK_ASXS={variant['bulk_asxs']}")
    if 'bulk_nm' in variant and variant['bulk_nm']:
        info_parts.append(f"BULK_NM={variant['bulk_nm']}")
    
    info = ';'.join(info_parts) if info_parts else '.'
    
    # Determine filter status
    filter_status = variant.get('filter', 'PASS')
    if filter_status == '':
        filter_status = 'PASS'
    
    # VCF line
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\n"

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
    
    # Parse variants
    variants = parse_variants_csv(args.variants_csv)
    
    print(f"Found {len(variants)} variants")
    
    # Sort variants by chromosome and position
    variants.sort(key=lambda v: (chr_sort_key(v.get('chr', '')), 
                                  int(v.get('pos', 0)) if v.get('pos', '0').isdigit() else 0))
    
    print(f"Writing VCF to {args.output}...")
    
    # Write VCF
    with gzip.open(args.output, 'wt') as vcf:
        # Write header
        vcf.write(vcfHeader(args.ref_fai, args.sample))
        
        # Write variants
        for variant in variants:
            vcf.write(variant_to_vcf_line(variant))
    
    print(f"Successfully wrote {len(variants)} variants to {args.output}")

if __name__ == '__main__':
    main()