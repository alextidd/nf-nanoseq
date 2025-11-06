import sys
import subprocess
import tempfile
import os
import re
import shutil

def file_chk(fn, idx_ext, msg_prefix, parser):
    if not os.path.isfile(fn):
        parser.error(f"{msg_prefix} file {fn} was not found!")
    if not os.path.isfile(fn + idx_ext):
        parser.error(f"{msg_prefix} index file {fn}{idx_ext} was not found!")

def check_dependencies(scripts):
    for icode in scripts:
        if shutil.which(icode) is None:
            raise ValueError("%s was not found in path!" % icode)

def runCommand(command):
    if command is None:
        return
    for ijob in command.rstrip(';').split(';'):
        print("\nExecuting: %s\n" % ijob)
        p = subprocess.Popen(
            ijob, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        std = p.communicate()
        if p.returncode != 0:
            error = std[1].decode()
            sys.stderr.write("\n!Error processing:  %s\n" % ijob)
            raise ValueError(error)

def getBAMcontigs(bam):
    contigs = []
    if bam is None:
        return contigs
    job = "samtools view -H %s" % bam
    p = subprocess.Popen(
        job, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while True:
        output = p.stdout.readline().decode()
        if output == '' and p.poll() is not None:
            break
        if re.match("@SQ", output):
            contigs.append(output.strip())
    rc = p.poll()
    if rc != 0:
        error = p.stderr.read().decode()
        sys.stderr.write("\n!Error processing:  %s\n" % job)
        raise ValueError(error)
    return contigs

def vcfHeader(args):
    header = '##fileformat=VCFv4.2\n'
    header += '##source=NanoSeq pipeline\n'
    header += '##FILTER=<ID=PASS,Description="All filters passed">\n'
    header += "##reference=file://%s\n" % (args.ref)

    from nanoseq_ginterval import GInterval
    contigs = []
    with open(args.ref + ".fai", 'r') as iofile:
        for iline in iofile:
            ichr = iline.split('\t')[0]
            ilength = int(iline.split('\t')[1])
            contigs.append(GInterval(ichr, 1, ilength))
    contigs.sort()
    for ii in contigs:
        ichr = ii.chr
        ilength = ii.end
        header += "##contig=<ID=%s,length=%s>\n" % (ichr, ilength)
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
    header += '##FILTER=<ID=dbsnp,Description="Common SNP site">\n'
    header += '##FILTER=<ID=shearwater,Description="Noisy site">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    return header