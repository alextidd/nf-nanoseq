process CHECK_CONTIGS {
  tag "${meta.id}"

  input:
  tuple val(meta), 
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)
  
  output:
  tuple val(meta), 
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)

  script:
  """
  # modules
  module load samtools-1.19.2/python-3.11.6

  # validate BAM files are not corrupted
  samtools quickcheck -v $duplex_bam
  samtools quickcheck -v $normal_bam

  # extract contig info from BAMs
  samtools idxstats $normal_bam | cut -f1,2 | grep -v "^\\*" > normal.contigs
  samtools idxstats $duplex_bam | cut -f1,2 | grep -v "^\\*" > duplex.contigs

  # extract contig info from reference
  cut -f1,2 $fai > ref.contigs

  # check that normal and duplex BAMs have matching contigs
  if ! diff normal.contigs duplex.contigs > /dev/null; then
    echo "ERROR: Normal and duplex BAMs have different contigs or lengths"
    diff normal.contigs duplex.contigs
    exit 1
  fi

  # check that the BAM contigs match the reference contigs
  if ! diff normal.contigs ref.contigs > /dev/null; then
    echo "ERROR: BAM contigs do not match reference contigs"
    diff normal.contigs ref.contigs
    exit 1
  fi

  echo "BAM validation passed for ${meta.id}"
  """
  stub:
  """
  echo "Stub: BAM validation skipped"
  """
}