process CHECK_CONTIGS {
  tag "${meta.id}"

  input:
tuple val(meta), 
      path(duplex_bam, stageAs: "duplex.bam"), 
      path(duplex_bai, stageAs: "duplex.bam.bai"),
      path(normal_bam, stageAs: "normal.bam"), 
      path(normal_bai, stageAs: "normal.bam.bai")
  tuple path(fasta), path(fai)
  
  output:
  tuple val(meta), 
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)

  script:
  def check_normal = params.dedup_as_normal ? "" : """
  # validate normal BAM is not corrupted
  samtools quickcheck -v $normal_bam

  # extract contig info from normal BAM
  samtools idxstats $normal_bam | cut -f1,2 | grep -v "^\\*" > normal.contigs

  # check that normal and duplex BAMs have matching contigs
  if ! diff normal.contigs duplex.contigs > /dev/null; then
    echo "ERROR: Normal and duplex BAMs have different contigs or lengths"
    diff normal.contigs duplex.contigs
    exit 1
  fi
  """
  """
  # modules
  module load samtools-1.19.2/python-3.11.6

  # validate duplex BAM is not corrupted
  samtools quickcheck -v $duplex_bam

  # extract contig info from duplex BAM
  samtools idxstats $duplex_bam | cut -f1,2 | grep -v "^\\*" > duplex.contigs

  # extract contig info from reference
  cut -f1,2 $fai > ref.contigs

  $check_normal

  # check that the BAM contigs match the reference contigs
  if ! diff duplex.contigs ref.contigs > /dev/null; then
    echo "ERROR: BAM contigs do not match reference contigs"
    diff duplex.contigs ref.contigs
    exit 1
  fi

  echo "BAM validation passed for ${meta.id}"
  """
  stub:
  """
  echo "Stub: BAM validation skipped"
  """
}