process COV {
  tag "${chr}_${meta.id}"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai),
        val(chr)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: dsa
  tuple val(meta), val(chr), path("cov_${chr}.bed.gz"), emit: part

  script:
  """
  # modules
  module load bcftools-1.19/python-3.11.6

  # get coverage per window
  bamcov \\
    --min-MQ ${params.cov_q} \\
    --region ${chr} \\
    --win ${params.cov_window} \\
    --output cov_${chr}.bed \\
    $duplex_bam

  # bgzip output
  bgzip --compress-level 2 --force cov_${chr}.bed
  """
  stub:
  """
  touch cov_${chr}.bed.gz
  """
}