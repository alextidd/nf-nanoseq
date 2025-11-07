process PREPROCESS {
  tag "${meta.id}"

  input:
  tuple val(meta), path(duplex_bam), path(duplex_bai)
  
  output:
  tuple val(meta),
        path("${meta.id}_bundled.bam"),
        path("${meta.id}_bundled.bam.bai"),
        path("${meta.id}_dedup.bam"),
        path("${meta.id}_dedup.bam.bai"),
        emit: effi
  tuple val(meta),
        path("${meta.id}_bundled.bam"),
        path("${meta.id}_bundled.bam.bai"),
        emit: cov

  script:
  """
  # modules
  module load biobambam2/2.0.180  
  module load samtools-1.19.2/python-3.11.6
  module load bcftools-1.19/python-3.11.6

  # add read bundles
  bamaddreadbundles -I $duplex_bam -O ${meta.id}_bundled.bam
  samtools index ${meta.id}_bundled.bam

  # deduplicate
  randomreadinbundle -I ${meta.id}_bundled.bam -O ${meta.id}_dedup.bam
  samtools index ${meta.id}_dedup.bam
  """
  stub:
  """
  touch ${meta.id}_bundled.bam
  touch ${meta.id}_bundled.bam.bai
  touch ${meta.id}_dedup.bam
  touch ${meta.id}_dedup.bam.bai
  """
}