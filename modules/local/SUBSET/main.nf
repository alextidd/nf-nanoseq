process SUBSET {
  tag "${meta.id}"

  input:
  tuple val(meta), path(bam), path(bai)
  path(contigs)

  output:
  tuple val(meta),
        path("${meta.id}_subset.bam"), path("${meta.id}_subset.bam.bai")

  script:
  """
  # modules
  module load samtools-1.19.2/python-3.11.6

  # filter bam to included contigs only
  samtools view --bam $bam \$(cat $contigs) --output ${meta.id}_subset.bam
  samtools index ${meta.id}_subset.bam
  """
  stub:
  """
  touch ${meta.id}_subset.bam
  touch ${meta.id}_subset.bam.bai
  """
}