// check bam is not truncated before proceeding + index
process samtools_index {
  tag "${meta.id}"

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path(bam), path("${bam}.bai")

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck ${bam}
  samtools index -@ ${task.cpus} ${bam}
  """
  stub:
  """
  touch ${bam}.bai
  """
}