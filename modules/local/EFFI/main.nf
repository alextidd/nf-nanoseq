process EFFI {
  tag "${meta.id}"

  input:
  tuple val(meta),
        path(bundled_bam), path(bundled_bai),
        path(dedup_bam), path(dedup_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta), path("efficiency/${meta.id}*")

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.38.0 
  module load bcftools-1.9/python-3.11.6 

  # dirs
  mkdir -p efficiency

  # calculate efficiency
  # TODO: convert this to bash?
  efficiency.pl \\
    -threads $task.cpus \\
    -duplex $bundled_bam \\
    -dedup $dedup_bam \\
    -ref $fasta \\
    -out efficiency/${meta.id}
  """
  stub:
  """
  mkdir -p efficiency
  touch efficiency/${meta.id}.RBs
  touch efficiency/${meta.id}.RBs.pdf
  """
}