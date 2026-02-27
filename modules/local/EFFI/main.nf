process EFFI {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta),
        path(bundled_bam), path(bundled_bai),
        path(dedup_bam), path(dedup_bai)
  tuple path(fasta), path(fai)
  path effi_panel_bed

  output:
  tuple val(meta), path("efficiency/${meta.id}*")

  script:
  def effi_panel_bed_opt = effi_panel_bed.name != 'NO_FILE' ? "-panel $effi_panel_bed" : ''
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
    -out efficiency/${meta.id} \\
    $effi_panel_bed_opt
  """
  stub:
  """
  mkdir -p efficiency
  touch efficiency/${meta.id}.RBs
  touch efficiency/${meta.id}.RBs.pdf
  """
}