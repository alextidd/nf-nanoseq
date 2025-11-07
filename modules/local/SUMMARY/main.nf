process SUMMARY {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(var_csvs)

  output:
  path("summary.txt")

  script:
  """
  summary.R "./"
  """
  stub:
  """
  touch summary.txt
  """
}