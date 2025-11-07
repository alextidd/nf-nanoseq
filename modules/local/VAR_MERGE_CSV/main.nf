process VAR_MERGE_CSV {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), val(partitions), path(var_tsvs)

  output:
  tuple val(meta), path("variants.csv"), emit: var_vcf
  tuple val(meta), path("*.csv"), emit: plot
  tuple val(meta), path("*.csv"), emit: summary


  script:
  """
  var_merge_csv.py \\
    --partitions ${partitions.join(',')}
  """
}