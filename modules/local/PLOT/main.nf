process PLOT {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(variants_csv)
  path post_triNuc

  output:
  path("plots/*")

  script:
  """
  mkdir -p plots
  plot.R \
    "./" \
    "./plots/results" \
    $post_triNuc
  """
  stub:
  """
  
  """
}