process VAR_MERGE_COV {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), val(partitions), path(vars_covs)

  output:
  tuple val(meta), path("results.cov.bed.gz"), path("results.cov.bed.gz.tbi")

  script:
  """
  # modules
  module load bcftools-1.19/python-3.11.6

  # concat all var_cov files
  > results.cov.bed
  for var_cov in ${vars_covs.join(' ')} ; do
    bgzip -dc \$var_cov >> results_unsorted.cov.bed
  done

  # sort bed
  sort -k1,1 -k2,2n results_unsorted.cov.bed > results.cov.bed

  # index and tabix
  bgzip -@ 2 -f results.cov.bed
  sleep 3
  bgzip -@ 2 -t results.cov.bed.gz
  tabix -f results.cov.bed.gz
  """
}