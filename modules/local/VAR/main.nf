process VAR {
  tag "${part_i}_${meta.id}"

  input:
  tuple val(meta), val(part_i), path(dsa_bed)

  output:
  tuple val(meta), val(part_i), path("var_${part_i}.tsv"), emit: var_vcf
  tuple val(meta), val(part_i), path("var_${part_i}.tsv"), emit: var_merge_csv
  tuple val(meta), val(part_i), path("var_cov_${part_i}.bed.gz"), emit: var_merge_cov

  script:
  """
  # run variantcaller
  variantcaller \\
    -B $dsa_bed \\
    -U var_cov_${part_i}.bed \\
    -O var_${part_i}.tsv \\
    -D var_discarded_${part_i}.tsv \\
    -a ${params.var_a} \\
    -b ${params.var_b} \\
    -c ${params.var_c} \\
    -d ${params.var_d} \\
    -f ${params.var_f} \\
    -i ${params.var_i} \\
    -m ${params.var_m} \\
    -n ${params.var_n} \\
    -p ${params.var_p} \\
    -q ${params.var_q} \\
    -r ${params.var_r} \\
    -v ${params.var_v} \\
    -x ${params.var_x} \\
    -z ${params.var_z}
  """
  stub:
  """
  touch var_${part_i}.tsv
  touch var_discarded_${part_i}.tsv
  touch var_cov_${part_i}.bed.gz
  """
}