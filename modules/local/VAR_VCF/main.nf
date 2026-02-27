process VAR_VCF {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(vars_csv)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta), path("results.muts.vcf.gz"), path("results.muts.vcf.gz.tbi")

  script:
  """
  module load bcftools-1.19/python-3.11.6
  
  var_vcf.py \\
    --variants_csv $vars_csv \\
    --ref_fai $fai \\
    --sample ${meta.id} \\
    --output results.muts.vcf.gz
  
  bcftools index -t -f results.muts.vcf.gz
  """
}