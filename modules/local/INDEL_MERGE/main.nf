process INDEL_MERGE {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), val(partitions),
        path(indel_filtered_vcfs), path(indel_filtered_tbis)

  output:
  tuple val(meta),
        path("results.indel.vcf.gz"), path("results.indel.vcf.gz.tbi")

  script:
  """
  # modules
  module load bcftools-1.19/python-3.11.6

  # initialise merged vcfs
  cp indel_1.filtered.vcf.gz merged.vcf.gz

  # merge indels
  echo -e "${partitions.join('\\n')}" |
  sort -n | sed 1d |
  while read part_i ; do
    bcftools concat \\
      --no-version -Oz \\
      -o tmp.vcf.gz \\
      merged.vcf.gz \\
      indel_\${part_i}.filtered.vcf.gz
    mv tmp.vcf.gz merged.vcf.gz
  done

  # sort and index
  bcftools \
    sort -Oz \
    -o results.indel.vcf.gz \
    merged.vcf.gz
  bcftools index -t -f results.indel.vcf.gz
  """
  stub:
  """
  touch results.indel.vcf.gz
  touch results.indel.vcf.gz.tbi
  """
}