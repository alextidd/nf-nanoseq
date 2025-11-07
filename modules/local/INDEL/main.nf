process INDEL {
  tag "${part_i}_${meta.id}"

  input:
  tuple val(meta), val(part_i),
        path(dsa_bed),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta), val(part_i),
        path("indel_${part_i}.filtered.vcf.gz"),
        path("indel_${part_i}.filtered.vcf.gz.tbi"),
        emit: indel_merge

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0
  module add bcftools-1.19/python-3.11.6

  # step 1
  indelCaller_step1.pl \\
    -out indel_${part_i}.bed.gz \\
    -reads-bundle ${params.indel_rb} \\
    -trim3 ${params.indel_t3} \\
    -trim5 ${params.indel_t5} \\
    -min-coverage ${params.indel_z} \\
    -max-vaf ${params.indel_v} \\
    -min-as-xs ${params.indel_a} \\
    -max-clip ${params.indel_c} \\
    $dsa_bed

  # step 2
  indelCaller_step2.pl \\
    -sort \\
    -out indel_${part_i} \\
    -ref $fasta \\
    -bam $duplex_bam \\
    indel_${part_i}.bed.gz

  # step 3
  indelCaller_step3.R \\
    $fasta \\
    indel_${part_i}.vcf.gz \\
    $normal_bam \\
    ${params.indel_v}
  """
  stub:
  """
  touch indel_${part_i}.bed.gz
  touch indel_${part_i}.vcf.gz
  touch indel_${part_i}.filtered.vcf.gz
  touch indel_${part_i}.filtered.vcf.gz.tbi
  """
}