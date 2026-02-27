process DSA {
  tag "${part_i}_${meta.id}"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai),
        val(part_i), path(part_bed)
  tuple path(fasta), path(fai)
  tuple path(dsa_noise_bed), path(dsa_noise_tbi)
  tuple path(dsa_snp_bed), path(dsa_snp_tbi)

  output:
  tuple val(meta), val(part_i),
        path("dsa_${part_i}.bed.gz"),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai),
        emit: indel
  tuple val(meta), val(part_i),
        path("dsa_${part_i}.bed.gz"),
        emit: var

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load bcftools-1.19/python-3.11.6

  # initialise dsa file
  > dsa_${part_i}.bed

  # run dsa per partition segment
  cat $part_bed |
  while read -r chr start end ; do

    # write dsa output per region
    dsa \\
      -A $normal_bam \\
      -B $duplex_bam \\
      -C $dsa_snp_bed \\
      -D $dsa_noise_bed \\
      -R $fasta \\
      -Q ${params.dsa_q} \\
      -M ${params.dsa_M} \\
      -r "\$chr" -b \$start -e \$end \\
      -d ${params.dsa_d} \\
      -O dsa_${part_i}_\${chr}_\${start}_\${end}.bed ;
    
    # check number of fields for truncation - should be 45
    awk 'END{ if (NF != 45) print "Truncated dsa output file for region \$chr:\$start-\$end !" > "/dev/stderr"}{ if (NF != 45) exit 1 }' \
      dsa_${part_i}_\${chr}_\${start}_\${end}.bed

    # append to final output
    cat dsa_${part_i}_\${chr}_\${start}_\${end}.bed >> dsa_${part_i}.bed

    # remove intermediate
    rm dsa_${part_i}_\${chr}_\${start}_\${end}.bed

  done

  # bgzip and test integrity
  bgzip -f -l 2 dsa_${part_i}.bed
  sleep 2
  bgzip -t dsa_${part_i}.bed.gz
  """
  stub:
  """
  touch dsa_${part_i}.bed.gz
  """
}