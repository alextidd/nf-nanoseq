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

    file="part_${part_i}_\${chr}_\${start}_\${end}.bed"
    echo \$file

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
      -O \$file ;
    
    # check number of fields for truncation - should be 45
    zcat \$file.gz |
    awk -v file="\$file.gz" '
    /^#/ { next }   # skip header/comment lines

    {
        data_lines++
        last_nf = NF
    }

    END {
        if (data_lines == 0) {
            print "WARNING: No data lines in " file > "/dev/stderr"
            exit 0
        }

        if (last_nf < 45) {
            print "ERROR: Final data line has " last_nf \\
                  " fields (<45) in " file > "/dev/stderr"
            exit 1
        }
    }'

    # append to final output
    zcat \$file.gz >> dsa_${part_i}.bed

    # remove intermediate
    rm \$file.gz

  done

  # bgzip, test integrity, index
  bgzip -f -l 2 dsa_${part_i}.bed
  sleep 2
  bgzip -t dsa_${part_i}.bed.gz
  """
  stub:
  """
  touch dsa_${part_i}.bed.gz
  """
}