#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { SAMTOOLS_INDEX } from './modules/local/SAMTOOLS_INDEX'

process CHECK_CONTIGS {
  tag "${meta.id}"

  input:
  tuple val(meta), 
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)
  
  output:
  tuple val(meta), 
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)

  script:
  """
  # modules
  module load samtools-1.19.2/python-3.11.6

  # validate BAM files are not corrupted
  samtools quickcheck -v $duplex_bam
  samtools quickcheck -v $normal_bam

  # extract contig info from BAMs
  samtools idxstats $normal_bam | cut -f1,2 | grep -v "^\\*" > normal.contigs
  samtools idxstats $duplex_bam | cut -f1,2 | grep -v "^\\*" > duplex.contigs

  # extract contig info from reference
  cut -f1,2 $fai > ref.contigs

  # check that normal and duplex BAMs have matching contigs
  if ! diff normal.contigs duplex.contigs > /dev/null; then
    echo "ERROR: Normal and duplex BAMs have different contigs or lengths"
    diff normal.contigs duplex.contigs
    exit 1
  fi

  # check that the BAM contigs match the reference contigs
  if ! diff normal.contigs ref.contigs > /dev/null; then
    echo "ERROR: BAM contigs do not match reference contigs"
    diff normal.contigs ref.contigs
    exit 1
  fi

  echo "BAM validation passed for ${meta.id}"
  """
  stub:
  """
  echo "Stub: BAM validation skipped"
  """
}

process PREPROCESS {
  tag "${meta.id}"

  input:
  tuple val(meta), path(duplex_bam), path(duplex_bai)
  
  output:
  tuple val(meta),
        path("${meta.id}_bundled.bam"),
        path("${meta.id}_bundled.bam.bai"),
        path("${meta.id}_dedup.bam"),
        path("${meta.id}_dedup.bam.bai"),
        emit: effi
  tuple val(meta),
        path("${meta.id}_bundled.bam"),
        path("${meta.id}_bundled.bam.bai"),
        emit: cov

  script:
  """
  # modules
  module load biobambam2/2.0.180  
  module load samtools-1.19.2/python-3.11.6
  module load bcftools-1.19/python-3.11.6

  # add read bundles
  bamaddreadbundles -I $duplex_bam -O ${meta.id}_bundled.bam
  samtools index ${meta.id}_bundled.bam

  # deduplicate
  randomreadinbundle -I ${meta.id}_bundled.bam -O ${meta.id}_dedup.bam
  samtools index ${meta.id}_dedup.bam
  """
  stub:
  """
  touch ${meta.id}_bundled.bam
  touch ${meta.id}_bundled.bam.bai
  touch ${meta.id}_dedup.bam
  touch ${meta.id}_dedup.bam.bai
  """
}

process EFFI {
  tag "${meta.id}"

  input:
  tuple val(meta),
        path(bundled_bam), path(bundled_bai),
        path(dedup_bam), path(dedup_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta), path("efficiency/${meta.id}*")

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.38.0 
  module load bcftools-1.9/python-3.11.6 

  # dirs
  mkdir -p efficiency

  # calculate efficiency
  efficiency_nanoseq.pl \
    -threads $task.cpus \
    -duplex $bundled_bam \
    -dedup $dedup_bam \
    -ref $fasta \
    -out efficiency/${meta.id}
  """
  stub:
  """
  mkdir -p efficiency
  touch efficiency/${meta.id}.RBs
  touch efficiency/${meta.id}.RBs.pdf
  """
}

process INTERVALS {
  input:
  tuple path(fasta), path(fai)

  output:
  tuple path("gIntervals.dat"), path("intervals.bed"), emit: intervals
  path "contigs.txt", emit: contigs

  script:
  def int_include = params.int_include ? "--include ${params.int_include}" : ""
  def int_exclude = params.int_exclude ? "--exclude ${params.int_exclude}" : ""
  def int_larger = params.int_larger ? "--larger ${params.int_larger}" : ""
  """
  nanoseq_intervals.py \
    --ref $fasta \
    $int_include $int_exclude $int_larger
  cut -f1 intervals.bed > contigs.txt
  """
  stub:
  """
  touch gIntervals.dat
  touch intervals.bed
  """
}

process COV {
  tag "${chr}_${meta.id}"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai),
        val(chr)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: dsa
  tuple val(meta), val(chr), path("cov_${chr}.bed.gz"), emit: part

  script:
  """
  # modules
  module load bcftools-1.19/python-3.11.6

  # get coverage per window
  bamcov \\
    --min-MQ ${params.cov_q} \\
    --region ${chr} \\
    --win ${params.cov_window} \\
    --output cov_${chr}.bed \\
    $duplex_bam

  # bgzip output
  bgzip --compress-level 2 --force cov_${chr}.bed
  """
  stub:
  """
  touch cov_${chr}.bed.gz
  """
}

process PART {
  tag "${meta.id}"

  input:
  tuple val(meta), val(chrs), path(cov_beds),
        path(gintervals), path(intervals_bed)

  output:
  tuple val(meta), path("part_*.bed", arity: params.jobs)

  script:
  def excludeBED = params.part_excludeBED ? "--excludeBED ${params.part_excludeBED}" : ""
  """
  # TODO: move excludeCov and excludeBED functionalities of nanoseq_part.py to command line
  # excludeCov is RENanoSeq-specific
  nanoseq_part.py \\
    --jobs ${params.jobs} \\
    --excludeCov ${params.part_excludeCov} \\
    --chrs ${chrs.join(',')} \\
    --gintervals $gintervals \\
    ${excludeBED}
  """
  stub:
  """
  for i in \$(seq 1 ${params.jobs}); do
    touch part_\${i}.bed
  done
  """
}

process DSA {
  tag "${part_i}_${meta.id}"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai),
        val(part_i), path(part_bed)
  tuple path(fasta), path(fai)
  tuple path(dsa_noise_bed), path(dsa_noise_tbi)

  output:
  tuple val(meta), val(part_i),
        path("dsa_${part_i}.bed.gz"), path("dsa_${part_i}.bed.gz.tbi"),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai),
        emit: indel
  tuple val(meta), val(part_i),
        path("part_${part_i}.bed.gz"), path("dsa_${part_i}.bed.gz.tbi"),
        emit: var

  script:
  // TODO: change `> dsa.bed` to `> dsa_${part_i}.bed`
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load bcftools-1.19/python-3.11.6

  # run dsa per partition segment
  > dsa_${part_i}.bed
  cat $part_bed |
  while read -r chr start end ; do
    dsa \\
      -A $normal_bam \\
      -B $duplex_bam \\
      -D $dsa_noise_bed \\
      -R $fasta \\
      -d ${params.dsa_d} \\
      -Q ${params.dsa_q} \\
      -M ${params.dsa_M} \\
      -t \\
      -r "\$chr" -b \$start -e \$end \\
      >> dsa_${part_i}.bed ;
  done

  # check number of fields for truncation - should be 45
  awk 'END{ if (NF != 45) print "Truncated dsa output file for partition $part_i !" > "/dev/stderr"}{ if (NF != 45) exit 1 }' \
    dsa_${part_i}.bed

  # bgzip and index
  bgzip -f -l 2 dsa_${part_i}.bed
  sleep 2
  bgzip -t dsa_${part_i}.bed.gz
  """
  stub:
  """
  touch dsa_${part_i}.bed.gz
  touch dsa_${part_i}.bed.gz.tbi
  """
}

process VAR {
  tag "${part_i}_${meta.id}"

  input:
  tuple val(meta), val(part_i), path(dsa_bed), path(dsa_tbi)

  output:
  tuple val(meta), val(part_i), path("var_${part_i}.tsv"), emit: merge_vars

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
  touch var_cov_${part_i}.bed
  """
}

process MERGE_VARS {
  tag "${meta.id}"

  input:
  tuple val(meta), val(partitions), path(var_tsvs)

  output:
  tuple val(meta),
        path("coverage.csv"), path("callvsqpos.csv"), path("pyrvsmask.csv"),
        path("readbundles.csv"), path("burdens.csv"), path("variants.csv"),
        path("discardedvariants.csv"), path("mismatches.csv")

  script:
  """
  nanoseq_merge_vars.py \\
    --partitions ${partitions.join(',')}
  """
}

process INDEL {
  tag "${part_i}_${meta.id}"

  input:
  tuple val(meta), val(part_i),
        path(dsa_bed), path(dsa_tbi),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path("indel_${part_i}.bed.gz"),
        path("indel_${part_i}.vcf.gz"),
        path("indel_${part_i}.filtered.vcf.gz"),
        path("indel_${part_i}.filtered.vcf.gz.tbi"),
        emit: merge_indels

  script:
  """
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

// process MERGE_INDELS
// process MERGE_VAR_COV
// process SUMMARY

process POST {
  tag "${meta.id}"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)
  path(post_triNuc)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: done
  tuple val(meta), path("tmpNanoSeq/post/*"), emit: post

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.30.0 
  module load bcftools-1.19/python-3.11.6

  # run
  runNanoSeq.py \
    --threads $task.cpus  \
    --normal $normal_bam \
    --duplex $duplex_bam \
    --ref $fasta \
    post \
    --triNuc $post_triNuc
  """
}

workflow {

  // TODO: keep id with normal_bam (can have different matched normals for the same donor)

  // get duplex bams
  ch_duplex =
    channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
            def meta = [donor_id: row.donor_id, id: row.id, type: "duplex"]
            [meta, file(row.duplex_bam, checkIfExists: true)]
    }

  // get normal bams
  ch_normal =
    channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
            def meta = [donor_id: row.donor_id, id: row.donor_id + "_normal", type: "normal"]
            [meta, file(row.normal_bam, checkIfExists: true)]
    }
    .unique()

  // get reference files
  fasta = [file(params.fasta, checkIfExists: true),
           file(params.fasta + ".fai", checkIfExists: true)]
  post_triNuc = file(params.post_triNuc, checkIfExists: true)
  dsa_noise_bed = [file(params.dsa_noise_bed, checkIfExists: true),
                   file(params.dsa_noise_bed + ".tbi", checkIfExists: true)]

  // index bams
  SAMTOOLS_INDEX(ch_duplex.concat(ch_normal))

  // split indexed channels again
  ch_bams =
    SAMTOOLS_INDEX.out
    .branch { meta, bam, bai ->
        duplex: meta.type == "duplex"
        normal: meta.type == "normal"
    }

  // check that bam and ref contigs match
  // TODO: keep meta.id with the donor, use in join
  ch_check =
    ch_bams.normal
    .map { meta, bam, bai -> [ meta.donor_id, [meta.id, bam, bai] ] }
    .join(
        ch_bams.duplex.map { meta, bam, bai -> [ meta.donor_id, [meta, bam, bai] ] }
    )
    .map { donor_id, normal, duplex ->
        def (normal_meta, normal_bam, normal_bai) = normal
        def (duplex_meta, duplex_bam, duplex_bai) = duplex
        tuple(duplex_meta, duplex_bam, duplex_bai, normal_bam, normal_bai)
    }
  CHECK_CONTIGS(ch_check, fasta)

  // get contig intervals
  INTERVALS(fasta)

  // create channel of contigs
  ch_contigs =
    INTERVALS.out.contigs
    .splitCsv()

  // preprocess duplex
  PREPROCESS(ch_bams.duplex)

  // effi
  EFFI(PREPROCESS.out.effi, fasta)

  // rejoin duplex and normal channels by donor id
  // TODO: join by id
  ch_input =
    ch_bams.normal
    .map { meta, bam, bai -> [ meta.donor_id, [meta.id, bam, bai] ] }
    .join(
        PREPROCESS.out.cov.map { meta, bam, bai -> [ meta.donor_id, [meta, bam, bai] ] }
    )
    .map { _donor_id, normal, duplex ->
        def (_normal_meta, normal_bam, normal_bai) = normal
        def (duplex_meta, duplex_bam, duplex_bai) = duplex
        tuple(duplex_meta, duplex_bam, duplex_bai, normal_bam, normal_bai)
    }
  ch_input_per_chr = ch_input.combine(ch_contigs)
  
  // get coverage per 100bp bin per chromosome
  COV(ch_input_per_chr, fasta)

  // partition the genome into even chunks by coverage
  // TODO: fix size to use number of contigs
  ch_partition =
    COV.out.part
    .groupTuple(size: 27)
    .combine(INTERVALS.out.intervals)
  PART(ch_partition)

  // run dsa per partition
  ch_partitions =
    PART.out
    .transpose()
    .map { meta, part_bed ->
      // extract partition number from filename (e.g., "part_1.bed" -> 1)
      def partition_num = part_bed.name.replaceAll(/part_(\d+)\.bed/, '$1').toInteger()
      tuple(meta, partition_num, part_bed)
    }
  DSA(ch_input.combine(ch_partitions, by: 0), fasta, dsa_noise_bed)

  // run variantcaller per partition, merge outputs
  VAR(DSA.out.var)
  MERGE_VARS(VAR.out.merge_vars.groupTuple(size: params.jobs))

  // run indelcaller per partition, merge outputs
  INDEL(DSA.out.indel, fasta)
  // POST(INDEL.out.done, fasta, post_triNuc)

}