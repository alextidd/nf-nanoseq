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
  publishDir "${params.out_dir}/${meta.donor_id}/", mode: 'copy'

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
  publishDir "${params.out_dir}/${meta.donor_id}/", mode: 'copy'

  input:
  tuple val(meta),
        path(bundled_bam), path(bundled_bai),
        path(filtered_bam), path(filtered_bai)
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
    -dedup $filtered_bam \
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
  tuple path("gIntervals.dat"), path("intervals.bed")

  script:
  def int_include = params.int_include ? "--include ${params.int_include}" : ""
  def int_exclude = params.int_exclude ? "--exclude ${params.int_exclude}" : ""
  def int_larger = params.int_larger ? "--larger ${params.int_larger}" : ""
  """
  nanoseq_intervals.py \
    --ref $fasta \
    $int_include $int_exclude $int_larger
  """
  stub:
  """
  touch gIntervals.dat
  touch intervals.bed
  """
}

process COV {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/analysis/${meta.id}/cov/",
    mode: 'copy', pattern: "tmpNanoSeq/cov/*"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)
  tuple path(intervals_dat), path(intervals_bed)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: done
  tuple val(meta), path("tmpNanoSeq/cov/*"), emit: cov

  script:
  def int_include = params.int_include ? "--include ${params.int_include}" : ""
  def int_exclude = params.int_exclude ? "--exclude ${params.int_exclude}" : ""
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.30.0 
  module load bcftools-1.19/python-3.11.6

  # run
  # TODO: handle -k and -j
  runNanoSeq.py \
    --threads $task.cpus \
    --duplex $duplex_bam \
    --normal $normal_bam \
    --ref $fasta \
    cov \
    --gintervals $intervals_dat \
    -Q ${params.cov_q} \
    --larger ${params.int_larger} \
    ${int_include} \
    ${int_exclude}
  """
  stub:
  """
  mkdir -p tmpNanoSeq/cov
  touch cov/1.done
  touch cov/1.cov.bed.gz
  touch cov/args.json
  touch cov/nfiles
  """
}

process PART {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/analysis/${meta.id}/part/",
    mode: 'copy', pattern: "tmpNanoSeq/part/*"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: done
  tuple path(fasta), path(fai)
  tuple val(meta), path("tmpNanoSeq/part/*"), emit: part

  script:
  def excludeBED = params.part_excludeBED ? "--excludeBED ${params.part_excludeBED}" : ""
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.30.0 
  module load bcftools-1.19/python-3.11.6

  # run
  # TODO: handle -n 100
  runNanoSeq.py \
    --threads $task.cpus  \
    --normal $normal_bam \
    --duplex $duplex_bam \
    --ref $fasta \
    part \
    --jobs ${params.jobs} \
    --excludeCov ${params.part_excludeCov} \
    ${excludeBED}
  """
  stub:
  """
  mkdir -p tmpNanoSeq/part
  touch tmpNanoSeq/part/1.done
  touch tmpNanoSeq/part/intervalsPerCPU.dat
  touch tmpNanoSeq/part/args.json
  """
}

process DSA {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/analysis/${meta.id}/dsa/",
    mode: 'copy', pattern: "tmpNanoSeq/dsa/*"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: done
  tuple val(meta), path("tmpNanoSeq/dsa/*"), emit: dsa

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
    dsa \
    --max_index ${params.jobs} \
    --mask ${params.dsa_noise_bed} \
    -d ${params.dsa_d} \
    -q ${params.dsa_q} \
    -M ${params.dsa_M} \
    --no_test
  """
  stub:
  """
  mkdir -p tmpNanoSeq/dsa

  """
}

process VAR {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/analysis/${meta.id}/var/",
    mode: 'copy', pattern: "tmpNanoSeq/var/*"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: done
  tuple val(meta), path("tmpNanoSeq/var/*"), emit: var

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.30.0 
  module load bcftools-1.19/python-3.11.6

  # run
  # TODO: handle -k $LSB_JOBINDEX_END and -j $LSB_JOBINDEX
  runNanoSeq.py \
    --duplex $duplex_bam \
    --normal $normal_bam \
    --ref $fasta \
    var \
    --max_index ${params.jobs} \
    -a $params.var_a \
    -b $params.var_b \
    -c $params.var_c \
    -d $params.var_d \
    -f $params.var_f \
    -i $params.var_i \
    -m $params.var_m \
    -n $params.var_n \
    -p $params.var_p \
    -q $params.var_q \
    -r $params.var_r \
    -v $params.var_v \
    -x $params.var_x \
    -z $params.var_z
  """
}

process INDEL {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/analysis/${meta.id}/indel/",
    mode: 'copy', pattern: "tmpNanoSeq/indel/*"

  input:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta),
        path(duplex_bam), path(duplex_bai),
        path(normal_bam), path(normal_bai), emit: done
  tuple val(meta), path("tmpNanoSeq/indel/*"), emit: indel

  script:
  """
  # modules
  module add samtools-1.19/python-3.12.0 
  module load perl-5.30.0 
  module load bcftools-1.19/python-3.11.6

  # run
  # TODO: handle -k $LSB_JOBINDEX_END and -j $LSB_JOBINDEX
  runNanoSeq.py \
    --threads $task.cpus  \
    --normal $normal_bam \
    --duplex $duplex_bam \
    --ref $fasta \
    indel \
    --sample ${meta.id} \
    --rb ${params.indel_rb} \
    --t3 ${params.indel_t3} \
    --t5 ${params.indel_t5} \
    -z ${params.indel_z} \
    -v ${params.indel_v}
    -a ${params.indel_a} \
    -c ${params.indel_c}
  """
}

process POST {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/analysis/${meta.id}/post/",
    mode: 'copy', pattern: "tmpNanoSeq/post/*"

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

  // get chromosome intervals
  INTERVALS(fasta)

  // preprocess duplex
  PREPROCESS(ch_bams.duplex)

  // effi
  PREPROCESS.out.effi.view()
  EFFI(PREPROCESS.out.effi, fasta)

  // rejoin duplex and normal channels by donor id
  ch_input =
    ch_bams.normal
    .map { meta, bam, bai -> [ meta.donor_id, [meta.id, bam, bai] ] }
    .join(
        PREPROCESS.out.cov.map { meta, bam, bai -> [ meta.donor_id, [meta, bam, bai] ] }
    )
    .map { donor_id, normal, duplex ->
        def (normal_meta, normal_bam, normal_bai) = normal
        def (duplex_meta, duplex_bam, duplex_bai) = duplex
        tuple(duplex_meta, duplex_bam, duplex_bai, normal_bam, normal_bai)
    }
  
  // run nanoseq analysis
  COV(ch_input, fasta, INTERVALS.out)
  PART(COV.out.done, fasta)
  // DSA(PART.out.done, fasta)
  // VAR(DSA.out.done, fasta)
  // INDEL(VAR.out.done, fasta)
  // POST(INDEL.out.done, fasta, post_triNuc)

}