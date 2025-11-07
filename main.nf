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
  intervals.py \
    --ref $fasta \
    $int_include $int_exclude $int_larger
  cut -f1 intervals.bed > contigs.txt
  """
  stub:
  """
  touch gIntervals.dat
  touch intervals.bed
  touch contigs.txt
  """
}

process SUBSET {
  tag "${meta.id}"

  input:
  tuple val(meta), path(bam), path(bai)
  path(contigs)

  output:
  tuple val(meta),
        path("${meta.id}_subset.bam"), path("${meta.id}_subset.bam.bai")

  script:
  """
  # modules
  module load samtools-1.19.2/python-3.11.6

  # filter bam to included contigs only
  samtools view --bam $bam \$(cat $contigs) --output ${meta.id}_subset.bam
  samtools index ${meta.id}_subset.bam
  """
  stub:
  """
  touch ${meta.id}_subset.bam
  touch ${meta.id}_subset.bam.bai
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
  efficiency.pl \\
    -threads $task.cpus \\
    -duplex $bundled_bam \\
    -dedup $dedup_bam \\
    -ref $fasta \\
    -out efficiency/${meta.id}
  """
  stub:
  """
  mkdir -p efficiency
  touch efficiency/${meta.id}.RBs
  touch efficiency/${meta.id}.RBs.pdf
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
  # TODO: move excludeCov and excludeBED functionalities of part.py to command line
  # excludeCov is RENanoSeq-specific
  part.py \\
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

process VAR_MERGE_CSV {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), val(partitions), path(var_tsvs)

  output:
  tuple val(meta), path("variants.csv"), emit: var_vcf
  tuple val(meta), path("*.csv"), emit: plot
  tuple val(meta), path("*.csv"), emit: summary


  script:
  """
  var_merge_csv.py \\
    --partitions ${partitions.join(',')}
  """
}

process VAR_VCF {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(vars_csv)
  tuple path(fasta), path(fai)

  output:
  tuple val(meta), path("results.muts.vcf.gz")

  script:
  """
  var_vcf.py \\
    --variants_csv $vars_csv \\
    --ref_fai $fai \\
    --sample ${meta.id} \\
    --output results.muts.vcf.gz
  """
}

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

process SUMMARY {
  tag "${meta.id}"
  publishDir "${params.out_dir}/${meta.donor_id}/${meta.id}", mode: 'copy', overwrite: true

  input:
  tuple val(meta), path(var_csvs)

  output:
  path("summary.txt")

  script:
  """
  summary.R "./"
  """
  stub:
  """
  touch summary.txt
  """
}

workflow {

  // get reference files
  fasta = [file(params.fasta, checkIfExists: true),
           file(params.fasta + ".fai", checkIfExists: true)]
  post_triNuc = file(params.post_triNuc, checkIfExists: true)
  dsa_noise_bed = [file(params.dsa_noise_bed, checkIfExists: true),
                   file(params.dsa_noise_bed + ".tbi", checkIfExists: true)]
  
  // get bams
  ch_input =
    channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
            def meta = [donor_id: row.donor_id, id: row.id]
            [meta,
             ["duplex", "normal"],
             [file(row.duplex_bam, checkIfExists: true),
              file(row.normal_bam, checkIfExists: true)]]
    }
    .transpose()
    .map { meta, types, bams ->
            def meta2 = meta + [type: types]
            [meta2, bams]
    }

  // index bams
  SAMTOOLS_INDEX(ch_input)

  // get contig intervals
  INTERVALS(fasta)

  // create channel of contigs
  ch_contigs =
    INTERVALS.out.contigs
    .splitCsv()

  // check that bam and ref contigs match
  ch_bams =
    SAMTOOLS_INDEX.out
    .branch { meta, bam, bai ->
        duplex: meta.type == "duplex"
        normal: meta.type == "normal"
    }
  ch_check =
    ch_bams.duplex
    .map { meta, bam, bai -> [ meta.subMap('donor_id', 'id'), bam, bai ] }
    .join(
        ch_bams.normal
          .map { meta, bam, bai -> [ meta.subMap('donor_id', 'id'), bam, bai ] }
    )
  CHECK_CONTIGS(ch_check, fasta)

  // subset bams to included contigs only
  ch_subset = SAMTOOLS_INDEX.out
  SUBSET(ch_subset, INTERVALS.out.contigs)

  // preprocess duplex
  ch_preprocess_duplex =
    SUBSET.out
    .filter { meta, bam, bai -> meta.type == "duplex" }
  PREPROCESS(ch_preprocess_duplex)

  // effi
  EFFI(PREPROCESS.out.effi, fasta)

  // rejoin preprocessed duplex and normal bams
  ch_nanoseq_duplex =
    PREPROCESS.out.cov
    .map { meta, bam, bai -> [meta.subMap('donor_id', 'id'), bam, bai] }
  ch_nanoseq_normal =
    SUBSET.out
    .filter { meta, bam, bai -> meta.type == "normal" }
    .map { meta, bam, bai -> [meta.subMap('donor_id', 'id'), bam, bai] }
  ch_nanoseq = ch_nanoseq_duplex.join(ch_nanoseq_normal)

  // get coverage per 100bp bin per chromosome
  ch_nanoseq_per_chr =
    ch_nanoseq
    .combine(ch_contigs.collect().map { [it] })
    .map { meta, duplex_bam, duplex_bai, normal_bam, normal_bai, chrs ->
            [groupKey(meta, chrs.size()),
             duplex_bam, duplex_bai, normal_bam, normal_bai, chrs]
    }
    .transpose()
  COV(ch_nanoseq_per_chr, fasta)

  // partition the genome into even chunks by coverage
  ch_partition =
    COV.out.part
    .groupTuple()
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
  DSA(ch_nanoseq.combine(ch_partitions, by: 0), fasta, dsa_noise_bed)

  // run variantcaller per partition, merge outputs
  VAR(DSA.out.var)
  VAR_MERGE_CSV(VAR.out.var_merge_csv.groupTuple(size: params.jobs))
  VAR_VCF(VAR_MERGE_CSV.out.var_vcf, fasta)
  VAR_MERGE_COV(VAR.out.var_merge_cov.groupTuple(size: params.jobs))

  // run indelcaller per partition, merge outputs
  INDEL(DSA.out.indel, fasta)
  INDEL_MERGE(INDEL.out.indel_merge.groupTuple(size: params.jobs))

  // plot results
  PLOT(VAR_MERGE_CSV.out.plot, post_triNuc)

  // summary statistics
  SUMMARY(VAR_MERGE_CSV.out.summary)

}