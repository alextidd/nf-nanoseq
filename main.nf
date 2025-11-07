#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// params set in nextflow.config

// import modules
include { SAMTOOLS_INDEX  } from './modules/local/SAMTOOLS_INDEX'
include { CHECK_CONTIGS   } from './modules/local/CHECK_CONTIGS'
include { INTERVALS       } from './modules/local/INTERVALS'
include { SUBSET          } from './modules/local/SUBSET'
include { PREPROCESS      } from './modules/local/PREPROCESS'
include { EFFI            } from './modules/local/EFFI'
include { COV             } from './modules/local/COV'
include { PART            } from './modules/local/PART'
include { DSA             } from './modules/local/DSA'
include { VAR             } from './modules/local/VAR'
include { VAR_MERGE_CSV   } from './modules/local/VAR_MERGE_CSV'
include { VAR_VCF         } from './modules/local/VAR_VCF'
include { VAR_MERGE_COV   } from './modules/local/VAR_MERGE_COV'
include { INDEL           } from './modules/local/INDEL'
include { INDEL_MERGE     } from './modules/local/INDEL_MERGE'
include { PLOT            } from './modules/local/PLOT'
include { SUMMARY         } from './modules/local/SUMMARY'

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