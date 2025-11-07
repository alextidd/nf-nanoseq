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