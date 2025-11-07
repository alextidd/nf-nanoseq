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