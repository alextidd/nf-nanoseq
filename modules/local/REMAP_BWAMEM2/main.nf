process REMAP_BWAMEM2 {
    tag "${meta.id}"

    input:
        tuple val(meta), path(bwamem), path(cram)
        path index_dir
        val min_mapQ

    output:
        tuple val(meta), path("sort/${meta.id}.cram"), path("sort/${meta.id}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: "-n -tags BC,QT,mb,rb -b \'-T $min_mapQ -Y\'"
        def args2 = task.ext.args ?: '-n'
        """
        mkdir -p sort
        bwa_mem.pl -p bwamem $args -bwamem2 -cram -t $task.cpus -mt $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.id} $cram
        # need the mark process so the final cram file gets placed in the correct location
        bwa_mem.pl -p mark -n $args -bwamem2  -cram -t $task.cpus -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.id} $cram
        bwa_mem.pl -p stats -n $args -bwamem2  -cram -o ./bwamem -r ${index_dir}/genome.fa -s ${meta.id} $cram
        samtools sort -@ $task.cpus $args2 -O cram -m 2G -o ./sort/${meta.id}.cram ./bwamem/${meta.id}.cram
        touch sort/${meta.id}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
            bwa_mem.pl: \$(bwa_mem.pl -v )
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args ?: ''
        """
        sleep \$[ ( \$RANDOM % 20 )  + 1 ]s
        mkdir -p bwamem
        mkdir -p sort
        touch ./sort/${meta.id}.cram
        touch ./sort/${meta.id}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
            bwa_mem.pl: \$(bwa_mem.pl -v )
        END_VERSIONS
        """
}