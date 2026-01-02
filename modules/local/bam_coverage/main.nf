/*
 * This process takes a BAM file as input and crate the bigwig.
 */
process BAM_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0 ' :
        'biocontainers/deeptools:3.5.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(scale_f)               
    val bin_size                      

    output:
    tuple val(meta), path("${prefix}.bw"), emit: bigwig
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bamCoverage \\
        -b ${bam} \\
        -o ${prefix}.bw \\
        --numberOfProcessors ${task.cpus} \\
        --scaleFactor ${scale_f} \\
        --binSize ${params.bin_size} \\
         ${args}
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$( bamCoverage --version | sed -e 's/bamCoverage //g' )
    END_VERSIONS
    """

}

