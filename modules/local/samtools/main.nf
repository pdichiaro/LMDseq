process SAMTOOLS {
    tag "$meta.id"
    label 'process_medium'

    container 'docker://pdichiaro/lmdseq:latest' 

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), stdout, emit: scale_f
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def scale_to = task.ext.scale_to ?: 10000000

    """
    factor=\$(samtools view -c ${bam})
    factor_pe=\$(echo "\$factor / 2" | bc)
    scale_f=\$(echo "scale=7 ; ${scale_to} / \$factor_pe" | bc)

    echo \$scale_f

    cat <<EOF > versions.yml
    "${task.process}":
    samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    EOF
    """
}

