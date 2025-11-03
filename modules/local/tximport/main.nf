/*
 * This process takes the kallisto quant output as input and generates the gene expressio nmatrix.
 */
process TXIMPORT {
    tag "tximport_all_samples"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://pdichiaro/lmdseq:latest'

    input:
    path all_folders   // this will be a list of kallisto result folders
    path gtf                            
    path reference                      // Optional: can be empty list []                       

    output:
    path("EX_reads_RAW.txt"), emit: matrix   // Gene expression matrix
    path("versions.yml"), emit: versions

    script:
    def reference_arg = (reference.toString() != "[]" && reference.toString() != "") ? "-r ${reference}" : ""
    def folders_list = all_folders.join(' ')
    """
    #!/bin/bash
    set -e  # Exit on any error
    
    echo "=== TXIMPORT MODULE ==="
    echo "Shell: \$0"
    echo "Available tools test:"
    which R || echo "WARNING: R not found"
    which Rscript || echo "WARNING: Rscript not found"
    which wc || echo "WARNING: wc not found"
    which head || echo "WARNING: head not found"
    which tr || echo "WARNING: tr not found"
    which cut || echo "WARNING: cut not found"
    echo "GTF file: ${gtf}"
    echo "Reference: ${reference}"
    echo "Number of kallisto folders: ${all_folders.size()}"
    echo "Kallisto directories: ${folders_list}"
    
    # Validate inputs exist
    [ -f "${gtf}" ] || { echo "ERROR: GTF file not found: ${gtf}"; exit 1; }
    
    # Check kallisto outputs (more robust approach)
    echo "Validating kallisto output directories..."
    ls -la
    for dir_path in *; do
        if [ -d "\$dir_path" ]; then
            echo "Checking directory: \$dir_path"
            if [ -f "\$dir_path/abundance.h5" ]; then
                echo "✓ Valid: \$dir_path"
            else
                echo "WARNING: abundance.h5 not found in \$dir_path"
            fi
        fi
    done

    # Run tximport
    echo "Running tximport R script..."
    Rscript ${projectDir}/bin/create_reference_db.R \\
        -g ${gtf} \\
        ${reference_arg} \\
        -i "${folders_list}" \\
        -o EX_reads_RAW.txt || { echo "ERROR: R script failed"; exit 1; }

    # Validate output
    [ -f "EX_reads_RAW.txt" ] || { echo "ERROR: Output file not created"; exit 1; }
    
    # Count genes and samples safely
    GENE_COUNT=\$(wc -l < EX_reads_RAW.txt)
    SAMPLE_COUNT=\$(head -n1 EX_reads_RAW.txt | tr '\\t' '\\n' | wc -l)
    SAMPLE_COUNT=\$((SAMPLE_COUNT - 1))  # Subtract 1 for gene_id column
    echo "✓ Success: \${GENE_COUNT} genes, \${SAMPLE_COUNT} samples"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | head -n1 | cut -d' ' -f3 2>/dev/null || echo "unknown")
        tximport: \$(Rscript -e "cat(packageVersion('tximport'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
