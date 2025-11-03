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
    set -e  # Exit on any error
    
    echo "=== TXIMPORT MODULE ==="
    echo "GTF file: ${gtf}"
    echo "Reference: ${reference}"
    echo "Number of kallisto folders: ${all_folders.size()}"
    echo "Kallisto directories: ${folders_list}"
    
    # Validate inputs exist
    [ -f "${gtf}" ] || { echo "ERROR: GTF file not found: ${gtf}"; exit 1; }
    
    # Create array of directories for validation
    IFS=' ' read -ra DIRS <<< "${folders_list}"
    
    # Check kallisto outputs
    for dir in "\${DIRS[@]}"; do
        [ -d "\$dir" ] || { echo "ERROR: Directory not found: \$dir"; exit 1; }
        [ -f "\$dir/abundance.h5" ] || { echo "ERROR: abundance.h5 not found in \$dir"; exit 1; }
        echo "✓ Valid: \$dir"
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
    
    echo "✓ Success: \$(wc -l < EX_reads_RAW.txt) genes, \$(awk -F'\t' 'NR==1{print NF-1}' EX_reads_RAW.txt) samples"

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        R: \$(R --version 2>&1 | head -n1 | sed 's/R version //' | sed 's/ .*//')
        tximport: \$(Rscript -e "cat(packageVersion('tximport'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
