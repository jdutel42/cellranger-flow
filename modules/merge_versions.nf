process MERGE_VERSIONS {
    tag "MrgV_${params.bcl_id}_Merge_Versions_${params.protocol_prefix}" // Log the run ID for traceability for "Mrg" = "Merge Versions"

    publishDir (
        path    : params.run_traceability_log_dir,
        mode    : 'copy',
        pattern : "*_versions.yaml"
    )

    input:
        path ch_versions

    output:
        path "Merge_versions_${params.bcl_id}.yaml", emit: merged_versions

    script:
        """
        set -euo pipefail

        # Load shared logging helpers
        source "${params.logging_script}"

        log_start "Step 5: Merging versions information for run_id = ${params.run_id}..."

        tmp_file="merged_versions.tmp"
        : > "\$tmp_file"

        for vf in *.yml *.yaml; do
            [ -e "\$vf" ] || continue
            if [ "\$vf" = "Merge_versions_${params.bcl_id}.yaml" ]; then
                continue
            fi
            cat "\$vf" >> "\$tmp_file"
            echo "" >> "\$tmp_file"
        done

        mv "\$tmp_file" "Merge_versions_${params.bcl_id}.yaml"

        log_success "Versions information merged successfully for run_id = ${params.run_id} and saved to ${params.run_traceability_log_dir}/Merge_versions_${params.bcl_id}.yaml !"
        """

    stub:
    """
    cat > "Merge_versions_${params.bcl_id}.yaml" <<EOF
"MERGE_VERSIONS":
    "run_id": "${params.run_id}"
    "timestamp": "${params.today_date}"
    "all_versions_merged": true
EOF
    """
}
