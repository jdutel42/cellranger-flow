process MERGE_VERSIONS {
    tag "${task.process.toLowerCase()}_${run_id}"

    publishDir (
        path    : params.log_dir,
        mode    : 'copy',
        pattern : "*_versions.yaml"
    )

    input:
        val run_id
        path ch_versions
        val log_dir
        val today_date

    output:
        path "${today_date}_versions.yaml", emit: merged_versions

    script:
        """
        set -euo pipefail

        tmp_file="merged_versions.tmp"
        : > "\$tmp_file"

        for vf in *.yml *.yaml; do
            [ -e "\$vf" ] || continue
            if [ "\$vf" = "${today_date}_versions.yaml" ]; then
                continue
            fi
            cat "\$vf" >> "\$tmp_file"
            echo "" >> "\$tmp_file"
        done

        mv "\$tmp_file" "${today_date}_versions.yaml"
        """

    stub:
    """
    cat > "${today_date}_versions.yaml" <<EOF
"MERGE_VERSIONS":
    "run_id": "${run_id}"
    "timestamp": "${today_date}"
    "all_versions_merged": true
EOF
    """
}
