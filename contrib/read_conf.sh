#!/usr/bin/env bash
# Usage:   source contrib/read_conf.sh pipeline.conf section key
# Example: source contrib/read_conf.sh pipeline.conf dirs proteome_dir
conf_get () {
    local file=$1 section=$2 key=$3
    awk -F' *= *' -v sec="[$section]" -v k="$key" '
        $0==sec {found=1; next}
        /^\[/{found=0}
        found && $1==k {print $2; exit}
    ' "$file"
}

