#!/usr/bin/env bash

set -o errexit
set -o nounset

JSON_FILE="${1}"

jq -r 'keys[]' "${JSON_FILE}" | 
    while read -r species; do
        ./run_star_index.sh "${species}" ./genomeFileUrls.json "./indexed_genome/${species}"
    done
