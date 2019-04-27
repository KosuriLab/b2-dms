#!/bin/bash
#
# bash strictmode
set -euo pipefail
IFS=$'\n\t'

allreads='false'
while getopts 'a' flag; do
    case "${flag}" in
        a) allreads='true' ;;
        *) error "Unexpected option ${flag}" ;;
    esac
done

if [[ "${allreads}" = 'true' ]]; then
    # shift positional vars to account for -a
    join -j1 -o 2.2 1.2 \
        <(zcat "$2" | awk '{print $3":"$4":"$5":"$6, $9}' | sort -k 1b,1) \
        <(zcat "$3" | awk '{print $3":"$4":"$5":"$6, $9}' | sort -k 1b,1)
else
    join -j1 -o 2.2 1.2 \
        <(zcat "$1" | awk '$11 == 1{print $3":"$4":"$5":"$6, $9}' | sort -k 1b,1) \
        <(zcat "$2" | awk '$11 == 1{print $3":"$4":"$5":"$6, $9}' | sort -k 1b,1)
fi
