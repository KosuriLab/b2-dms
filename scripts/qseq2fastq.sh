#!/bin/bash
#
# bash strictmode
set -euo pipefail
IFS=$'\n\t'

paste -d '_' \
    <(zcat $2 | awk '{print "@"$9}') \
    <(zcat $1 | awk '{print $3":"$4":"$5":"$6, $9, $10}') | \
    awk '{print $1"\n"$2"\n+\n"$3}'
