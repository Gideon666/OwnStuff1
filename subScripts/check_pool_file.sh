#!/bin/bash

if [[ $1 != *"_"*  ]]; then
    echo $1
    exit
fi

declare -A pool_files
pool_ids=$(echo $1 | tr "_" "\n")

#if [[ ${#pool_ids[@]} < 1 ]]; then
#    echo ${#pool_ids[@]}
#    exit
#fi

for id in $pool_ids; do
    pool_files[${#pool_files[@]}]="${id}_pool.fasta"
done

new_file_name=$1

cat ${pool_files[@]} > "${new_file_name}_pool.fasta"
    
