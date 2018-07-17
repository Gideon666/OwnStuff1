#!/bin/bash

work_dir=$2

if [[ $1 != *"_"*  ]]; then
    echo "nothing to do "$1
    exit
fi

echo -e "creating Pool : "$1

declare -A pool_files
pool_ids=$(echo $1 | tr "_" "\n")

#if [[ ${#pool_ids[@]} < 1 ]]; then
#    echo ${#pool_ids[@]}
#    exit
#fi

for id in $pool_ids; do
    pool_files[${#pool_files[@]}]="${work_dir}${id}_pool.fasta"
done

new_file_name=$1

cat ${pool_files[@]} > "${work_dir}${new_file_name}_pool.fasta"
    
