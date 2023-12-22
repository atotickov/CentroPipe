#!/bin/bash

print_usage() {
        echo "  -i      Assemblies file in fasta format ("assembly1.fasta assembly2.fasta assembly3.fasta")."
        echo "  -n      1n for each assembly ("22 21 16")" 
}

assemblies=();
ChrN=();

while getopts 'i:n:' flag; do
        case "${flag}" in
                i) IFS=' ' read -ra assemblies <<< "${OPTARG}" ;;
                n) IFS=' ' read -ra ChrN <<< "${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

workflow_path=$(pwd)
echo "Workflow Path: ${workflow_path}"

for ((i=0; i<${#assemblies[@]}; i++)); do

        prefix=${assemblies[i]%.*}
        echo "Assembly: ${assemblies[i]} | Chr number (1n): ${ChrN[i]} | Prefix: ${prefix}"
        
        mkdir -p ${workflow_path}/${prefix}
        mv ${assemblies[i]} ${workflow_path}/${prefix}/
        cd ${workflow_path}/${prefix}/
        
        samtools faidx ${assemblies[i]}
        
        cat ${assemblies[i]}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 > ${prefix}.len 
        
        cat ${assemblies[i]}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 | head -n ${ChrN[i]} > ${prefix}.lengths 
        
        cat ${prefix}.lengths | awk '{print $1}' > ${prefix}.whitelist 
        
        cat ${prefix}.whitelist | awk '{print $1"\t"$1}' > ${prefix}.syn 
        
        cat ${prefix}.syn | awk '{print $2}' > ${prefix}.orderlist
        
        cd ${workflow_path}/
done
