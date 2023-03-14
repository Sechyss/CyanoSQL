for i in $(seq 1 110); do bbsplit.sh ref=${i}_*.fna build=${i} path=../rna_reads/; done
