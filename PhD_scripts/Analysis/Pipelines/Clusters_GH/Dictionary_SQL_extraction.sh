#!/usr/bin/env bash

cd ~/Documents/Pipeline/Cluster_gene_extraction || exit

touch Cluster.tab
touch LocusProtID.tab
touch Manually_curated_cyanorak.tab
touch Dict_Cluster.pickle
touch Dict_LocusProtID.pickle

python ~/Documents/Scripts/My_Scripts/pickle_creation.py

for file in ~/Documents/Pipeline/Cluster_gene_extraction/*faa

do

echo "Working on ${file}"

python ~/Documents/Scripts/My_Scripts/Gene_cluster_extraction.py ${file}

done

sed -e "s/\/home\/sls_alberto\/Documents\/Pipeline\/Cluster_gene_extraction\///g" Cluster.tab > ClusterArrayNonDef.tab
sed -e "s/ID://g" ClusterArrayNonDef.tab > ClusterArray.tab

sed -e "s/\/home\/sls_alberto\/Documents\/Pipeline\/Cluster_gene_extraction\///g" Dict_Cluster.pickle > Dict_ClusterNonDef.pickle
sed -e "s/ID://g" Dict_ClusterNonDef.pickle > Dict_Cluster_Def.pickle


echo "Done"


for file in ~/Documents/Pipeline/Cya_clusterPipeline/*gbk

do

python ~/Documents/Scripts/My_Scripts/Gene_ProtIDLocus_extraction.py $file

done

sed -e "s/ /_/g" Dict_LocusProtID.pickle > Dict_LocusProtID_Def.pickle

python ~/Documents/Scripts/My_Scripts/Gene_Cluster_ProtID_comparison.py

echo "Done"
