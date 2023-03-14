import pandas as pd
from tqdm import tqdm

from CyanoScripts.MyPackage import translate_cluster_locustag_id


def amt_merging_genes_clusters(table: pd.DataFrame, translator):
    sql_table_full = pd.read_csv(table, index_col=0, low_memory=True)
    temporal_dict = {}
    headers = sql_table_full.columns
    for index, row in tqdm(sql_table_full.iterrows()):
        key = translate_cluster_locustag_id(index, translator)
        if key in temporal_dict.keys():
            addition = {
                key: [temporal_dict[key][0] + row['AMT22_18_18_CT.TPM'],
                      temporal_dict[key][1] + row['AMT22_3_4_CT.TPM'],
                      temporal_dict[key][2] + row['AMT22_53_53_CT.TPM'],
                      temporal_dict[key][3] + row['AMT22_75_75_CT.TPM'],
                      temporal_dict[key][4] + row['AMT22_75_75_SS.TPM'], temporal_dict[key][5] + row['AMT22_CTD4.TPM'],
                      temporal_dict[key][6] + row['AMT22_CTD75.TPM'], temporal_dict[key][7] + row['AMT23_11_15_CT.TPM'],
                      temporal_dict[key][8] + row['AMT23_11_15_DS.TPM'],
                      temporal_dict[key][9] + row['AMT23_37_46_CT.TPM'],
                      temporal_dict[key][10] + row['AMT23_37_46_DS.TPM'],
                      temporal_dict[key][11] + row['AMT23_39_49_DS.TPM'],
                      temporal_dict[key][12] + row['AMT23_3_5_CT.TPM'],
                      temporal_dict[key][13] + row['AMT23_3_5_DS.TPM'],
                      temporal_dict[key][14] + row['AMT23_54_65_CT.TPM'],
                      temporal_dict[key][15] + row['AMT23_54_65_DS.TPM']]}
            temporal_dict.update(addition)
        else:
            addition = {key: [row['AMT22_18_18_CT.TPM'], row['AMT22_3_4_CT.TPM'], row['AMT22_53_53_CT.TPM'],
                              row['AMT22_75_75_CT.TPM'], row['AMT22_75_75_SS.TPM'], row['AMT22_CTD4.TPM'],
                              row['AMT22_CTD75.TPM'], row['AMT23_11_15_CT.TPM'], row['AMT23_11_15_DS.TPM'],
                              row['AMT23_37_46_CT.TPM'], row['AMT23_37_46_DS.TPM'], row['AMT23_39_49_DS.TPM'],
                              row['AMT23_3_5_CT.TPM'], row['AMT23_3_5_DS.TPM'], row['AMT23_54_65_CT.TPM'],
                              row['AMT23_54_65_DS.TPM']]}
            temporal_dict.update(addition)
    table_clust = pd.DataFrame.from_dict(temporal_dict, orient='index', columns=headers)

    return table_clust
