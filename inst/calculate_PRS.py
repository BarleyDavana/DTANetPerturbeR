import os
import sys
import numpy as np
import networkx as nx
import pandas as pd
from enm.Enm import *

def calculate_PRS(network_file, affinity_file, output_folder):
    """
    Calculating drug scores using PRS methods for gene-target-drug network and saves results to output directory.

    Args:
        network_file (str): The path to the PPI network file.
        output_dir (str): The path to the output directory where the results will be saved.
        affinity_file(str): The path to the predicted binding affinity file.
    Returns:
        None
    """
    # Read and process network data
    enm = Enm('PPIN')
    enm.read_network(network_file, sep='\t')
    enm.gnm_analysis(normalized=False)
    enm.cluster_matrix(enm.prs_mat)
    enm.df.to_csv(os.path.join(output_folder, 'pcc_df.csv'), index=True, index_label='orf_name_id')
    np.savetxt(os.path.join(output_folder, 'prs_df.txt'), enm.prs_mat)
    np.savetxt(os.path.join(output_folder, 'prs_mat_df.txt'), enm.prs_mat_df)

    # Read and process binding affinity data
    with open(affinity_file, "r") as f:
        lines = f.readlines()
    data = []
    for line in lines[3:-1]:
        row = [field.strip() for field in line.split("|")][1:-1]
        data.append(row)
    df1 = pd.DataFrame(data, columns=["Rank", "Drug Name", "Target Name", "Binding Score"])
    df1['Binding Score'] = pd.to_numeric(df1['Binding Score'], errors='coerce')
    df1 = df1.dropna(subset=['Binding Score'])

    # Merge network and drug data
    df2 = pd.read_csv(os.path.join(output_folder, 'pcc_df.csv'))
    df = pd.merge(df1, df2, left_on='Target Name', right_on='orf_name')
    df = df[['Drug Name', 'Target Name', 'Binding Score', 'sens']]
    df['ps'] = df['Binding Score'] * df['sens']
    df['ps'] = (df['ps'] - df['ps'].min()) / (df['ps'].max() - df['ps'].min())
    df = df.sort_values(by='ps', ascending=False)
    new_df = df.groupby(['Drug Name'], as_index=False).agg({'Target Name': lambda x: ', '.join(sorted(set(x))), 'ps': 'sum'})
    new_df = new_df.rename(columns={'ps': 'PS'})
    new_df['Target Name'] = new_df['Target Name'].apply(lambda x: ', '.join(sorted(set(x.split(', ')))))
    new_df = new_df.sort_values(by='PS', ascending=False)

    # Save output files
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    df.to_csv(os.path.join(output_folder, 'prs_dti_ps.csv'), index=False)
    new_df.to_csv(os.path.join(output_folder, 'drug_PS.csv'), index=False)

if __name__ == '__main__':
    network_file = sys.argv[1]
    affinity_file = sys.argv[2]
    output_folder = sys.argv[3]
    calculate_PRS(network_file, affinity_file, output_folder)
