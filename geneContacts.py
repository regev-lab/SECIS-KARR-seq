#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:15:48 2025

@author: frankliu
"""
import os
import pandas as pd
import gzip as gz
from tqdm import tqdm

import numpy as np

import matplotlib.pyplot as plt

from scipy.ndimage import zoom

import matplotlib


#%%
df = pd.read_csv("gene_canonical_transcript.csv", index_col = 0)
#%%
gene_transcript_ids = df['Canonical Transcript']

#%%
gene_lengths = df['Length']

#%%
gene_names = list(gene_transcript_ids.index)

#%%
secis_df = pd.read_csv("secis_info.csv", index_col = 0)

#%%
secis_df["transcript_length"] =secis_df.index.map(gene_lengths)
secis_df["dist_stop"] = secis_df["dist_stop"].str.split(" ").str[0].astype("int64")
secis_df["secis_end"] = secis_df["transcript_length"] - secis_df["dist_end"]
secis_df["secis_start"] = secis_df["secis_end"] - secis_df["secis_length"]
secis_df["stop_codon"] = secis_df["secis_end"] - secis_df["dist_stop"]

secis_df = secis_df[["secis_start", "secis_end", "stop_codon", "transcript_length"]]
#%%
secis_df.head()

#%%
def create_transcript_contact_file(file_path, gene_name, total_lines = None):
    transcriptId = gene_transcript_ids[gene_name]

    # Create directory named geneData if it does not exist
    os.makedirs("geneData", exist_ok=True)

    output_file_path = os.path.join("geneData", f"{gene_name}.pairs.gz")
    
    contact_ct = 0
    with gz.open(file_path, "rt") as f, gz.open(output_file_path, "wt") as output:
        for line in tqdm(f, total=total_lines, desc=f"Searching for {gene_name}"):
            if transcriptId in line:
                output.write(line)
                contact_ct += 1
    return contact_ct


#%% Plotting
from matplotlib.backends.backend_pdf import PdfPages

def generate_contact_maps(df, m=100, pdf_output=None, norm = "linear"):
    pdf = PdfPages(pdf_output) if pdf_output else None  # Open a PDF if output is provided

    for gene in gene_names:
        # Extract positions for the current gene
        gene_contact_positions = df[df["gene_name"] == gene][["position1", "position2"]].values
        if len(gene_contact_positions) == 0:
            continue
        
        #gene_length = gene_contact_positions.max()
        gene_length = gene_lengths[gene]
        # Initialize contact map
        contact_map = np.zeros((m, m), dtype=int)
        
        # Populate contact map
        for pos1, pos2 in gene_contact_positions:
            if 1 <= pos1 <= gene_length and 1 <= pos2 <= gene_length:
                a = int((pos1-1) * m / gene_length)
                b = int((pos2-1) * m / gene_length)
                contact_map[a, b] += 1
                contact_map[b, a] += 1
            else:
                print(f"out of bounds: {gene}, ({pos1}, {pos2}), Length: {gene_length}")
        
        #tick_positions = [int(i/m*gene_length) for i in range (m+1)]
        tick_positions = [int(i/m*gene_length) for i in range(0, m+1, 5)]
        tick_labels = [str(pos) for pos in tick_positions]
        tick_locations = [i * (m / gene_length) - 0.5 for i in tick_positions]
        # Plot the contact map
        fig, ax = plt.subplots(figsize=(6, 6))
        im = ax.imshow(contact_map, cmap="hot", origin="lower", interpolation="nearest", norm = norm)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(f"Summed {norm.capitalize()} Contact Map for {gene} (Count: {len(gene_contact_positions)})")
        ax.set_xlabel("Binned Position 1 (nt)")
        ax.set_ylabel("Binned Position 2 (nt)")
        ax.set_xticks(tick_locations)
        ax.set_xticklabels(tick_labels, rotation=45, ha="right", rotation_mode="anchor")  
        ax.set_yticks(tick_locations)
        ax.set_yticklabels(tick_labels)
            
        if gene in secis_df.index:
            for line_label in ["secis_start", "secis_end", "stop_codon"]:
                line_pos = secis_df[line_label][gene]/gene_length*m
                # Add vertical line (x-axis)
                ax.axvline(x=line_pos, color="blue", linestyle="--", linewidth=1)
                ax.text(line_pos, 0, line_label, ha="center", va="bottom", rotation=90, fontsize=10, color="blue")
        
                # Add horizontal line (y-axis)
                ax.axhline(y=line_pos, color="blue", linestyle="--", linewidth=1)
                ax.text(0, line_pos, line_label, ha="left", va="center", fontsize=10, color="blue")
            
        
        if pdf:
            pdf.savefig(fig)  # Save to PDF
        else:
            plt.show()  # Show on screen

        plt.close(fig)  # Close figure to free memory

    if pdf:
        pdf.close()  # Close PDF file
        
#%% Optimized version
def create_transcript_contact_files(file_path, gene_names, dirName, total_lines = None):
    
    transcript_ids = {gene_transcript_ids[gene_name]: gene_name for gene_name in gene_names}

    if total_lines == None:
        with gz.open(file_path, "rt") as f:
            total_lines = sum(1 for _ in f)
            
        print(f"Total Lines: {total_lines}")
    else if total_lines == -1:
        total_lines = None

    gene_contact_lines = {gene_name: [] for gene_name in gene_names}

    gene_contact_rows = []
    

    with gz.open(file_path, "rt") as f:
        for line in tqdm(f, total=total_lines, desc=f"Processing: {os.path.basename(file_path)}"):
            fields = line.strip().split()
            read_id, chr1, pos1, chr2, pos2, strand1, strand2 = fields[:7]
            
            # Check if transcriptID matches one of the genes in our dictionary
            if chr1 in transcript_ids and chr1 == chr2:
                gene_contact_lines[transcript_ids[chr1]].append(line)
                gene_contact_rows.append({
                    "read_id": read_id,
                    "gene_name": transcript_ids[chr1],
                    "transcript_id": chr1,
                    "position1": pos1,
                    "position2": pos2,
                    "strand1": strand1,
                    "strand2": strand2
                    })
    # Create directory named geneData if it does not exist
    output_dir = os.path.join("geneData", dirName)
    os.makedirs(output_dir, exist_ok=True)
    
    for gene_name in gene_names:
        output_file_path = os.path.join(output_dir, f"{gene_name}.pairs")
        with open(output_file_path, "wt") as output:
            # write all of the gene_contact_lines[gene_name] into output
            output.writelines(gene_contact_lines[gene_name])
    
    gene_contact_counts = {gene: len(contact_lines)  for gene, contact_lines in gene_contact_lines.items()}
    
    df = pd.DataFrame([
        {"Gene": gene, 
         "Transcript_Id": gene_transcript_ids[gene],
         "Length": gene_lengths[gene],
         "Contact_Counts": count}
        for gene, count in gene_contact_counts.items()
        ])
    
    df.to_csv(os.path.join(output_dir, "genes_summary.csv"), index=False)
    columns = [
    "read_id", "gene_name", "transcript_id",
    "position1", "position2", "strand1", "strand2"
    ]    
    gene_contact_df = pd.DataFrame(gene_contact_rows, columns=columns)
    gene_contact_df = gene_contact_df.astype({"position1": int, "position2": int})
    
    gene_contact_df.to_csv(os.path.join(output_dir, "genes_contact_rows.csv"), index=False)
    
    gene_contact_df_negative = gene_contact_df[(gene_contact_df["strand1"]+gene_contact_df["strand2"])=="--"]
    generate_contact_maps(gene_contact_df_negative, 100, os.path.join(output_dir, "contact_maps.pdf"))
    
    norm = "log"
    generate_contact_maps(gene_contact_df_negative, 100, os.path.join(output_dir, f"{norm}_contact_maps.pdf"), norm = norm)

    
    return gene_contact_counts, gene_contact_df

#%%
#gunzip -c GSM5064767_G1_kethoxal-K562_M15_R01.dedup.pairs.gz | grep NM_206926 | wc
#find contacts
'''
data_dir = "data/pairs/GSE166155_RAW"

file = "GSM5064767_G1_kethoxal-K562_M15_R01.dedup.pairs.gz"

file_path = os.path.join(data_dir,file)

with gz.open(file_path, "rt") as f:
    total_lines = sum(1 for _ in f)

print(f"Total Lines: {total_lines}")
gene_contact_counts = {}
for gene in gene_names:
    gene_contact_counts[gene] = create_transcript_contact_file(file_path, gene, total_lines)
'''
#%%
data_dir = "data/pairs/GSE166155_RAW"

file = "GSM5064767_G1_kethoxal-K562_M15_R01.dedup.pairs.gz"

file_path = os.path.join(data_dir,file)

dir_name = file.split(".")[0]


gene_contact_counts, contact_df1 = create_transcript_contact_files(file_path, gene_names, dir_name)
#%%
data_dir = "data/pairs/GSE166155_RAW"

file = "GSM5064768_G1_kethoxal-K562_M15_R02.dedup.pairs.gz"

file_path = os.path.join(data_dir,file)
dir_name = file.split(".")[0]


gene_contact_counts2, contact_df2 = create_transcript_contact_files(file_path, gene_names, dir_name)


#%% Generate Gene Files and Contact Maps for all Files
data_dir = "data/pairs/GSE166155_RAW"
for file in os.listdir(data_dir):
    file_path = os.path.join(data_dir,file)
    dir_name = file.split(".")[0]
    gene_contact_counts, contact_df = create_transcript_contact_files(file_path, gene_names, dir_name)
#%% Created Combined Contact DF and find contact map
data_dir = "data/pairs/GSE166155_RAW"
columns = [
    "read_id", "gene_name", "transcript_id",
    "position1", "position2", "strand1", "strand2"
]   

gene_contact_row_dfs = []

def get_file_contact_rows(file):
    file_path = os.path.join(data_dir, file)
    dir_name = file.split(".")[0]
    contact_data_path = os.path.join("geneData", dir_name, "genes_contact_rows.csv")
    curr_df = pd.read_csv(contact_data_path)
    curr_df["file"] = dir_name # Adding the file column
    return curr_df
    
contacts_df = pd.concat( [get_file_contact_rows(file) for file in os.listdir(data_dir)])
print(contacts_df.head())
print(contacts_df.info())

#%% Reorder columns and add derived columns
contacts_df["position1"] = contacts_df["position1"].astype(int)
contacts_df["position2"] = contacts_df["position2"].astype(int)

contacts_df["strand1"] = contacts_df["strand1"].str == '+'
contacts_df["strand2"] = contacts_df["strand2"].str == '+'

contacts_df["distance"] = abs(contacts_df["position1"] - contacts_df["position2"])

# Reorder columns to place "file" after "read_id"
column_order = ["read_id", "file", "gene_name", "transcript_id", 
                "position1", "position2", "distance", "strand1", "strand2"]
contacts_df = contacts_df[column_order]

#%%
matplotlib.scale.get_scale_names()
#%% Generate Overall Contact Map
norm = "log"
generate_contact_maps(contacts_df, 100, norm = norm)

#%% Generate Export PDF
norm = "log"
generate_contact_maps(contacts_df, 100, pdf_output = f"combined_{norm}_contact_maps.pdf", norm = norm)

#%%
contacts_df = contacts_df.sort_values(by=["gene_name", "distance"], ascending=[True, False])

contacts_df.to_csv("combined_gene_contacts.csv", index = False)

#%%
contacts_df["gene_name"].value_counts().to_csv("combined_contact_counts.csv")

#%%
def get_file_contact_counts(file):
    file_path = os.path.join(data_dir, file)
    dir_name = file.split(".")[0]
    contact_data_path = os.path.join("geneData", dir_name, "genes_contact_rows.csv")
    curr_df = pd.read_csv(contact_data_path)
    return dir_name, curr_df["gene_name"].value_counts()

# Get all value_counts and store them in a dictionary
data_dict = {}
data_dir = "data/pairs/GSE166155_RAW"
for file in os.listdir(data_dir):
    dir_name, value_counts = get_file_contact_counts(file)
    data_dict[dir_name] = value_counts

# Create a DataFrame from the dictionary
summary_df = pd.DataFrame.from_dict(data_dict, orient="index").T
# Fill missing values with 0
summary_df = summary_df.fillna(0).astype(int)

# Ensure all gene names are included
summary_df = summary_df.reindex(gene_names, fill_value=0)
#%%
summary_df.head()

#%%
summary_df.to_csv("gene_contacts_file_summary.csv")



