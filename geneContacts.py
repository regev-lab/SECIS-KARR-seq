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




#%%
df = pd.read_csv("gene_canonical_transcript.csv", index_col = 0)
#%%
gene_transcript_ids = df['Canonical Transcript']

#%%
gene_lengths = df['Length']

#%%
gene_names = list(gene_transcript_ids.index)
#%%
gene_names

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

#%% Optimized version
def create_transcript_contact_files(file_path, gene_names, dirName, total_lines = None):
    
    transcript_ids = {gene_transcript_ids[gene_name]: gene_name for gene_name in gene_names}

    if total_lines == None:
        with gz.open(file_path, "rt") as f:
            total_lines = sum(1 for _ in f)
            
        print(f"Total Lines: {total_lines}")

    gene_contact_lines = {gene_name: [] for gene_name in gene_names}

    gene_contact_rows = []
    

    with gz.open(file_path, "rt") as f:
        for line in tqdm(f, total=total_lines, desc=f"Processing: {os.path.basename(file_path)}"):
            fields = line.strip().split()
            read_id, chr1, pos1, chr2, pos2, strand1, strand2 = fields[:7]
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
        
    gene_contact_df = pd.DataFrame(gene_contact_rows) 
    gene_contact_df = gene_contact_df.astype({"position1": int, "position2": int})
    
    gene_contact_df.to_csv(os.path.join(output_dir, "genes_contact_rows.csv"), index=False)
    
    gene_contact_df_negative = gene_contact_df[(gene_contact_df["strand1"]+gene_contact_df["strand2"])=="--"]
    generate_contact_maps(gene_contact_df_negative, 100, os.path.join(output_dir, "contact_maps.pdf"))
    
    return gene_contact_counts, gene_contact_df

#%% Plotting
from matplotlib.backends.backend_pdf import PdfPages

def generate_contact_maps(df, m=30, pdf_output=None):
    pdf = PdfPages(pdf_output) if pdf_output else None  # Open a PDF if output is provided

    for gene in gene_names:
        # Extract positions for the current gene
        gene_contact_positions = df[df["gene_name"] == gene][["position1", "position2"]].values
        if len(gene_contact_positions) == 0:
            continue
        
        gene_length = gene_contact_positions.max()
        
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

        # Plot the contact map
        fig, ax = plt.subplots(figsize=(6, 6))
        im = ax.imshow(contact_map, cmap="hot", origin="lower", interpolation="nearest")
        fig.colorbar(im, label="Contact Count")
        ax.set_title(f"Summed Contact Map for {gene} ({m}x{m})")
        ax.set_xlabel("Binned Position 1")
        ax.set_ylabel("Binned Position 2")

        if pdf:
            pdf.savefig(fig)  # Save to PDF
        else:
            plt.show()  # Show on screen

        plt.close(fig)  # Close figure to free memory

    if pdf:
        pdf.close()  # Close PDF file
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
contact_df1.info()
#%%
data_dir = "data/pairs/GSE166155_RAW"

file = "GSM5064767_G1_kethoxal-K562_M15_R01.dedup.pairs.gz"

file_path = os.path.join(data_dir,file)

dir_name = "_".join(file.split(".")[0].split("_")[1:])


gene_contact_counts, contact_df1 = create_transcript_contact_files(file_path, gene_names, dir_name)
#%%
data_dir = "data/pairs/GSE166155_RAW"

file = "GSM5064768_G1_kethoxal-K562_M15_R02.dedup.pairs.gz"

file_path = os.path.join(data_dir,file)
dir_name = "_".join(file.split(".")[0].split("_")[1:])


gene_contact_counts2, contact_df2 = create_transcript_contact_files(file_path, gene_names, dir_name)

#%%
dir_name = "_".join(file.split(".")[0].split("_")[1:])

#%%
gene = "TXNRD1"
#%%
contact_df2.info()

#%%
(contact_df1["strand1"]+contact_df1["strand2"]).value_counts()

#%%
contact_df1["gene_name"].value_counts()
#%%
contact_df1[(contact_df1["strand1"]+contact_df1["strand2"])=="--"]
#%%

contact_df1.info()
#%%
gene_contact_positions = contact_df1[contact_df1["gene_name"]==gene][["position1", "position2"]]
#%%
contact_df_negative = contact_df1[(contact_df1["strand1"]+contact_df1["strand2"])=="--"]
generate_contact_maps(contact_df_negative,100, "output.pdf")