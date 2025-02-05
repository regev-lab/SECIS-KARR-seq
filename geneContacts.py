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

#%%
df = pd.read_csv("gene_canonical_transcript.csv", index_col = 0)
#%%
gene_transcript_ids = df['Canonical Transcript']

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
    
    with gz.open(file_path, "rt") as f:
        for line in tqdm(f, total=total_lines, desc=f"Processing: {os.path.basename(file_path)}"):
            fields = line.strip().split()
            read_id, chr1, pos1, chr2, pos2, strand1, strand2 = fields[:7]
            if chr1 in transcript_ids:
                gene_contact_lines[transcript_ids[chr1]].append(line)
            if chr2 != chr1 and chr2 in transcript_ids:
                gene_contact_lines[transcript_ids[chr2]].append(line)
                
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
         "TranscriptID": gene_transcript_ids[gene], 
         "Counts": count}
        for gene, count in gene_contact_counts.items()
        ])
    
    df.to_csv(os.path.join(output_dir, "genes_summary.csv"), index=False)
    return gene_contact_counts
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

dir_name = "_".join(file.split(".")[0].split("_")[1:])


gene_contact_counts = create_transcript_contact_files(file_path, gene_names, dir_name)
#%%
data_dir = "data/pairs/GSE166155_RAW"

file = "GSM5064768_G1_kethoxal-K562_M15_R02.dedup.pairs.gz"

file_path = os.path.join(data_dir,file)
dir_name = "_".join(file.split(".")[0].split("_")[1:])


gene_contact_counts2 = create_transcript_contact_files(file_path, gene_names, dir_name)

#%%
gene_contact_counts