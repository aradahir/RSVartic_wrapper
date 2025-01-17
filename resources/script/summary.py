import glob
import os
import pandas as pd

#combine figures and consensus
if os.path.isdir('result/artic/Bysubtypes/RSVA'):
   
    os.system('pdfunite result/artic/Bysubtypes/RSVA/figures/figure_annotation/*.pdf result/artic/Bysubtypes/RSVA/figures/coverage_depth_RSVA.pdf')
    os.system('cat result/artic/Bysubtypes/RSVA/consensus/*.fa > result/artic/Bysubtypes/RSVA/RSVA_consensus.fasta')

if os.path.isdir('result/artic/Bysubtypes/RSVB'):
    
    os.system('pdfunite result/artic/Bysubtypes/RSVB/figures/figure_annotation/*.pdf result/artic/Bysubtypes/RSVB/figures/coverage_depth_RSVB.pdf')
    os.system('cat result/artic/Bysubtypes/RSVB/consensus/*.fa > result/artic/Bysubtypes/RSVB/RSVB_consensus.fasta')

#summary F gene mutation
file_counts = {}

# Iterate through each CSV file in the folder
if os.path.isdir('result/F_proteins/RSVA/mutation'):
    for filename in os.listdir('result/F_proteins/RSVA/mutation'):
        if filename.endswith('.csv'):  # Process only CSV files
            file_path = os.path.join('result/F_proteins/RSVA/mutation', filename)
            
            df = pd.read_csv(file_path)

            y_count = df['Has_mutation'].value_counts().get('Y', 0)  # Default to 0 if 'Y' not found
            
            file_counts[filename] = y_count
    
    with open('result/F_proteins/RSVA/RSVA_mutation.txt', 'w') as file:
        file.write("Summary mutations in F protiens of RSVA \nsample name: count of mutation\n")
        for filename, count in file_counts.items():
            file.write(f"{filename}: {count}\n")

if os.path.isdir('result/F_proteins/RSVB/mutation'):
    for filename in os.listdir('result/F_proteins/RSVB/mutation'):
        if filename.endswith('.csv'):  # Process only CSV files
            file_path = os.path.join('result/F_proteins/RSVB/mutation', filename)
            
            df = pd.read_csv(file_path)

            y_count = df['Has_mutation'].value_counts().get('Y', 0)  # Default to 0 if 'Y' not found
            
            file_counts[filename] = y_count
    
    with open('result/F_proteins/RSVB/RSVB_mutation.txt', 'w') as file:
        file.write("Summary mutations in F protiens of RSVB \nsample name: count of mutation\n")
        for filename, count in file_counts.items():
            file.write(f"{filename}: {count}\n")