from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import pysam
import pandas as pd
import os

#Connect with command line 
parser = argparse.ArgumentParser(description="A script to find the mutation in F genes.")
parser.add_argument("-n", "--name", type=str,
                    help="given name for plotting")
parser.add_argument("-t", "--type", type=str,
                    help="RSV subtype. Optional: RSVA or RSVB")
parser.add_argument("-s", "--segment", type=str,
                    help="Optional: all or F")
parser.add_argument("-m", "--min", type=int,
                    help="minimum threshold for softmasking")
args = parser.parse_args()

name = args.name
subtype = args.type
segment = args.segment
coverage_threshold = args.min

samfile = pysam.AlignmentFile(f"result/artic/{name}_polished.primertrimmed.rg.sorted.bam", "rb")
ref_seq = SeqIO.read(f"resources/mutation/{subtype}_ref.fasta", "fasta")

#F genes position
if subtype == 'RSVA' and segment == 'F':
    start =  5619
    end = 7521
elif subtype == 'RSVB' and segment == 'F': #check F frame of type B
    start = 5692
    end = 7416
elif segment == 'all':
    start = 0
    end = len(ref_seq)

ref_slice = str(ref_seq.seq[start:end])

# Create an array to hold the most frequent base at each position
coverage = [0] * (end - start)
consensus = list(ref_slice)
total_coverage = list(ref_slice)

for pileupcolumn in samfile.pileup(subtype, start, end, min_base_quality=0):
    pos = pileupcolumn.pos - start  # Relative position
    if 0 <= pos < len(coverage):
        base_counts = {"A": 0, "T": 0, "G": 0, "C": 0}
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            base = pileupread.alignment.query_sequence[pileupread.query_position]
            base_counts[base.upper()] += 1

        # Calculate total coverage at this position
        total_coverage[pos] = sum(base_counts.values())
        
        # Determine the consensus base if coverage meets the threshold
        if total_coverage[pos] >= coverage_threshold:
            consensus[pos] = max(base_counts, key=base_counts.get)
        else:
            consensus[pos] = 'N'  # Assign 'N' if below the threshold

# Save consensus sequence
    
consensus_sequence = "".join(consensus)
if segment == 'F':
    if 'N' not in consensus_sequence:
        with open(f"result/F_proteins/{subtype}/consensus/{name}_F_consensus_nucleotide.fa" , "w") as f:
            renamed = name.split("_")[1]
            f.write(f">F_{renamed}\n{consensus_sequence}\n")

        translated_F_genes = Seq(consensus_sequence).translate()
        if subtype == 'RSVA':
            F_genes_without_signals = translated_F_genes[translated_F_genes.find('M'):]
        elif subtype == 'RSVB':
            F_genes_without_signals = translated_F_genes[translated_F_genes.find('M'):]

        mutation_RSV = pd.read_csv(f'resources/mutation/Mutations_{subtype}.csv')
        RSV_F_ref = SeqIO.read(f"resources/mutation/{subtype}_F_protien.fasta", "fasta")
        mutation_RSV['mutation_position'] = mutation_RSV['Mutation'].str.findall(r'(\d+)')

        #check if the subtype and protien are matching, then write the table
        for index, row in mutation_RSV.iterrows():
            if row['Protein'] == 'F' :
                mutation_array = []
                for mutation in row['mutation_position']:
                    mutation_in_sample = str(RSV_F_ref[int(mutation)-1])+ mutation + str(F_genes_without_signals[int(mutation)-1])
                    mutation_array.append(mutation_in_sample)
                mutation_array = '+'.join(mutation_array)

                mutation_RSV.at[index, 'sample_mutation'] = mutation_array
                if row['Mutation'] == mutation_array:
                    mutation_RSV.at[index, 'Has_mutation'] = 'Y'
                else: mutation_RSV.at[index, 'Has_mutation'] = 'N'
                
        # Drop the 'mutation_position' column
        df = mutation_RSV.drop(columns=['mutation_position'])

        # Write the modified DataFrame to a CSV file
        df.to_csv(f'result/F_proteins/{subtype}/mutation/{name}_mutation.csv', index=False) 
    else:
        print(f"Error: Found gap or N in F proteins in {name}, skipping")

elif segment == 'all':
    with open(f"result/artic/Bysubtypes/{subtype}/consensus/{name}.fa" , "w") as f:
        renamed = name.split("_")[1]
        f.write(f">{renamed}\n{consensus_sequence}\n")
    os.system(f'mkdir result/artic/Bysubtypes/{subtype}/artic_consensus/')
    os.system(f'cp result/artic/{name}_polished.consensus.fasta result/artic/Bysubtypes/{subtype}/artic_consensus/{name}.fasta')

