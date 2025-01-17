import pysam
import matplotlib.pyplot as plt
import pandas as pd
import argparse

#Connect with command line 
parser = argparse.ArgumentParser(description="A script to plot the depth results.")
parser.add_argument("-n", "--name", type=str,
                    help="given name for plotting")
parser.add_argument("-t", "--type", type=str,
                    help="RSV subtype. Optional: RSVA or RSVB")
parser.add_argument("-m", "--min", type=int,
                    help="minimum threshold for softmasking")
args = parser.parse_args()

name = args.name
subtype = args.type
coverage_threshold = args.min
# Load BAM file
bam_file = f"result/artic/{name}_polished.primertrimmed.rg.sorted.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")

# Load gene annotations from a GTF file
annotations = pd.read_csv(f"resources/mutation/{subtype}.gtf", sep='\t', comment='#', header=None)
annotations.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Filter for CDS features only
genes = annotations[annotations['feature'] == 'CDS']
region_length = samfile.get_reference_length(subtype)

# Compute depth coverage
depths = [0] * region_length
for pileup_column in samfile.pileup(subtype):
    depths[pileup_column.reference_pos] = pileup_column.nsegments

depth_data = pd.DataFrame({'Position': range(region_length), 'Depth': depths})
depth_data.to_csv(f"result/artic/Bysubtypes/{subtype}/table/{name}_depth_table.csv", index=False)

# Plot coverage
plt.figure(figsize=(15, 6))
plt.plot(range(region_length), depths, color='blue', label='Coverage')

# Add gene annotations
for _, gene in genes[genes['chrom'] == subtype].iterrows():
    plt.axvspan(gene['start'], gene['end'], color='orange', alpha=0.3)
    plt.axhline(y = coverage_threshold, color = 'r')
    mid = (gene['start'] + gene['end']) // 2
    plt.text(mid, max(depths) , gene['attribute'].split(';')[2].split('"')[1], ha='center', fontsize=8)

plt.title(f"Coverage Plot with Gene Annotations: {name}")
plt.xlabel("Genomic Position")
plt.ylabel("Depth")
plt.legend()
plt.show()
plt.savefig(f"result/artic/Bysubtypes/{subtype}/figures/figure_annotation/coverage_annotation_{name}.pdf")