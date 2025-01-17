import pysam
import matplotlib.pyplot as plt
import argparse

#Connect with command line 
parser = argparse.ArgumentParser(description="A script to plot the comparison between before and after primer trimming.")
parser.add_argument("-n", "--name", type=str,
                    help="given name for plotting")
parser.add_argument("-t", "--type", type=str,
                    help="RSV subtype. Optional: RSVA or RSVB")
args = parser.parse_args()

name = args.name
subtype = args.type
# Load BAM file
bam_file= f"result/artic/{name}_polished.sorted.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")

bam_file_primer = f"result/artic/{name}_polished.primertrimmed.rg.sorted.bam"
samfile_primer = pysam.AlignmentFile(bam_file_primer, "rb")

def get_position(samfile):
    coverage_data = []
    chrom_names = samfile.references  

    for chrom in chrom_names:
        chrom_length = samfile.get_reference_length(chrom)
        depths = [0] * chrom_length
        
        for pileup_column in samfile.pileup(chrom):
            depths[pileup_column.reference_pos] = pileup_column.nsegments
        
        coverage_data.extend(depths)
    return coverage_data

coverage = get_position(samfile)
coverage_without_primer = get_position(samfile_primer)

# Plot the combined coverage: before and after trimming
plt.figure(figsize=(15, 6))
plt.plot(range(len(coverage)), coverage, color='blue')
plt.plot(range(len(coverage_without_primer)), coverage_without_primer, color='red')
plt.title(f"Depth Coverage Plot {name}: before and after trimming")
plt.xlabel("Genomic Position")
plt.ylabel("Depth")
plt.show()
plt.savefig(f"result/artic/Bysubtypes/{subtype}/figures/figure_compare/coverage_primer_trimming_{name}.png")