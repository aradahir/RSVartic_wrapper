import pandas as pd
import os, glob, snakemake
from collections import Counter

depth = 40 #change this line for minimum depth coverage plot
df = pd.read_csv("samplesheet.csv", header=None) 

for index, row in df.iterrows():
	if index == 0:  
		continue 
	else:
		os.makedirs("data/concatenated_files/", exist_ok=True)
		os.system('cat data/raw/%s/*.fastq.gz > data/concatenated_files/%s.fastq.gz' %(row[0],  row[1]))
		os.makedirs("data/fasta_blastn/", exist_ok=True)
		os.system("seqkit fq2fa data/concatenated_files/%s.fastq.gz -o data/fasta_blastn/%s.fasta" %(row[1], row[1]))
		os.system("blastn -query data/fasta_blastn/%s.fasta -db resources/db/RSV_db -out data/fasta_blastn/%s_blastn.txt -outfmt 6 " %(row[1],row[1]))

		file_path = "data/fasta_blastn/" + row[1] + "_blastn.txt"

		with open(file_path) as f:
			subtypes = [line.split('\t')[1] for line in f if line.strip()]

		most_common = []
		count = Counter(subtypes)
		most_common = count.most_common(1)

		if most_common:
			subtype, frequency = most_common[0]
			os.makedirs(os.path.join("result/data/", subtype), exist_ok=True)
			os.makedirs(os.path.join("result/log/"), exist_ok=True)
			os.makedirs(os.path.join("result/artic/"), exist_ok=True)
			os.makedirs(os.path.join("result/F_proteins/"), exist_ok=True)
			os.makedirs(os.path.join("result/status/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/artic/Bysubtypes/{subtype}/figures/figure_compare/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/artic/Bysubtypes/{subtype}/figures/figure_annotation/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/artic/Bysubtypes/{subtype}/consensus/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/artic/Bysubtypes/{subtype}/artic_consensus/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/F_proteins/{subtype}/mutation/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/F_proteins/{subtype}/consensus/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/artic/Bysubtypes/{subtype}/variants_annotation/"), exist_ok=True)
			os.makedirs(os.path.join(f"result/artic/Bysubtypes/{subtype}/table/"), exist_ok=True)
			os.system("cp data/concatenated_files/%s.fastq.gz result/data/%s/%s.fastq.gz" %(row[1], subtype, row[1]))

		else:
			os.makedirs(os.path.join("result/data/", "unknown"), exist_ok=True)
			print(f"Error: Cannot determine subtype")

SAMPLE= glob_wildcards('result/data/{sub}/{sample}.fastq.gz')

rule all:
	input:
		expand('result/status/{sub}_{samples}_rearrange.txt', zip , sub = SAMPLE.sub, samples = SAMPLE.sample),
					
#change the artic command for artic 1.5
rule artic:
	input:
		'result/data/{sub}/{samples}.fastq.gz'
	output:
		'result/status/{sub}_{samples}_artic.txt'
	params:
		scheme = '{sub}/V1',
		output = 'result/artic/{sub}_{samples}_polished',
		log = 'result/log/{sub}_{samples}_artic.txt', 
		bed = 'resources/{sub}/V1/{sub}.scheme.bed',
		ref = 'resources/{sub}/V1/{sub}.reference.fasta'
	message:
		'Primer trimming using Artic on {wildcards.samples}..'
	shell:'''
		artic minion --normalise 0 --threads 45 \
		--bed {params.bed} \
		--ref {params.ref} \
		--read-file {input} \
		{params.output} 2>> {params.log}

		touch {output}
		'''

rule depth:
	input:
		'result/status/{sub}_{samples}_artic.txt'
	output:
		'result/status/{sub}_{samples}_depth.txt'
	params:
		name = '{sub}_{samples}',
		subtype = '{sub}',
		min_depth = depth,
		variant_annotate = 'result/artic/Bysubtypes/{sub}/variants_annotation/{sub}_{samples}_variants.txt'
	shell:'''
		python resources/script/depth_compare.py -n {params.name} -t {params.subtype}
		python resources/script/depth_annotated.py -n {params.name} -t {params.subtype} -m {params.min_depth}
		python resources/script/mutation_analysis.py -n {params.name} -t {params.subtype} -s F -m {params.min_depth}
		python resources/script/mutation_analysis.py -n {params.name} -t {params.subtype} -s all -m {params.min_depth}
		bedtools intersect -a result/artic/{params.name}_polished.pass.vcf -b resources/mutation/{params.subtype}.gtf -wa -wb > {params.variant_annotate}
		touch {output}	
		'''

rule rearrange:
	input:
		'result/status/{sub}_{samples}_depth.txt'
	output:
		'result/status/{sub}_{samples}_rearrange.txt'
	params:
		name = '{samples}',
		subtype = '{sub}', 
		folder = 'result/artic/{samples}'
	shell:'''
		mkdir -p {params.folder}
		mv result/artic/{params.subtype}_{params.name}* {params.folder}
		touch {output}	
		'''

#create the summary report and rearrange files
onsuccess:
	shell("python resources/script/summary.py")


