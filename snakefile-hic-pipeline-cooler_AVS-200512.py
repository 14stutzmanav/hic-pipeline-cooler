import pandas as pd
import os
import re
import glob
import getpass

configfile: './slurmConfig.json'

# Set Genome and Bin Size:
GenomeAssembly = 'dm6'
binsize = '1'
chrsizePATH = str('/proj/mckaylab/users/astutzman/hic-analyses/juicer/AS1-21_juicer-output/dm6.chrom.sizes')

########################

# Set Module Versions:

bwaVer = str('bwa/0.7.17')
samtoolsVer = str('samtools/1.9')
pairtoolsVer = str('pairtools/0.3.0')
coolerVer = str('cooler/0.8.6')

#########################

# Load Required Input Files:

BWAIndexPATH = str('/proj/seq/data/' + GenomeAssembly + '_UCSC/Sequence/')
sampleSheetPath = str('/proj/mckaylab/users/astutzman/samplesheet-hic.csv')
sampleDF = pd.read_csv(sampleSheetPath, comment = '#')

techList = list(set(sampleDF.techName))
fastqList= list(set(sampleDF.fastq))
baseNameList= list(set(sampleDF.baseName))
readNumList = ['R1', 'R2']

#########################
localrules: all

rule all:
	input:
		expand("Fastq/{fastq}_{readNum}.fastq.gz", fastq=fastqList, readNum=readNumList),
		#expand("Fastq/{base}_{readNum}.fastq.gz", base=baseNameList, readNum=readNumList),
		expand("Bam/{tech}.bam", tech=techList),
		expand("Pairs/{tech}_parsed.pairsam", tech=techList),
		expand("Pairs/{tech}_sorted.pairsam", tech=techList),
		expand("Pairs/{tech}_dedup.pairsam", tech=techList),
		expand("Pairs/{tech}.pairs", tech=techList),
		expand("Cool/{tech}.cool", tech=techList)

rule copyFiles:
	input:
		lambda x: list(sampleDF.htsfFile)
	output:
		expand("Fastq/{fastq}_{readNum}.fastq.gz", fastq=list(sampleDF.fastq), readNum=readNumList)
	run:
		for htsf in list(sampleDF.htsfFile):
			outFileFilt = sampleDF [ sampleDF.htsfFile == htsf ] 
			outFileBase = list(outFileFilt.fastq)[0]
			outFile = 'Fastq/{fastq}.fastq.gz'.format(fastq = outFileBase)
			#print('THERE')
			#print(outFile)
			shutil.copyfile(htsf, outFile)
			print('copied file')


rule map:
	input:
		#fastq_R1 = "/proj/mckaylab/users/astutzman/hic-analyses/hic-pipeline-cooler/yw-3LW-wing-HiC/{sample}_R1.fastq.gz",
		#fastq_R2 = "/proj/mckaylab/users/astutzman/hic-analyses/hic-pipeline-cooler/yw-3LW-wing-HiC/{sample}_R2.fastq.gz",
		fastq_R1 = "Fastq/{fastq}_R1.fastq.gz",
		fastq_R2 = "Fastq/{fastq}_R2.fastq.gz",
		BWAIndex = BWAIndexPATH
	
	output:
		bam = "Bam/{sample}.bam"
		
	params:
		module1 = bwaVer,
		module2 = samtoolsVer
	
	shell:
		"""
		module purge && module load {params.module1} && module load {params.module2}
		bwa mem -SP5M {input.BWAIndex} {input.fastq_R1} {input.fastq_R2} | samtools view -bhS - > {output.bam}
		"""	


rule filter:
	input:
		bam = "Bam/{sample}.bam",
		chromSize = chrsizePATH
	output:
		parsed = "Pairs/{sample}_parsed.pairsam",
		sorte = "Pairs/{sample}_sorted.pairsam",
		dedup = "Pairs/{sample}_dedup.pairsam",
		filtered = "Pairs/{sample}_filtered.pairsam",
		pairs = "Pairs/{sample}.pairs"
	params:
		moduleVer = pairtoolsVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		samtools view -h {input.bam} | \ pairtools parse -c {input.chromSize} -o {output.parsed}
		pairtools sort --nproc 8 -o {output.sorte} {output.parsed}
		pairtools dedup --mark-dups -o {output.dedup} {output.sorte}
		pairtools select \ '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU)' \ -o {output.filtered} {output.dedup}
		pairtools split --output-pairs {output.pairs} {output.filtered}
		"""

rule bin:
	input:
		binsize = "{binsize}",
		chrsize = "{chrsizePATH}",
		pairs = "Pairs/{sample}.pairs"
	output:
		cool = "Cool/{sample}.cool"
	params:
		moduleVer = coolerVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		cooler cload pairix {input.chrsize}:{input.binsize} {input.pairs} {output.cool}
		cooler balance {output.cool}
		"""

rule aggregate:
	input:
		cool="Cool/{sample}.cool"
	output:
		mcool="Cool/{sample}.mcool"
	params:
		moduleVer=coolerVer
	shell:
		"""
		module purge && module load {params.moduleVer}
		cooler zoomify {cool.input}
		"""



