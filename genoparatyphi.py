#!/usr/bin/env python
#
# Input BAM (recommended) or VCF (if highly trusted SNP data) relative to Salmonella Paratyphi A AKU_12601 (FM200053) for detection of amr mutations.
#
# Authors - Zoe Anne Dyson (zad24@medschl.cam.ac.uk), Kat Holt (drkatholt@gmail.com)
#
# Documentation - https://github.com/katholt/genoparatyphi
#
# Dependencies:
#	 SAMtools (v1.2) and bcftools (v1.2) are required to genotype from BAMs.
#
# Last modified - Feb 6, 2020
#

from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, copy
import gzip
import logging
import time
import datetime
from subprocess import call, check_output, CalledProcessError, STDOUT
from Bio.Seq import Seq


def parse_args():
	"Parse the input arguments, use '-h' for help"
	parser = ArgumentParser(description='VCF to Paratyphi A AMR mutations')
	parser.add_argument('--mode', required=True,
						help='Mode to run in based on input files (vcf, bam, or vcf_parsnp)')
	parser.add_argument(
		'--vcf', nargs='+', type=str, required=False,
		help='VCF file(s) to genotype (Mapping MUST have been done using AKU_12601 as a reference sequence)')
	parser.add_argument('--bam', nargs='+', type=str, required=False,
						help='BAM file(s) to genotype (Mapping MUST have been done using AKU_12601 as a reference sequence')
	parser.add_argument('--ref_id', type=str, required=False,
						help='Name of the reference in the VCF file (#CHROM column)')
	parser.add_argument('--phred', type=int, required=False, default=20,
						help='Minimum phred quality to count a variant call vs AKU_12601 as a true SNP (default 20)')
	parser.add_argument('--min_prop', type=float, required=False, default=0.1,
						help='Minimum proportion of reads required to call a SNP (default 0.1)')
	parser.add_argument('--ref', type=str, required=False,
						help='Reference sequence in fasta format. Required if bam files provided.')
	parser.add_argument('--output', type=str, required=False, default='amr_mutations.txt',
						help='Location and name for output file.')
	parser.add_argument('--samtools_location', type=str, required=False, default='',
						help='Location of folder containing samtools installation if not standard/in path.')
	parser.add_argument('--bcftools_location', type=str, required=False, default='',
						help='Location of folder containing bcftools installation if not standard/in path.')
	return parser.parse_args()



# check if this SNP defines a amr group
def checkamrSNP(vcf_line_split, strain_amr_bases, amr_loci_list, amr_loci, args):
	# get the snp position
	amr_snp = int(vcf_line_split[1])
	# if the snp position is in our dict of amr snp loci, parse
	if amr_snp in amr_loci_list:
		# if our snp passes the phred cutoff, continue
		if float(vcf_line_split[5]) > args.phred:
			# get the dp4 read counts
			m = re.search("DP4=(\d+),(\d+),(\d+),(\d+)", vcf_line_split[7])
			# calculate proportion of reads that support this snp position
			if m != None:
				alt_read_count = int(m.group(3)) + int(m.group(4))
				total_read_count = alt_read_count + int(m.group(1)) + int(m.group(2))
				if float(total_read_count) == 0:
					amr_snp_proportion = float(-1)
				else:
					amr_snp_proportion = float(alt_read_count) / total_read_count
			# if no dp4, then calculate proportion of reads supporting snp a different way
			else:
				if vcf_line_split[4] != '.':  # if the ALT is not '.' i.e. if the alt is not same as ref
					try:
						ad = vcf_line_split[9].split(':')[1].split(',')	 # get the AD ratio
						alt_read_count = int(ad[1])
						total_read_count = int(ad[0]) + alt_read_count
						amr_snp_proportion = float(alt_read_count) / total_read_count
					except IndexError:
						amr_snp_proportion = float(-1)
				# otherwise we don't know the proportion of reads supporting snp
				else:
					amr_snp_proportion = float(
						-1)	 # set unknowns to negative so that we know this is not a real proportion

			# get the actually snp base call in the strain
			amr_snp_allele = vcf_line_split[4]

			# work out what codon this snp position is in
			# get the codon number of the snp (is it the 1st, 2nd or 3rd base in the codon?)
			for loci in amr_loci:
				if amr_snp in amr_loci[loci]:
					loci_of_interest = loci
					base_num = amr_loci[loci].index(amr_snp)
			# if our snp meets the proportion threshold for inclusion, then we alter our amr dict list at the appropriate place
			if amr_snp_proportion > args.min_prop:
				strain_amr_bases[loci_of_interest][base_num] = amr_snp_allele


	return (strain_amr_bases)


# exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
	pass


def run_command(command, **kwargs):
	'Execute a shell command and check the exit status and any O/S exceptions'
	command_str = ' '.join(command)
	logging.info('Running: {}'.format(command_str))
	try:
		exit_status = call(command, **kwargs)
	except OSError as e:
		message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
		raise CommandError({"message": message})
	if exit_status != 0:
		message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
		raise CommandError({"message": message})


# main function
def main():
	args = parse_args()

	### AMR SNP definitions

	amr_wt_aa = {'gyrA-87': 'D', 'gyrA-83': 'S', 'parC-80': 'S','parC-84': 'E', 'acrB-717':'R'}
	# base positions for each codon, in order of base 1, 2, 3 on the FORWARD STRAND
	# note that the read codons and bases for parC are on the reverse complement strand
	amr_loci = {'gyrA-87': [669213, 669214, 669213], 'gyrA-83': [669201,669202,669203],'parC-80': [3145220,3145221,3145222],'parC-84': [3145208,3145209,3145210], 'acrB-717': [2338395,2338396,2338397]}
	amr_loci_list = [669213, 669214, 669213, 669201,669202,669203, 3145220,3145221,3145222, 3145208,3145209,3145210,2338395,2338396,2338397]
	amr_ref_alleles = {'gyrA-87': ['G','A','C'], 'gyrA-83': ['T','C','C'], 'parC-80': ['G','C','T'], 'parC-84': ['T','T','C'], 'acrB-717': ['C','G','A']}

	if (((args.mode == 'vcf') and args.vcf and args.ref_id) or (
			(args.mode == 'bam') and args.bam and args.ref and args.ref_id) or (
			(args.mode == 'vcf_parsnp') and args.vcf)):

		# Initialise output file and timestamp
		timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d_%H.%M.%S_')
		output_file = open(timestamp + args.output, 'w')

		# GENERATE VCFS (1 per strain) FROM BAMS

		if args.mode == 'bam':

			with open(args.ref, 'r') as fasta_file:	 # create SAMtools compatible fasta file
				#sequence = fasta_file.read()
				sequences = fasta_file.read()
				for sequence in sequences.split('>'):
					print("got to here 1")
					if args.ref_id in sequence:
						print("got to here 2")
						new_header = '>' + str(args.ref_id)
						replacement_index = sequence.find('\n')
						#print("replacement index is" + str(replacement_index)
						with open('temp_reference.fasta', 'w') as temp_fasta_file:
							temp_fasta_file.write(new_header + sequence[replacement_index:])
							print("got to here 3")


			vcfFiles = []

			# coordinates in zero-base, half-open for SAMtools compatible bed file
			ordered_loci = list(amr_loci_list)
			sorted(ordered_loci)
			temp_bed_file = open(args.ref_id + '.bed', 'w')	 # create temporary bed file for SAMtools
			for locus in ordered_loci:
				temp_bed_file.write(
					args.ref_id + '\t' + str(locus - 1) + '\t' + str(locus) + '\n')	 # write bed file from matrix
			temp_bed_file.close()  # close bedFile

			run_command(['samtools', 'faidx', 'temp_reference.fasta'])	# index fasta file

			for bam in args.bam:
				print('bam files supplied, generating vcf file for ' + bam)
				if os.path.exists(bam + '.bai') == False:  # index bam file if indexed bam not provided
					run_command([args.samtools_location + 'samtools', 'index', bam])

				run_command(
					[args.samtools_location + 'samtools', 'mpileup', '-q', str(args.phred), '-ugB', '-f',
					 'temp_reference.fasta',
					 '-l', args.ref_id + '.bed', bam, '-o', bam[:-4] + '.output', '-I'])  # detect SNPs

				run_command(
					[args.bcftools_location + 'bcftools', 'call', '-c', bam[:-4] + '.output', '-o',
					 bam[:-4] + '.vcf'])  # generate vcf files
				run_command(['rm', bam[:-4] + '.output'])

				vcfFiles.append(bam[:-4] + '.vcf')	# supply generated vcf file to script

			run_command(['rm', args.ref_id + '.bed'])  # remove temp files
			run_command(['rm', 'temp_reference.fasta'])
			run_command(['rm', 'temp_reference.fasta.fai'])

			args.vcf = vcfFiles

		# PRINT OUTPUT HEADER

		if args.mode == 'bam':
			output_file.write('\t'.join(
				['File', 'AMR mutations', '\n']))
		elif args.mode == 'vcf':
			output_file.write('\t'.join(
				['File', 'AMR mutations', '\n']))
		else:
			output_file.write('\t'.join(
				['File', 'AMR mutations', '\n']))

		# PARSE MAPPING BASED VCFS (1 per strain)

		if (args.vcf and (args.mode != 'vcf_parsnp')):
			for vcf in args.vcf:
				this_amr_groups = []  # list of amr SNPs found
				amr_proportions = {}  # proportion of reads supporting each defining SNP; key = group, value = proportion
				# initailse the amr bases to be the same as the reference
				strain_amr_bases = copy.deepcopy(amr_ref_alleles)

				# read file
				(file_name, ext) = os.path.splitext(vcf)

				if ext == '.gz':
					f = gzip.open(vcf, 'r')
				else:
					f = open(vcf, 'r')

				any_ref_line = 0

				for line in f:
					if not line.startswith('#'):
						x = line.rstrip().split()
						if x[0] == args.ref_id:
							# parse this SNP line
							any_ref_line = 1
							strain_amr_bases = checkamrSNP(x, strain_amr_bases, amr_loci_list, amr_loci, args)

				f.close()

				# work out what the amr mutations are by looking at the codons that have been generated, and comparing
				# these to WT
				# in this instance our codons need to be reverse complemented before they can be compared
				for amr_codon in strain_amr_bases:
					strain_codon = Seq(''.join(strain_amr_bases[amr_codon]))
					# reverse complement it, then translate
					print(str(amr_codon))
					if str(amr_codon)=="parC-80" or str(amr_codon)=="parC-84":
						strain_codon = strain_codon.reverse_complement().translate()
					else:
						strain_codon = strain_codon.translate()
					# compare it to the WT aa for that codon, store if it's NOT a match
					if strain_codon != amr_wt_aa[amr_codon]:
						gene_name = amr_codon.split('-')[0]
						codon_num = amr_codon.split('-')[1]
						mutation = gene_name + '-' + amr_wt_aa[amr_codon] + codon_num + str(strain_codon)
						this_amr_groups.append(mutation)

				if any_ref_line > 0:
					if args.bam:
						output_file.write(
							vcf + '\t' + ','.join(this_amr_groups) + '\t' + '\n')


					else:
						output_file.write(vcf + '\t' + ','.join(this_amr_groups) + '\n')

				else:
					output_file.write(
						vcf + '\tNo SNPs encountered against expected reference. Wrong reference or no SNP calls?\n')

		# PARSE PARSNP VCF (multiple strains)

		if args.mode == 'vcf_parsnp':

			if not args.ref_id:
				args.ref_id = '1'

			for vcfm in args.vcf:

				# read file
				(file_name, ext) = os.path.splitext(vcfm)

				if ext == '.gz':
					f = gzip.open(vcfm, 'r')
				else:
					f = open(vcfm, 'r')

				any_ref_line = 0

				strains = []

				for line in f:
					x = line.rstrip().split()
					if x[0] == '#CHROM':
						strains = x[10:]
					if not line.startswith('#'):
						if x[0] == args.ref_id:
							any_ref_line = 1  # parse this SNP line

				f.close()

				# collate by strain
				if any_ref_line > 0:
					for strain in this_groups:
						output_file.write(strains[strain] + '\n')

				else:
					output_file.write(strains[strain] + '\tNo SNPs encountered against expected reference. Wrong reference or no SNP calls?\n')

		output_file.close()
	else:
		print('Missing or incomplete input parameters, please check these and try again.')


# call main function
if __name__ == '__main__':
	main()
