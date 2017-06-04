#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
patho_typing.py - In silico pathogenic typing directly from raw
Illumina reads
<https://github.com/B-UMMI/patho_typing/>

Copyright (C) 2017 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: June 04, 2017

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import time
import sys

import modules.utils as utils
import modules.run_rematch as run_rematch
import modules.typing as typing

version = '0.2'


def set_reference(species, outdir, script_path, trueCoverage):
	reference_file = None
	trueCoverage_file = None
	trueCoverage_sequences = None
	trueCoverage_headers = None
	trueCoverage_config = None
	typing_file = None
	typing_sequences = None
	typing_headers = None
	typing_rules = None
	typing_config = None

	species_folder = os.path.join(os.path.dirname(script_path), 'modules', 'seq_conf', '_'.join([s.lower() for s in species]), '')

	if os.path.isdir(species_folder):
		typing_rules = os.path.join(species_folder, 'typing_rules.tab')
		typing_file = os.path.join(species_folder, 'typing.fasta')
		typing_sequences, ignore = utils.get_sequence_information(typing_file, 0)
		typing_sequences, typing_headers = utils.clean_headers_sequences(typing_sequences)
		typing_sequences = utils.simplify_sequence_dict(typing_sequences)
		typing_config = os.path.join(species_folder, 'typing.config')
		if trueCoverage:
			trueCoverage_file = os.path.join(species_folder, 'trueCoverage.fasta')
			trueCoverage_sequences, ignore = utils.get_sequence_information(trueCoverage_file, 0)
			trueCoverage_sequences, trueCoverage_headers = utils.clean_headers_sequences(trueCoverage_sequences)
			trueCoverage_sequences = utils.simplify_sequence_dict(trueCoverage_sequences)
			trueCoverage_config = os.path.join(species_folder, 'trueCoverage.config')

			trueCoverage_typing_sequences = trueCoverage_sequences.copy()
			for header in typing_sequences:
				if header not in trueCoverage_sequences:
					trueCoverage_typing_sequences[header] = typing_sequences[header]
				else:
					print 'Sequence {header} of typing.fasta already present in trueCoverage.fasta'.format(header=header)

			reference_file = os.path.join(outdir, 'trueCoverage_typing.fasta')
			write_sequeces(reference_file, trueCoverage_typing_sequences)
		else:
			reference_file = os.path.join(outdir, 'typing.fasta')
			write_sequeces(reference_file, typing_sequences)
	else:
		species_present = []
		seq_conf_folder = os.path.join(os.path.dirname(script_path), 'modules', 'seq_conf', '')
		species_folder = [d for d in os.listdir(seq_conf_folder) if not d.startswith('.') and os.path.isdir(os.path.join(seq_conf_folder, d, ''))]
		for species in species_folder:
			species = species.split('_')
			species[0] = species[0][0].upper() + species[0][1:]
			species_present.append(' '.join(species))
		sys.exit('Only these species are available:' + '\n' + '\n'.join(species_present))

	return reference_file, trueCoverage_file, trueCoverage_sequences, trueCoverage_headers, trueCoverage_config, typing_file, typing_sequences, typing_headers, typing_rules, typing_config


def confirm_genes_fasta_rules(typing_fasta_headers, typing_rules_file):
	with open(typing_rules_file, 'rtU') as reader:
		genes = []
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				line = line.split('\t')
				if line[0].startswith('#'):
					genes = line[1:]
					break
		missing_genes = []
		for gene in genes:
			if gene.lower() not in typing_fasta_headers:
				missing_genes.append(gene)
		if len(missing_genes) > 0:
			sys.exit('The following genes required for typing are not present in typing fasta file: {missing_genes}'.format(missing_genes=', '.join(missing_genes)))


def index_fasta_samtools(fasta, region_None, region_outfile_none, print_comand_True):
	command = ['samtools', 'faidx', fasta, '', '', '']
	shell_true = False
	if region_None is not None:
		command[3] = region_None
	if region_outfile_none is not None:
		command[4] = '>'
		command[5] = region_outfile_none
		shell_true = True
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, shell_true, None, print_comand_True)
	return run_successfully, stdout


def indexSequenceBowtie2(referenceFile, threads):
	if os.path.isfile(str(referenceFile + '.1.bt2')):
		run_successfully = True
	else:
		command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


def run_bowtie(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc):
	sam_file = os.path.join(outdir, str('alignment.sam'))

	run_successfully = indexSequenceBowtie2(referenceFile, threads)
	if run_successfully:
		command = ['bowtie2', '-k', str(numMapLoc), '-q', '', '--threads', str(threads), '-x', referenceFile, '', '--no-unal', '-S', sam_file]

		if len(fastq_files) == 1:
			command[9] = '-U ' + fastq_files[0]
		elif len(fastq_files) == 2:
			command[9] = '-1 ' + fastq_files[0] + ' -2 ' + fastq_files[1]
		else:
			return False, None

		if conserved_True:
			command[4] = '--sensitive'
		else:
			command[4] = '--very-sensitive-local'

		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

	if not run_successfully:
		sam_file = None

	return run_successfully, sam_file


def sortAlignment(alignment_file, output_file, sortByName_True, threads):
	outFormat_string = os.path.splitext(output_file)[1][1:].lower()
	command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
	if sortByName_True:
		command[6] = '-n'
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	if not run_successfully:
		output_file = None
	return run_successfully, output_file


def indexAlignment(alignment_file):
	command = ['samtools', 'index', alignment_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


def mapping_reads(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc):
	print '\n' + 'Mapping the reads' + '\n'
	run_successfully, sam_file = run_bowtie(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc)
	if run_successfully:
		run_successfully, bam_file = sortAlignment(sam_file, str(os.path.splitext(sam_file)[0] + '.bam'), False, threads)
		if run_successfully:
			os.remove(sam_file)
			run_successfully = indexAlignment(bam_file)
			if run_successfully:
				index_fasta_samtools(referenceFile, None, None, True)
	return run_successfully, bam_file


def include_rematch_dependencies_path():
	command = ['which', 'rematch.py']
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
	if run_successfully:
		rematch = stdout.splitlines()[0]
		utils.setPATHvariable(False, rematch)
	return rematch


def split_bam(bam_file, list_sequences, outdir, threads):
	new_bam = os.path.join(outdir, 'partial.bam')
	command = ['samtools', 'view', '-b', '-u', '-h', '-o', new_bam, '-@', str(threads), bam_file, ' '.join(list_sequences)]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully, new_bam


def parse_config(config_file):
	config = {'reference_file': None, 'length_extra_seq': None, 'maximum_number_absent_genes': None, 'maximum_number_genes_multiple_alleles': None, 'minimum_read_coverage': None, 'minimum_depth_presence': None, 'minimum_depth_call': None, 'minimum_depth_frequency_dominant_allele': None, 'minimum_gene_coverage': None, 'minimum_gene_identity': None}

	with open(config_file, 'rtU') as reader:
		field = None
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				line = line.split(' ')[0]
				if line.startswith('#'):
					line = line[1:].split(' ')[0]
					field = line
				else:
					if field is not None:
						if field in ['length_extra_seq', 'maximum_number_absent_genes', 'maximum_number_genes_multiple_alleles', 'minimum_read_coverage', 'minimum_depth_presence', 'minimum_depth_call', 'minimum_gene_coverage', 'minimum_gene_identity']:
							line = int(line)
							if field in ['minimum_gene_coverage', 'minimum_gene_identity']:
								if line < 0 or line > 100:
									sys.exit('minimum_gene_coverage in trueCoverage_rematch config file must be an integer between 0 and 100')
						elif field == 'minimum_depth_frequency_dominant_allele':
							line = float(line)
							if line < 0 or line > 1:
								sys.exit('minimum_depth_frequency_dominant_allele in trueCoverage_rematch config file must be a double between 0 and 1')
						config[field] = line
						field = None

	for field in config:
		if config[field] is None:
			sys.exit(field + ' in trueCoverage_rematch config file is missing')

	return config


def clean_pathotyping_folder(outdir, reference_file, debug_mode_true):
	if not debug_mode_true:
		files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
		for file_found in files:
			if file_found.startswith(('alignment.', os.path.splitext(os.path.basename(reference_file))[0])):
				file_found = os.path.join(outdir, file_found)
				os.remove(file_found)


def write_sequeces(out_file, sequences_dict):
	with open(out_file, 'wt') as writer:
		for header in sequences_dict:
			writer.write('>' + header + '\n')
			writer.write('\n'.join(utils.chunkstring(sequences_dict[header]['sequence'], 80)) + '\n')


def main():
	parser = argparse.ArgumentParser(prog='patho_typing.py', description='In silico pathogenic typing directly from raw Illumina reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_required = parser.add_argument_group('Required options')
	parser_required.add_argument('-f', '--fastq', nargs='+', action=utils.required_length((1, 2), '--fastq'), type=argparse.FileType('r'), metavar=('/path/to/input/file.fq.gz'), help='Path to single OR paired-end fastq files. If two files are passed, they will be assumed as being the paired fastq files', required=True)
	parser_required.add_argument('-s', '--species', nargs=2, type=str, metavar=('Yersinia', 'enterocolitica'), help='Species name', required=True)

	parser_optional_general = parser.add_argument_group('General facultative options')
	parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/', help='Path to the directory where the information will be stored', required=False, default='.')
	parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use', required=False, default=1)
	parser_optional_general.add_argument('--trueCoverage', action='store_true', help='Assess true coverage before continue typing')
	parser_optional_general.add_argument('--noCheckPoint', action='store_true', help='DeBug Mode: do not remove temporary files')
	parser_optional_general.add_argument('--minGeneCoverage', type=int, metavar='N', help='Minimum typing percentage of target reference gene sequence covered to consider a gene to be present (value between [0, 100])', required=False)
	parser_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N', help='Minimum typing percentage of identity of reference gene sequence covered to consider a gene to be present (value between [0, 100]). One INDEL will be considered as one difference', required=False)
	parser_optional_general.add_argument('--minGeneDepth', type=int, metavar='N', help='Minimum typing gene average coverage depth of present positions to consider a gene to be present (default 15, or 1/3 of average sample coverage assessed by true coverage analysis)', required=False)
	parser_optional_general.add_argument('--doNotRemoveConsensus', action='store_true', help='Do not remove ReMatCh consensus sequences')
	parser_optional_general.add_argument('--debug', action='store_true', help='DeBug Mode: do not remove temporary files')

	args = parser.parse_args()

	if args.minGeneCoverage is not None and (args.minGeneCoverage < 0 or args.minGeneCoverage > 100):
		parser.error('--minGeneCoverage should be a value between [0, 100]')
	if args.minGeneIdentity is not None and (args.minGeneIdentity < 0 or args.minGeneIdentity > 100):
		parser.error('--minGeneIdentity should be a value between [0, 100]')

	start_time = time.time()

	args.outdir = os.path.abspath(args.outdir)
	if not os.path.isdir(args.outdir):
		os.makedirs(args.outdir)

	# Start logger
	logfile, time_str = utils.start_logger(args.outdir)

	script_path = utils.general_information(logfile, version, args.outdir, time_str)
	print '\n'

	rematch = include_rematch_dependencies_path()

	args.fastq = [fastq.name for fastq in args.fastq]

	reference_file, trueCoverage_file, trueCoverage_sequences, trueCoverage_headers, trueCoverage_config, typing_file, typing_sequences, typing_headers, typing_rules, typing_config = set_reference(args.species, args.outdir, script_path, args.trueCoverage)
	original_reference_file = str(reference_file)

	confirm_genes_fasta_rules(typing_headers, typing_rules)

	run_successfully, bam_file = mapping_reads(args.fastq, reference_file, args.threads, args.outdir, False, 1)
	if run_successfully:
		rematch_dir = os.path.join(args.outdir, 'rematch', '')
		if not os.path.isdir(rematch_dir):
			os.makedirs(rematch_dir)

		if args.trueCoverage:
			trueCoverage_dir = os.path.join(rematch_dir, 'trueCoverage', '')
			if not os.path.isdir(trueCoverage_dir):
				os.makedirs(trueCoverage_dir)

			print '\n'
			run_successfully, trueCoverage_bam = split_bam(bam_file, trueCoverage_headers, trueCoverage_dir, args.threads)
			if run_successfully:
				run_successfully = indexAlignment(trueCoverage_bam)
				if run_successfully:
					reference_file = os.path.join(trueCoverage_dir, 'reference.fasta')
					write_sequeces(reference_file, trueCoverage_sequences)
					index_fasta_samtools(reference_file, None, None, True)
					config = parse_config(trueCoverage_config)
					runtime, run_successfully, sample_data_general, data_by_gene = run_rematch.run_rematch(rematch, trueCoverage_dir, reference_file, trueCoverage_bam, args.threads, config['length_extra_seq'], config['minimum_depth_presence'], config['minimum_depth_call'], config['minimum_depth_frequency_dominant_allele'], config['minimum_gene_coverage'], config['minimum_gene_identity'], args.debug, args.doNotRemoveConsensus)

					if run_successfully and sample_data_general['mean_sample_coverage'] is not None and sample_data_general['number_absent_genes'] is not None and sample_data_general['number_genes_multiple_alleles'] is not None:
						if args.minGeneDepth is None:
							args.minGeneDepth = sample_data_general['mean_sample_coverage'] / 3

						exit_info = []
						if sample_data_general['mean_sample_coverage'] < config['minimum_read_coverage']:
							exit_info.append('Sample coverage ({mean_sample_coverage}) lower than the minimum required ({minimum_read_coverage})'.format(mean_sample_coverage=sample_data_general['mean_sample_coverage'], minimum_read_coverage=config['minimum_read_coverage']))
						if sample_data_general['number_absent_genes'] > config['maximum_number_absent_genes']:
							exit_info.append('Number of absent genes ({number_absent_genes}) higher than the maximum allowed ({maximum_number_absent_genes})'.format(number_absent_genes=sample_data_general['number_absent_genes'], maximum_number_absent_genes=config['maximum_number_absent_genes']))
						if sample_data_general['number_genes_multiple_alleles'] > config['maximum_number_genes_multiple_alleles']:
							exit_info.append('Number of genes with multiple alleles ({number_genes_multiple_alleles}) higher than the maximum allowed ({maximum_number_genes_multiple_alleles})'.format(number_genes_multiple_alleles=sample_data_general['number_genes_multiple_alleles'], maximum_number_genes_multiple_alleles=config['maximum_number_genes_multiple_alleles']))

						if len(exit_info) > 0:
							print '\n' + '\n'.join(exit_info) + '\n'
							e = 'TrueCoverage requirements not fullfill'
							print '\n' + e + '\n'
							if not args.noCheckPoint:
								clean_pathotyping_folder(args.outdir, original_reference_file, args.debug)
								time_taken = utils.runTime(start_time)
								sys.exit(e)
					else:
						e = 'TrueCoverage module did not run successfully'
						print '\n' + e + '\n'
						if not args.noCheckPoint:
							clean_pathotyping_folder(args.outdir, original_reference_file, args.debug)
							time_taken = utils.runTime(start_time)
							sys.exit(e)

					print '\n'
					typing_dir = os.path.join(rematch_dir, 'typing', '')
					if not os.path.isdir(typing_dir):
						os.makedirs(typing_dir)
					run_successfully, bam_file = split_bam(bam_file, typing_headers, typing_dir, args.threads)
					if run_successfully:
						run_successfully = indexAlignment(bam_file)
						if run_successfully:
							reference_file = os.path.join(typing_dir, 'reference.fasta')
							write_sequeces(reference_file, typing_sequences)
							index_fasta_samtools(reference_file, None, None, True)
							rematch_dir = str(typing_dir)
			if not run_successfully:
				if args.noCheckPoint:
					clean_pathotyping_folder(args.outdir, original_reference_file, args.debug)
					time_taken = utils.runTime(start_time)
					sys.exit('Something in the required TrueCoverage analysis went wrong')

		if run_successfully:
			config = parse_config(typing_config)
			if args.minGeneCoverage is not None:
				config['minimum_gene_coverage'] = args.minGeneCoverage
			if args.minGeneIdentity is not None:
				config['minimum_gene_identity'] = args.minGeneIdentity

			runtime, run_successfully, sample_data_general, data_by_gene = run_rematch.run_rematch(rematch, rematch_dir, reference_file, bam_file, args.threads, config['length_extra_seq'], config['minimum_depth_presence'], config['minimum_depth_call'], config['minimum_depth_frequency_dominant_allele'], config['minimum_gene_coverage'], config['minimum_gene_identity'], args.debug, args.doNotRemoveConsensus)
			if run_successfully and data_by_gene is not None:
				if args.minGeneDepth is None:
					args.minGeneDepth = 15

				runtime, ignore, ignore = typing.typing(data_by_gene, typing_rules, config['minimum_gene_coverage'], config['minimum_gene_identity'], args.minGeneDepth, args.outdir)
			else:
				clean_pathotyping_folder(args.outdir, original_reference_file, args.debug)
				time_taken = utils.runTime(start_time)
				sys.exit('ReMatCh run for pathotyping did not run successfully')
		else:
			clean_pathotyping_folder(args.outdir, original_reference_file, args.debug)
			time_taken = utils.runTime(start_time)
			sys.exit('Something did not run successfully')

	clean_pathotyping_folder(args.outdir, original_reference_file, args.debug)

	print '\n'
	time_taken = utils.runTime(start_time)
	del time_taken


if __name__ == "__main__":
	main()
