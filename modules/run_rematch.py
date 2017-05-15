import functools
import os
import sys

import utils


# {'noMatter': '/home/ubuntu/NGStools/patho_typing/mpmachado_stuff.out_test/rematch/sample.noMatter.fasta', 'correct': '/home/ubuntu/NGStools/patho_typing/mpmachado_stuff.out_test/rematch/sample.correct.fasta', 'alignment': '/home/ubuntu/NGStools/patho_typing/mpmachado_stuff.out_test/rematch/sample.alignment.fasta'}


def remove_alignment(alignment_file):
	directory = os.path.dirname(alignment_file)
	files = [f for f in os.listdir(directory) if not f.startswith('.') and os.path.isfile(os.path.join(directory, f))]
	for file_found in files:
		if file_found.startswith(os.path.splitext(os.path.basename(alignment_file))[0]):
			file_found = os.path.join(directory, file_found)
			os.remove(file_found)


def remove_reference_stuff(outdir, reference_file):
	files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
	for file_found in files:
		if file_found.startswith(os.path.splitext(os.path.basename(reference_file))[0]):
			file_found = os.path.join(outdir, file_found)
			os.remove(file_found)


def clean_rematch_folder(consensus_files, bam_file, reference_file, outdir, doNotRemoveConsensus, debug_mode_true):
	if not debug_mode_true:
		if not doNotRemoveConsensus:
			for consensus_type, file_path in consensus_files.items():
				if os.path.isfile(file_path):
					os.remove(file_path)
		if bam_file is not None:
			remove_alignment(bam_file)
		remove_reference_stuff(outdir, reference_file)


def sequence_data(sample, reference_file, bam_file, outdir, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, debug_mode_true, rematch):
	sequence_data_outdir = os.path.join(outdir, 'sequence_data', '')
	utils.removeDirectory(sequence_data_outdir)
	os.mkdir(sequence_data_outdir)

	sequences, headers = utils.get_sequence_information(reference_file, length_extra_seq)

	threads_2_use = rematch.determine_threads_2_use(len(sequences), threads)

	import multiprocessing

	pool = multiprocessing.Pool(processes=threads)
	for sequence_counter in sequences:
		sequence_dir = os.path.join(sequence_data_outdir, str(sequence_counter), '')
		utils.removeDirectory(sequence_dir)
		os.makedirs(sequence_dir)
		pool.apply_async(rematch.analyse_sequence_data, args=(bam_file, sequences[sequence_counter], sequence_dir, sequence_counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, threads_2_use,))
	pool.close()
	pool.join()

	run_successfully, sample_data, consensus_files, consensus_sequences = rematch.gather_data_together(sample, sequence_data_outdir, sequences, outdir.rsplit('/', 2)[0], debug_mode_true, length_extra_seq)

	return run_successfully, sample_data, consensus_files, consensus_sequences


def write_report(outdir, sample_data, minimum_gene_coverage, minimum_gene_identity):
	print 'Writing report file'
	number_absent_genes = 0
	number_genes_multiple_alleles = 0
	mean_sample_coverage = 0
	with open(os.path.join(outdir, 'rematchModule_report.txt'), 'wt') as writer:
		writer.write('\t'.join(['#gene', 'percentage_gene_coverage', 'gene_mean_read_coverage', 'percentage_gene_low_coverage', 'number_positions_multiple_alleles', 'percentage_gene_identity']) + '\n')
		for i in range(1, len(sample_data) + 1):
			writer.write('\t'.join([sample_data[i]['header'], str(round(sample_data[i]['gene_coverage'], 2)), str(round(sample_data[i]['gene_mean_read_coverage'], 2)), str(round(sample_data[i]['gene_low_coverage'], 2)), str(sample_data[i]['gene_number_positions_multiple_alleles']), str(round(sample_data[i]['gene_identity'], 2))]) + '\n')

			if sample_data[i]['gene_coverage'] < minimum_gene_coverage or sample_data[i]['gene_identity'] < minimum_gene_identity:
				number_absent_genes += 1
			else:
				mean_sample_coverage += sample_data[i]['gene_mean_read_coverage']
				if sample_data[i]['gene_number_positions_multiple_alleles'] > 0:
					number_genes_multiple_alleles += 1

		if len(sample_data) - number_absent_genes > 0:
			mean_sample_coverage = float(mean_sample_coverage) / float(len(sample_data) - number_absent_genes)
		else:
			mean_sample_coverage = 0

		writer.write('\n'.join(['#general', '>number_absent_genes', str(number_absent_genes), '>number_genes_multiple_alleles', str(number_genes_multiple_alleles), '>mean_sample_coverage', str(round(mean_sample_coverage, 2))]) + '\n')

		print '\n'.join([str('number_absent_genes: ' + str(number_absent_genes)), str('number_genes_multiple_alleles: ' + str(number_genes_multiple_alleles)), str('mean_sample_coverage: ' + str(round(mean_sample_coverage, 2)))]) + '\n'

	return number_absent_genes, number_genes_multiple_alleles, mean_sample_coverage


module_timer = functools.partial(utils.timer, name='Module ReMatCh')


@module_timer
def run_rematch(rematch, outdir, reference_file, bam_file, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, minimum_gene_coverage, minimum_gene_identity, debug_mode_true, doNotRemoveConsensus):
	module_dir = os.path.join(outdir, 'rematch', '')
	utils.removeDirectory(module_dir)
	os.makedirs(module_dir)

	sys.path.append(os.path.join(os.path.dirname(rematch), 'modules', ''))
	import rematch_module as rematch

	print 'Analysing alignment data'
	run_successfully, sample_data, consensus_files, consensus_sequences = sequence_data('sample', reference_file, bam_file, module_dir, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, debug_mode_true, rematch)

	if run_successfully:
		number_absent_genes, number_genes_multiple_alleles, mean_sample_coverage = write_report(outdir, sample_data, minimum_gene_coverage, minimum_gene_identity)

	if not debug_mode_true:
		utils.removeDirectory(module_dir)

	clean_rematch_folder(consensus_files, bam_file, reference_file, outdir, doNotRemoveConsensus, debug_mode_true)

	return run_successfully, {'number_absent_genes': number_absent_genes if 'number_absent_genes' in locals() else None, 'number_genes_multiple_alleles': number_genes_multiple_alleles if 'number_genes_multiple_alleles' in locals() else None, 'mean_sample_coverage': round(mean_sample_coverage, 2) if 'mean_sample_coverage' in locals() else None}, sample_data if 'sample_data' in locals() else None
