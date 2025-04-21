import os
import sys
import gzip
import traceback
from skDER import util
import subprocess
from collections import defaultdict

def runSkani(genome_listing_file, skani_result_file, skani_triangle_parameters, 
			 aligned_fraction_cutoff, selection_mode, test_cutoffs_flag, logObject, threads=1):	
	"""
	Run skani triangle command to find the best ANI cutoffs for the genomes.
		:param genome_listing_file: File containing the list of genomes.
		:param skani_result_file: File to save the results.
		:param skani_triangle_parameters: Parameters for skani triangle command.
		:param aligned_fraction_cutoff: Cutoff for aligned fraction.
		:param selection_mode: Mode for selecting genomes.
		:param test_cutoffs_flag: Flag to test cutoffs.
		:param logObject: Logger object for logging messages.
		:param threads: Number of threads to use for skani command.
	"""
	try:
		skani_triangle_cmd = ['skani', 'triangle', '-l', genome_listing_file, 
							'--min-af', str(aligned_fraction_cutoff), '-E', skani_triangle_parameters, 
							'-t', str(threads), '-o', skani_result_file]
		if test_cutoffs_flag:
			min_af_cutoff = 10.0
			if selection_mode == 'dynamic':
				min_af_cutoff = max([min_af_cutoff - 20.0, 0.0])
			skani_triangle_cmd = ['skani', 'triangle', '-l', genome_listing_file,
								'--min-af', str(min_af_cutoff), '-E', skani_triangle_parameters, 
								'-t', str(threads), '-o', skani_result_file]
		util.runCmd(skani_triangle_cmd, logObject, check_files=[skani_result_file])
	except Exception as e:
		raise RuntimeError('Error running skani triangle command: %s' % str(e))	

def dynamicDerep(skani_result_file, concat_n50_result_file, skder_result_file, outdir,
				 ani_cutoff, af_cutoff, max_af_dist, logObject, mge_proc_to_unproc_mapping=None):
	# perform representative selection using dynamic method (default)
	skder_core_prog = 'skDERcore'

	try:
		result = subprocess.run([skder_core_prog], capture_output=True, text=True)
		assert(result.stdout.startswith("Usage:\nskDERcore"))
	except:
		sys.stderr.write('skDERcore command not found or not executable.\n')
		sys.exit(1)

	if mge_proc_to_unproc_mapping is None:
		skder_core_cmd = [skder_core_prog, skani_result_file, concat_n50_result_file, str(ani_cutoff), 
							str(af_cutoff), str(max_af_dist), '>', skder_result_file]
		util.runCmd(skder_core_cmd, logObject, check_files=[skder_result_file])
	else:
		skder_result_file_mge_paths = outdir + 'tmp_skDER_Results.txt'
		skder_core_cmd = [skder_core_prog, skani_result_file, concat_n50_result_file, str(ani_cutoff), 
							str(af_cutoff), str(max_af_dist), '>', skder_result_file_mge_paths]
		util.runCmd(skder_core_cmd, logObject, check_files=[skder_result_file_mge_paths])
		
		srf_handle = open(skder_result_file, 'w')
		with open(skder_result_file_mge_paths) as osrfmp:
			for line in osrfmp:
				line = line.strip()
				unproc_genome = mge_proc_to_unproc_mapping[line]
				srf_handle.write(unproc_genome + "\n")
		srf_handle.close()

def greedyDerep(skani_result_file, concat_n50_result_file, skder_result_file, outdir,
				ani_cutoff, af_cutoff, logObject, mge_proc_to_unproc_mapping=None, threads=1):
	# perform representative selection using greedy method
	skder_sum_prog = 'skDERsum'
	genome_summary_file = outdir + 'Genome_Information_for_Greedy_Clustering.txt'
	skder_sum_cmd = [skder_sum_prog, skani_result_file, concat_n50_result_file, str(ani_cutoff), 
					str(af_cutoff), '>', genome_summary_file]
	util.runCmd(skder_sum_cmd, logObject, check_files=[genome_summary_file])

	sorted_genome_summary_file = outdir + 'Genome_Information_for_Greedy_Clustering.sorted.txt'
	sort_cmd = ['sort', '-k', '2', '--parallel=' + str(threads), '-gr', genome_summary_file, '>', sorted_genome_summary_file]
	util.runCmd(sort_cmd, logObject, check_files=[sorted_genome_summary_file])

	# greedy clustering
	skder_result_handle = open(skder_result_file, 'w')
	already_accounted = set([])
	with open(sorted_genome_summary_file) as osgsf:
		for line in osgsf:
			line = line.strip('\n')
			ls = line.split('\t')
			curr_g = ls[0]
			if curr_g in already_accounted: continue
			ls = line.split('\t')
			for g in ls[2].split('; '):
				already_accounted.add(g)
			if mge_proc_to_unproc_mapping != None:
				skder_result_handle.write(mge_proc_to_unproc_mapping[curr_g] + '\n')
			else:
				skder_result_handle.write(curr_g + '\n')
	skder_result_handle.close()


def determineClusters(skder_result_file, skani_result_file, mge_unproc_to_proc_mapping,
					  mge_proc_to_unproc_mapping, skder_cluster_result_file, specified_af_cutoff,
					  specified_ani_cutoff):
	
	scrf_handle = open(skder_cluster_result_file, 'w')
	scrf_handle.write('genome\tnearest_representative_genome\taverage_nucleotide_identity\talignment_fraction\tmatch_category\n')
	rep_genomes = set([])
	with open(skder_result_file) as osrf:
		for line in osrf:
			line = line.strip()
			if mge_unproc_to_proc_mapping != None:
				rep_genomes.add(mge_unproc_to_proc_mapping[line])
			else:
				rep_genomes.add(line)
			scrf_handle.write(line + '\t' + line + '\t100.0\t100.0\trepresentative_to_self\n')

	best_rep_match_at_default_af = defaultdict(lambda: [set(["NA"]), 0.0, 0.0])
	best_rep_match_at_loose_af = defaultdict(lambda: [set(["NA"]), 0.0, 0.0])
	with open(skani_result_file) as osrf:
		for i, line in enumerate(osrf):
			if i == 0: continue
			line = line.strip()
			ref, que, ani, raf, qaf, _, _ = line.split('\t')
			ani = float(ani)
			raf = float(raf)
			qaf = float(qaf)
			if que in rep_genomes and not ref in rep_genomes:
				if raf >= specified_af_cutoff:
					if ani > best_rep_match_at_default_af[ref][1]:
						best_rep_match_at_default_af[ref] = [set([que]), ani, raf]
					elif ani == best_rep_match_at_default_af[ref][1]:
						if raf > best_rep_match_at_default_af[ref][2]:
							best_rep_match_at_default_af[ref] = [set([que]), ani, raf]
						elif raf == best_rep_match_at_default_af[ref][2]:
							best_rep_match_at_default_af[ref][0].add(que)
				else:
					if ani > best_rep_match_at_loose_af[ref][1]:
						best_rep_match_at_loose_af[ref] = [set([que]), ani, raf]
					elif ani == best_rep_match_at_loose_af[ref][1]:
						if raf > best_rep_match_at_loose_af[ref][2]:
							best_rep_match_at_loose_af[ref] = [set([que]), ani, raf]
						elif raf == best_rep_match_at_loose_af[ref][2]:
							best_rep_match_at_loose_af[ref][0].add(que)

			if ref in rep_genomes and not que in rep_genomes:
				if qaf >= specified_af_cutoff:
					if ani > best_rep_match_at_default_af[que][1]:
						best_rep_match_at_default_af[que] = [set([ref]), ani, qaf]
					elif ani == best_rep_match_at_default_af[que][1]:
						if qaf > best_rep_match_at_default_af[que][2]:
							best_rep_match_at_default_af[que] = [set([ref]), ani, qaf]
						elif qaf == best_rep_match_at_default_af[que][2]:
							best_rep_match_at_default_af[que][0].add(ref)
				else:
					if ani > best_rep_match_at_loose_af[que][1]:
						best_rep_match_at_loose_af[que] = [set([ref]), ani, qaf]
					elif ani == best_rep_match_at_loose_af[que][1]:
						if qaf > best_rep_match_at_loose_af[que][2]:
							best_rep_match_at_loose_af[que] = [set([ref]), ani, qaf]
						elif qaf == best_rep_match_at_loose_af[que][2]:
							best_rep_match_at_loose_af[que][0].add(ref)				
										
	if mge_proc_to_unproc_mapping != None:
		strict_nearest_ref_found = set([])
		for gen in best_rep_match_at_default_af:
			if best_rep_match_at_default_af[gen][0] != 'NA':
				strict_nearest_ref_found.add(gen)
				if best_rep_match_at_default_af[gen][1] >= specified_ani_cutoff:
					scrf_handle.write('\t'.join([mge_proc_to_unproc_mapping[gen], 
									', '.join([mge_proc_to_unproc_mapping[x] for x in best_rep_match_at_default_af[gen][0]]), 
									str(best_rep_match_at_default_af[gen][1]), 
									str(best_rep_match_at_default_af[gen][2]), 
									'within_cutoffs_requested']) + '\n')
				else:
					scrf_handle.write('\t'.join([mge_proc_to_unproc_mapping[gen], 
							', '.join([mge_proc_to_unproc_mapping[x] for x in best_rep_match_at_default_af[gen][0]]), 
							str(best_rep_match_at_default_af[gen][1]), 
							str(best_rep_match_at_default_af[gen][2]), 
							'outside_cutoffs_requested']) + '\n')
		for gen in best_rep_match_at_loose_af:
			if gen in strict_nearest_ref_found: continue
			scrf_handle.write('\t'.join([mge_proc_to_unproc_mapping[gen], 
							', '.join([mge_proc_to_unproc_mapping[x] for x in best_rep_match_at_loose_af[gen][0]]), 
							str(best_rep_match_at_loose_af[gen][1]), 
							str(best_rep_match_at_loose_af[gen][2]), 
							'outside_cutoffs_requested']) + '\n')

	else:
		strict_nearest_ref_found = set([])
		for gen in best_rep_match_at_default_af:
			if best_rep_match_at_default_af[gen][0] != 'NA':
				strict_nearest_ref_found.add(gen)
				if best_rep_match_at_default_af[gen][1] >= specified_ani_cutoff:
					scrf_handle.write('\t'.join([gen, ', '.join(best_rep_match_at_default_af[gen][0]), 
									str(best_rep_match_at_default_af[gen][1]), 
									str(best_rep_match_at_default_af[gen][2]), 
									'within_cutoffs_requested']) + '\n')
				else:
					scrf_handle.write('\t'.join([gen, ', '.join(best_rep_match_at_default_af[gen][0]), 
					str(best_rep_match_at_default_af[gen][1]), 
					str(best_rep_match_at_default_af[gen][2]), 
					'outside_cutoffs_requested']) + '\n')
		for gen in best_rep_match_at_loose_af:
			if gen in strict_nearest_ref_found: continue
			scrf_handle.write('\t'.join([gen, ', '.join(best_rep_match_at_loose_af[gen][0]), 
							str(best_rep_match_at_loose_af[gen][1]), 
							str(best_rep_match_at_loose_af[gen][2]), 
							'outside_cutoffs_requested']) + '\n')

	scrf_handle.close()
