import os
import sys
import traceback
import multiprocessing
from Bio import SeqIO
from skDER import util
from collections import defaultdict
import shutil
from operator import itemgetter

def runPyrodigal(genome_listing_file, outdir, combined_proteome_faa, logObject, threads=1, metagenome_mode=False, include_edge_orfs=False):
	"""
	Run Pyrodigal for gene calling on a list of genomes.
	********************************************************************************************************************
	Parameters:
	- genome_listing_file: File containing the list of genomes.
	- outdir: Output directory for Pyrodigal results.
	- combined_proteome_faa: Output file for combined proteome.
	- logObject: Logger object for logging messages.
	- threads: Number of threads to use for Pyrodigal command.
	"""

	# Run pyrodigal for gene calling
	pyrodigal_resdir = outdir + 'pyrodigal_results/'
	util.setupDirectories([pyrodigal_resdir])

	number_of_genomes = None
	pyrodigal_cmds = []
	genome_name_to_path = {}
	try:
		assert(os.path.isfile(genome_listing_file))
		number_of_genomes = 0
		genome_names = set([])
		with open(genome_listing_file) as oaglf:
			for line in oaglf:
				genome_path = line.strip()
				genome_name = '.'.join(line.split('/')[-1].split('.')[:-1])
				genome_name_to_path[genome_name] = genome_path
				genome_names.add(genome_name) 
				result_path = pyrodigal_resdir + genome_name + '.faa'
				pyrodigal_cmd = ['pyrodigal', '-i', genome_path, '-a', result_path]
				if metagenome_mode:
					pyrodigal_cmd += ['-m']
				if not include_edge_orfs:
					pyrodigal_cmd += ['-c']
				pyrodigal_cmds.append(pyrodigal_cmd)
				
				number_of_genomes += 1
		assert(number_of_genomes >= 2 and number_of_genomes == len(genome_names))
	except Exception as e:
		msg = 'Fewer than 2 genomes downloaded / provided or multiple genome files with the same name but from different directories provided. Exiting ...'
		logObject.error(msg)
		sys.stderr.write(msg + '\n')
		sys.exit(1)

	try:
		p = multiprocessing.Pool(threads)
		p.map(util.multiProcess, pyrodigal_cmds)
		p.close()
	except Exception as e:
		msg = 'Issues with parallel running of pyrodigal commands.'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

	cpf_handle = open(combined_proteome_faa, 'w')
	for f in os.listdir(pyrodigal_resdir):
		genome_name = '.faa'.join(f.split('.faa')[:-1])
		with open(pyrodigal_resdir + f) as oprf:
			for rec in SeqIO.parse(oprf, 'fasta'):
				cpf_handle.write('>' + genome_name + '|' + rec.id + '\n' + str(rec.seq) + '\n')
	cpf_handle.close()
				
	# remove indiividual predicted proteome FASTA files
	shutil.rmtree(pyrodigal_resdir)

def runAndProcessCDHIT(combined_proteome_faa, cdhit_result_prefix, cd_hit_params, memory, threads, logObject):
	# run CD-HIT on combined proteomes file
	cdhit_cluster_file = cdhit_result_prefix + '.clstr'
	cdhit_mem = 1000.0*memory
	if cdhit_mem > 0 and cdhit_mem < 1:
		cdhit_mem = 1
	cdhit_mem = int(cdhit_mem)
	cdhit_cmd = ['cd-hit', '-d', '0', '-T', str(threads), '-M', str(cdhit_mem), cd_hit_params, '-i', combined_proteome_faa, '-o', cdhit_result_prefix]
	try:
		util.runCmd(cdhit_cmd, logObject, check_files=[cdhit_cluster_file])
	except:
		msg = 'Issue with running CD-HIT. This is likely because the memory requested was not sufficient, try increaseing memory using the `-m` argument in cidder.'
		sys.stderr.write(msg + '\n') 
		logObject.error(msg)
		sys.exit(1)

	# process CD-HIT clustering resulting results
	genome_protein_clusters = defaultdict(set)
	cluster_genomes = defaultdict(set)
	cluster_id = None
	c_count = 0
	with open(cdhit_cluster_file) as occf:
		for line in occf:
			line = line.strip()
			if line.startswith('>'):
				cluster_id = line[1:]
				c_count += 1
			else:
				genome = line.split()[2][1:].split('|')[0]
				assert(cluster_id != None)
				genome_protein_clusters[genome].add(cluster_id)
				cluster_genomes[cluster_id].add(genome)

	multi_genome_clusters = set([])
	for c in cluster_genomes:
		if len(cluster_genomes[c]) > 1:
			multi_genome_clusters.add(c)
	mgc_count = len(multi_genome_clusters)

	genome_cluster_counts = defaultdict(int)
	for g in genome_protein_clusters:
		gpc = len(genome_protein_clusters[g])
		genome_cluster_counts[genome] = gpc

	return([genome_protein_clusters, genome_cluster_counts, c_count, mgc_count, multi_genome_clusters])


def performRepresentativeGenomeSelection(genome_protein_clusters, genome_cluster_counts, c_count,
		   								 mgc_count, multi_genome_clusters, new_proteins_needed, saturation_cutoff, 
										 multigenome_saturation_cutoff, rep_appending_order_file, genome_name_to_path, 
										 proteome_name_to_path, cidder_drep_dir, logObject, symlink_flag=False):

	raof_handle = open(rep_appending_order_file, 'w')

	rep_genomes = set([])
	rep_genomes_protein_clusters = set([])
	rep_genomes_multigenome_protein_clusters = set([])
	# First, select (one of) the genome(s) with the most distinct protein clusters.
	for i, gc in sorted(genome_cluster_counts.items(), key=itemgetter(1), reverse=True):
		if i == 0: 
			rep_genomes.add(gc[0])
			msg = 'Starting genome: %s - %d distinct protein clusters' % (gc[0], gc[1])
			sys.stdout.write(msg + '\n')
			logObject.info(msg)
			raof_handle.write(gc[0] + '\t0\n')
			file_path = None
			if gc[0] in genome_name_to_path:
				file_path = genome_name_to_path[gc[0]]
			elif gc[0] in proteome_name_to_path:
				file_path = proteome_name_to_path[gc[0]]	

			try:
				assert(os.path.isfile(file_path))
			except:
				msg = 'Genome %s not found in genome_name_to_path or proteome_name_to_path mapping.' % gc[0]
				sys.stderr.write(msg + '\n')
				logObject.error(msg)
				sys.exit(1)
				
			if symlink_flag:
				symlink_file = cidder_drep_dir + file_path.split('/')[-1]
				
				os.symlink(file_path, symlink_file)
			else:
				shutil.copy2(file_path, cidder_drep_dir)
			rep_genomes_protein_clusters = genome_protein_clusters[gc[0]]
			rep_genomes_multigenome_protein_clusters = genome_protein_clusters[gc[0]].intersection(multi_genome_clusters)

	curr_saturation = (len(rep_genomes_protein_clusters)/c_count)*100.0
	curr_multigenome_saturation = (len(rep_genomes_multigenome_protein_clusters)/mgc_count)*100.0

	rep_index = 1
	if curr_saturation >= saturation_cutoff or curr_multigenome_saturation >= multigenome_saturation_cutoff:
		msg = 'Requirements met! Protein cluster saturation of representative genomes is: %0.2f%%\nMulti-genome protein cluster saturation of representative genomes is  %0.2f%%' % (curr_saturation, curr_multigenome_saturation)
		sys.stdout.write(msg + '\n')
		logObject.info(msg)
	else:
		# Start iterative process
		limits_hit = False
		while not limits_hit:
			genome_additions = defaultdict(int)
			for genome in genome_protein_clusters:
				if genome in rep_genomes: continue
				novelty_added = len(genome_protein_clusters[genome].difference(rep_genomes_protein_clusters))
				genome_additions[genome] = novelty_added

			not_enough_new_proteins = False
			new_rep = None
			for ga in sorted(genome_additions.items(), key=itemgetter(1), reverse=True):
				if ga[1] < new_proteins_needed:
					not_enough_new_proteins = True
					break
				else:
					new_rep = ga[0]
					break

			if new_rep != None:
				file_path = None
				if new_rep in genome_name_to_path:
					file_path = genome_name_to_path[new_rep]
				elif new_rep in proteome_name_to_path:
					file_path = proteome_name_to_path[new_rep]	
				if symlink_flag:
					symlink_file = cidder_drep_dir + file_path.split('/')[-1]
					os.symlink(file_path, symlink_file)
				else:
					shutil.copy2(file_path, cidder_drep_dir)
				rep_genomes.add(new_rep)
				raof_handle.write(new_rep + '\t' + str(rep_index) + '\n')
				rep_index += 1
				rep_genomes_protein_clusters = rep_genomes_protein_clusters.union(genome_protein_clusters[new_rep])
				rep_genomes_multigenome_protein_clusters = rep_genomes_multigenome_protein_clusters.union(genome_protein_clusters[new_rep].intersection(multi_genome_clusters))
				msg = 'Adding genome %s' % new_rep
				sys.stdout.write(msg + '\n')
				logObject.info(msg)
			
			curr_saturation = (len(rep_genomes_protein_clusters)/c_count)*100.0
			curr_multigenome_saturation = (len(rep_genomes_multigenome_protein_clusters)/mgc_count)*100.0
			
			if not_enough_new_proteins or curr_saturation >= saturation_cutoff or curr_multigenome_saturation >= multigenome_saturation_cutoff:
				msg = 'Requirements met! Protein cluster saturation of representative genomes is: %0.2f%%\nMulti-genome protein cluster saturation of representative genomes is %0.2f%%' % (curr_saturation, curr_multigenome_saturation)
				sys.stdout.write(msg + '\n')
				logObject.info(msg)
				limits_hit = True

	raof_handle.close()
	return rep_genomes

def appendAdditionalReps(rep_genomes, genome_protein_clusters, required_similarity, rep_appending_order_file,
					     genome_name_to_path, proteome_name_to_path, cidder_drep_dir, logObject, symlink_flag=False):
	"""
	Function to iteratively assess and, if needed, convert non-representative genomes to representatives if >= X%
	of their proteins are not shared with one of the individual representative genomes. Non-representatives 
	are added in order based on the count of unique proteins they have relative to their nearest representative 
	genome.
	"""

	rep_index = len(rep_genomes) + 1

	raof_handle = open(rep_appending_order_file, 'a+')

	total_genomes = len(genome_protein_clusters.keys())

	limits_hit = False
	accounted_nonrep_genomes = set([])
	while not limits_hit:
		genome_additions = defaultdict(int)
		for genome in genome_protein_clusters:
			if genome in rep_genomes: continue
			if genome in accounted_nonrep_genomes: continue

			gpc = genome_protein_clusters[genome]
			difference_to_nearest_rep = 1e10
			for rg in rep_genomes:
				rpc = genome_protein_clusters[rg]
				containment = 100.0*(len(gpc.intersection(rpc))/float(len(gpc)))
				if containment >= required_similarity:
					accounted_nonrep_genomes.add(genome)
				else:
					difference_to_nearest_rep = min([difference_to_nearest_rep, containment])
					genome_additions[genome] = difference_to_nearest_rep

		new_rep = None
		for ga in sorted(genome_additions.items(), key=itemgetter(1), reverse=True):
			new_rep = ga[0]
			break
		
		if new_rep != None:
			file_path = None
			if new_rep in genome_name_to_path:
				file_path = genome_name_to_path[new_rep]
			elif new_rep in proteome_name_to_path:
				file_path = proteome_name_to_path[new_rep]	
			if symlink_flag:
				symlink_file = cidder_drep_dir + file_path.split('/')[-1]
				os.symlink(file_path, symlink_file)
			else:
				shutil.copy2(file_path, cidder_drep_dir)

			msg = 'Adding genome %s' % new_rep
			sys.stdout.write(msg + '\n')
			logObject.info(msg)

			rep_genomes.add(new_rep)
			raof_handle.write(new_rep + '\t' + str(rep_index) + '\n')
			rep_index += 1

			rpc = genome_protein_clusters[new_rep]
			for genome in genome_protein_clusters:
				if genome in rep_genomes: continue
				if genome in accounted_nonrep_genomes: continue
				gpc = genome_protein_clusters[genome]
				containment = 100.0*(len(gpc.intersection(rpc))/float(len(gpc)))
				if containment >= required_similarity:
					accounted_nonrep_genomes.add(genome)

		if len(rep_genomes) + len(accounted_nonrep_genomes) == total_genomes: 
			msg = 'All non-representative genomes within an acceptable distance from representative genomes now ...'
			sys.stdout.write(msg + '\n')
			logObject.info(msg)
			limits_hit = True

	raof_handle.close()
	return rep_genomes

def secondaryClustering(rep_genomes, genome_protein_clusters, cidder_cluster_result_file):
	msg = 'Performing assignment of genomes to their nearest representative genomes based on protein cluster containment.'
	sys.stdout.write(msg + '\n')

	scrf_handle = open(cidder_cluster_result_file, 'w')
	scrf_handle.write('genome\tnearest_representative_genome(s)\tmax_containment_of_genome_protein_clusters\tgenome_protein_cluster_count\trepresentative_genome_protein_cluster_count\n')

	for rg in rep_genomes:
		printlist = [rg, rg, 100.0, len(genome_protein_clusters[rg]), len(genome_protein_clusters[rg])]
		scrf_handle.write('\t'.join([str(x) for x in printlist]) + '\n')

	for genome in genome_protein_clusters:
		if genome in rep_genomes: continue
		fpcs = genome_protein_clusters[genome]
		max_containment = 0.0
		rep_genome_distances = []
		for rep_genome in rep_genomes:
			rpcs = genome_protein_clusters[rep_genome]
			containment = len(fpcs.intersection(rpcs))/float(len(fpcs))
			rep_genome_distances.append([rep_genome, containment])
			if max_containment < containment:
				max_containment = containment

		superset_refs = set([])
		for rgr in sorted(rep_genome_distances, key=itemgetter(1), reverse=True):
			if rgr[1] == max_containment:
				superset_refs.add(rgr[0])
			else:
				break

		printlist = [genome, ', '.join(sorted(superset_refs)), 
						str(max_containment*100.0), str(len(fpcs)), str(len(rpcs))]
		scrf_handle.write('\t'.join(printlist) + '\n')

	scrf_handle.close()

def secondaryClusteringSkani(all_genomes_listing_file, rep_genomes, genome_protein_clusters, skani_result_file, cidder_cluster_result_file, logObject, threads=1):
	msg = 'Performing assignment of genomes to their nearest representative genomes based on skani ANI estimates.'
	sys.stdout.write(msg + '\n')
	
	# run skani triangle	
	skani_triangle_cmd = ['skani', 'triangle', '-l', all_genomes_listing_file, '--min-af', '80.0', '-E', '-t', str(threads), '-o', skani_result_file]
	util.runCmd(skani_triangle_cmd, logObject, check_files=[skani_result_file])

	scrf_handle = open(cidder_cluster_result_file, 'w')
	scrf_handle.write('genome\tnearest_representative_genome\taverage_nucleotide_identity\talignment_fraction\n')

	for rg in rep_genomes:
		scrf_handle.write(rg + '\t' + rg + '\t100.0\t100.0\n')

	best_rep_match_at_default_af = defaultdict(lambda: [set(["NA"]), 0.0, 0.0])
	with open(skani_result_file) as osrf:
		for i, line in enumerate(osrf):
			if i == 0: continue
			line = line.strip()
			ref, que, ani, raf, qaf, _, _ = line.split('\t')

			ref = '.'.join(ref.split('/')[-1].split('.')[:-1])
			que = '.'.join(que.split('/')[-1].split('.')[:-1])

			ani = float(ani)
			raf = float(raf)
			qaf = float(qaf)
			if que in rep_genomes and not ref in rep_genomes:
				if ani > best_rep_match_at_default_af[ref][1]:
					best_rep_match_at_default_af[ref] = [set([que]), ani, raf]
				elif ani == best_rep_match_at_default_af[ref][1]:
					if raf > best_rep_match_at_default_af[ref][2]:
						best_rep_match_at_default_af[ref] = [set([que]), ani, raf]
					elif raf == best_rep_match_at_default_af[ref][2]:
						best_rep_match_at_default_af[ref][0].add(que)

			if ref in rep_genomes and not que in rep_genomes:
				if ani > best_rep_match_at_default_af[que][1]:
					best_rep_match_at_default_af[que] = [set([ref]), ani, qaf]
				elif ani == best_rep_match_at_default_af[que][1]:
					if qaf > best_rep_match_at_default_af[que][2]:
						best_rep_match_at_default_af[que] = [set([ref]), ani, qaf]
					elif qaf == best_rep_match_at_default_af[que][2]:
						best_rep_match_at_default_af[que][0].add(ref)
										
	for gen in genome_protein_clusters:
		if gen in rep_genomes: continue
		scrf_handle.write('\t'.join([gen, ', '.join(sorted(best_rep_match_at_default_af[gen][0])), 
							str(best_rep_match_at_default_af[gen][1]), 
							str(best_rep_match_at_default_af[gen][2])]) + '\n')
			
	scrf_handle.close()