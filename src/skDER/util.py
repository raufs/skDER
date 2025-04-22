import os
import sys
import logging
import traceback
import subprocess
from numpy import where, cumsum
from Bio import SeqIO
import gzip
from collections import defaultdict
import asyncio
import aiohttp 
import aiofile
import math
import multiprocessing
import shutil
import tqdm
import resource

ACCEPTED_FASTA_SUFFICES = set(['fasta', 'fas', 'fna', 'fa'])
ACCEPTED_PROTEIN_FASTA_SUFFICES = set(['fasta', 'faa', 'fa'])

def memory_limit(mem):
	"""
	Description:
	Experimental function to limit memory.
	********************************************************************************************************************
	Parameters:
	- mem: The memory limit in GB.
	********************************************************************************************************************
	"""
	max_virtual_memory = mem*1000000000
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, hard))
	print(resource.getrlimit(resource.RLIMIT_AS))

def _download_files(urls, resdir):
	"""
	Download files from the given URLs and save them to the specified directory.
	Note, this function was taken from: 
	https://gist.github.com/darwing1210/c9ff8e3af8ba832e38e6e6e347d9047a
	********************************************************************
	Parameters:
	- urls: List of URLs to download.	
	- resdir: Directory to save the downloaded files.
	********************************************************************
	"""
	os.makedirs(resdir, exist_ok=True)
	sema = asyncio.BoundedSemaphore(5)

	async def fetch_file(session, url):
		fname = url.split("/")[-1]
		async with sema:
			async with session.get(url) as resp:
				assert resp.status == 200
				data = await resp.read()

		async with aiofile.async_open(
			os.path.join(resdir, fname), "wb"
		) as outfile:
			await outfile.write(data)

	async def main():
		async with aiohttp.ClientSession() as session:
			tasks = [fetch_file(session, url) for url in urls]
			await asyncio.gather(*tasks)

	loop = asyncio.get_event_loop()
	loop.run_until_complete(main())
	loop.close()

def downloadGTDBGenomes(taxa_name, gtdb_release, outdir, genome_listing_file, logObject, sanity_check=False, automated_download=False):
	"""
	Download GTDB genomes from NCBI Genbank using ncbi-genome-download.
	**********************************************************
	Parameters:
		- logObject: Logger object for logging messages.
		- taxa_name: Taxa name to search for in GTDB.
		:param gtdb_release: GTDB release version.
		:param genomes: List of genomes to download.
		:param outdir: Output directory for downloaded genomes.
		:param ngd_url: URL for ncbi-genome-download.
		:param sanity_check: Boolean flag for sanity check.	    
	"""
	
	""" 
	Download GTDB listing file from lsaBGC git repo, parse GTDB information 
	file, get list of Genbank accessions, and perform dry-run with ncbi-genome-download 
	if requested.
	"""
	msg = "Using wget to download GTDB listing."
	sys.stdout.write(msg + '\n')
	logObject.info(msg)
	wget_cmd = ['wget', '-q', 'https://github.com/raufs/gtdb_gca_to_taxa_mappings/raw/main/GTDB_' + gtdb_release + '_Information_with_Genome_URLs.txt.gz', '-P', outdir]
	gtdb_listing_file = outdir + "GTDB_" + gtdb_release + "_Information_with_Genome_URLs.txt.gz"
	runCmd(wget_cmd, logObject, check_files=[gtdb_listing_file])

	genbank_accession_listing_file = outdir + 'GTDB_Genomes_to_Download.txt'
	msg = "Beginning by assessing which genomic assemblies are available for the taxa %s in GTDB %s" % (taxa_name, gtdb_release)
	sys.stdout.write(msg + '\n')
	logObject.info(msg) 
	
	select_genome_listing_handle = open(genbank_accession_listing_file, 'w')
	genome_url_paths = []
	url_file_to_polished_name = {}
	with gzip.open(gtdb_listing_file, 'rt') as ogtdb:
		for i, line in enumerate(ogtdb):
			line = line.strip('\n')
			if i == 0:
				select_genome_listing_handle.write(line + '\n')
				continue
			ls = line.split('\t')
			gca, gtdb_genus, gtdb_species, genome_url, version_match = line.split('\t')
			if gca == 'none': continue
			if genome_url == "NA": continue
			if len(taxa_name.split()) == 1:
				if gtdb_genus == taxa_name:
					select_genome_listing_handle.write(line + '\n')
					genome_url_paths.append(genome_url)
					filename = genome_url.split('/')[-1]
					url_file_to_polished_name[filename] = '_'.join(gtdb_species.split()) + '_' + gca + '.fasta.gz'
			elif len(taxa_name.split()) == 2:
				if gtdb_species == taxa_name:
					select_genome_listing_handle.write(line + '\n')
					genome_url_paths.append(genome_url)
					filename = genome_url.split('/')[-1]
					url_file_to_polished_name[filename] = '_'.join(gtdb_species.split()) + '_' + gca + '.fasta.gz'

	select_genome_listing_handle.close()

	genome_count = len(genome_url_paths)
	if genome_count == 0:
		msg = "Error: no genomes found to belong the genus or species specified in GTDB."
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		sys.exit(1)
	else:
		if not automated_download:
			try:
				response = input("Will be downloading %d genomes for the taxon %s. Note, each\ngenome in compressed FASTA formatting is approximately 1-2 MB.\nIf downloading thousands of genomes, this can lead to significant disk space being\nused. Do you wish to continue? (yes/no): " % (genome_count, taxa_name))
				if response.lower() != 'yes':
					os.system('User does NOT want to download genomes for taxa from NCBI. Exiting ...')
					sys.exit(1)
			except:
				msg = 'Error: user did not respond to download genomes for taxa from NCBI. Exiting ...'
				sys.stderr.write(msg + '\n')
				logObject.info(msg)
				sys.exit(1)

		genomes_directory = outdir + 'gtdb_ncbi_genomes/'
		try:
			_download_files(genome_url_paths, genomes_directory)
		except Exception as e:
			msg = 'Error downloading genomes from NCBI. Exiting ...'
			sys.stderr.write(msg + '\n')
			logObject.info(msg)
			sys.exit(1)

		gf_listing_handle = open(genome_listing_file, 'a+')
		final_genome_count = 0
		for f in os.listdir(genomes_directory):
			genome_file = genomes_directory + f
			if not os.path.isfile(genome_file): continue
			if sanity_check:
				try:
					assert (is_fasta(genome_file))
					final_genome_count += 1
				except AssertionError as e:
					msg = 'Warning: genome %s is not a valid FASTA file. Removing ...' % f
					try:
						os.remove(genome_file)
					except:
						pass
					sys.stderr.write(msg + '\n')
					logObject.info(msg)
			else:
				final_genome_count += 1
			gca = '_'.join(f.split('_')[:2])
			polished_filename = url_file_to_polished_name[f]
			renamed_gfile = genomes_directory + polished_filename
			os.rename(genome_file, renamed_gfile)
			gf_listing_handle.write(renamed_gfile + '\n')
		gf_listing_handle.close()
		
		msg = 'Was able to download %d of %d genomes belonging to taxa "%s" in GTDB %s.' % (final_genome_count, genome_count, taxa_name, gtdb_release)
		sys.stdout.write(msg + '\n')
		logObject.info(msg)

def processInputProteomes(proteomes, combined_proteome_faa, logObject, sanity_check=False):
	combined_proteome_handle = open(combined_proteome_faa, 'a+')
	proteome_name_to_path = {}
	for p in proteomes:
		p_path = os.path.abspath(p)
		if os.path.isdir(p_path):
			for f in os.listdir(p_path):
				proteome_file = p_path + '/' + f
				proteome_name = '.'.join(proteome_file.split('/')[-1].split('.')[:-1])
				suffix = proteome_file.split('.')[-1].lower()
				if proteome_file.endswith('.gz'):
					proteome_name = '.'.join(proteome_file.split('/')[-1][:-3].split('.')[:-1])
					suffix = proteome_file[:-3].split('.')[-1].lower()
				
				try:
					assert (suffix in ACCEPTED_PROTEIN_FASTA_SUFFICES)
				except:
					msg = 'Warning: proteome %s does not have a valid FASTA suffix. Skipping ...' % proteome_file
					sys.stderr.write(msg + '\n')
					logObject.warning(msg)
					continue 
				
				if sanity_check:
					try:
						assert (is_fasta(proteome_file))
						with open(proteome_file) as of:
							for rec in SeqIO.parse(of, 'fasta'):
								combined_proteome_handle.write('>' + proteome_name + '|' + rec.id + '\n' + str(rec.seq) + '\n')
						proteome_name_to_path[proteome_name] = proteome_file
					except Exception as e:
						msg = 'Warning: proteome %s is not a valid FASTA file. Skipping ...' % proteome_file
						sys.stderr.write(msg + '\n')
						logObject.warning(msg)
						continue
				else:
					try:
						with open(proteome_file) as of:
							for rec in SeqIO.parse(of, 'fasta'):
								combined_proteome_handle.write('>' + proteome_name + '|' + rec.id + '\n' + str(rec.seq) + '\n')			
						proteome_name_to_path[proteome_name] = proteome_file
					except Exception as e:
						msg = 'Warning: unable to parse proteome %s. Skipping ...' % proteome_file
						sys.stderr.write(msg + '\n')
						logObject.warning(msg)
						continue
		elif os.path.isfile(p_path):
			proteome_file = p_path
			proteome_name = '.'.join(proteome_file.split('/')[-1].split('.')[:-1])
			suffix = proteome_file.split('/')[-1].split('.')[-1].lower()
			if proteome_file.endswith('.gz'):
				proteome_name = '.'.join(proteome_file.split('/')[-1][:-3].split('.')[:-1])
				suffix = proteome_file.split('/')[-1][:-3].split('.')[-1].lower()
			
			try:
				assert (suffix in ACCEPTED_PROTEIN_FASTA_SUFFICES)
			except:
				msg = 'Warning: proteome %s does not have a valid FASTA suffix. Skipping ...' % proteome_file
				sys.stderr.write(msg + '\n')
				logObject.warning(msg)
				continue 
			
			if sanity_check:
				try:
					assert (is_fasta(proteome_file))
					with open(proteome_file) as of:
						for rec in SeqIO.parse(of, 'fasta'):
							combined_proteome_handle.write('>' + proteome_name + '|' + rec.id + '\n' + str(rec.seq) + '\n')
					proteome_name_to_path[proteome_name] = proteome_file
				except Exception as e:
					msg = 'Warning: proteome %s is not a valid FASTA file. Skipping ...' % proteome_file
					sys.stderr.write(msg + '\n')
					logObject.warning(msg)
					continue
			else:
				try:
					with open(proteome_file) as of:
						for rec in SeqIO.parse(of, 'fasta'):
							combined_proteome_handle.write('>' + proteome_name + '|' + rec.id + '\n' + str(rec.seq) + '\n')			
					proteome_name_to_path[proteome_name] = proteome_file
				except Exception as e:
					msg = 'Warning: unable to parse proteome %s. Skipping ...' % proteome_file
					sys.stderr.write(msg + '\n')
					logObject.warning(msg)
					continue
		else:
			msg = 'Warning: unable to validate potential proteome %s exists as a file/directory. Skipping ...' % p
			sys.stderr.write(msg + '\n')
			logObject.info(msg)
	combined_proteome_handle.close()
	return proteome_name_to_path

def processInputGenomes(genomes, genome_listing_file, logObject, sanity_check=False):
	gf_listing_handle = open(genome_listing_file, 'a+')
	for g in genomes:
		g_path = os.path.abspath(g)
		if os.path.isdir(g_path):
			for f in os.listdir(g_path):
				genome_file = g_path + '/' + f
				suffix = genome_file.split('/')[-1].split('.')[-1].lower()
				if genome_file.endswith('.gz'):
					suffix = genome_file.split('/')[-1][:-3].split('.')[-1].lower()
				
				try:
					assert (suffix in ACCEPTED_FASTA_SUFFICES)
				except:
					msg = 'Warning: genome %s does not have a valid FASTA suffix. Skipping ...' % genome_file
					sys.stderr.write(msg + '\n')
					logObject.warning(msg)
					continue

				if sanity_check:
					try:
						assert (is_fasta(genome_file))
						gf_listing_handle.write(genome_file + '\n')
					except:
						msg = 'Warning: genome %s is not a valid FASTA file. Skipping ...' % genome_file
						sys.stderr.write(msg + '\n')
						logObject.warning(msg)
						continue
				else:
					gf_listing_handle.write(genome_file + '\n')
		elif os.path.isfile(g_path):
			genome_file = os.path.abspath(g_path)
			suffix = genome_file.split('/')[-1].split('.')[-1].lower()
			if genome_file.endswith('.gz'):
				suffix = genome_file.split('/')[-1][:-3].split('.')[-1].lower()
			
			try:
				assert (suffix in ACCEPTED_FASTA_SUFFICES)
			except:
				msg = 'Warning: genome %s does not have a valid FASTA suffix. Skipping ...' % genome_file
				sys.stderr.write(msg + '\n')
				logObject.warning(msg)
				continue

			if sanity_check:
				try:
					assert (is_fasta(genome_file))
					gf_listing_handle.write(genome_file + '\n')
				except:
					msg = 'Warning: genome %s is not a valid FASTA file. Skipping ...' % genome_file
					sys.stderr.write(msg + '\n')
					logObject.warning(msg)
					continue
			else:
				gf_listing_handle.write(genome_file + '\n')
		else:
			msg = 'Warning: unable to validate genome %s exists as a file. Skipping ...' % g
			sys.stderr.write(msg + '\n')
			logObject.info(msg)
	gf_listing_handle.close()

def copyRepresentativeGenomesToDirectory(skder_result_file, outdir, logObject, symlink_flag=False):
	# copy over genomes which are non-redundant to a separate directory
	skder_drep_dir = outdir + 'Dereplicated_Representative_Genomes/'	
	setupDirectories([skder_drep_dir])

	with open(skder_result_file) as osrf:
		for line in osrf:
			genome_path = line.strip()
			try:
				if symlink_flag:
					symlink_file = skder_drep_dir + genome_path.split('/')[-1]
					os.symlink(genome_path, symlink_file)
				else:
					shutil.copy2(genome_path, skder_drep_dir)				
			except:
				sys.stderr.write('Warning: issues copying over representative genome %s to final dereplicated sub-directory.\n' % genome_path)
				logObject.warning('Issues copying over representative genome %s to final dereplicated sub-directory.' % genome_path)

def determineN50(genome_listing_file, outdir, logObject, threads=1):
	genomes = []
	with open(genome_listing_file) as oaglf:
		for line in oaglf:
			genome_path = line.strip()
			genomes.append(genome_path)

	n50_dir = outdir + 'N50/'
	setupDirectories([n50_dir])

	genome_count = len(genomes)
	chunk_size = math.ceil(genome_count/threads)
	genome_chunks = divide_chunks(genomes, chunk_size)

	n50_inputs = []
	for i, gc in enumerate(genome_chunks):
		n50_chunk_file = n50_dir + 'chunk_' + str(i) + '.txt'
		n50_inputs.append([gc, n50_chunk_file])

	try:
		msg = "Calculating assembly N50 for %d genomic assemblies" % genome_count
		sys.stdout.write(msg + '\n')
		logObject.info(msg)
		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(compute_n50, n50_inputs), total=len(n50_inputs)):
			pass
		p.map(compute_n50, n50_inputs)
		p.close()
	except Exception as e:
		msg = "An error occurred during multiprocessing: %s" % str(e)
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
	
	# store results in a dictionary
	n50s = {}
	for f in os.listdir(n50_dir):
		with open(n50_dir + f) as on50f:
			for line in on50f:
				line = line.strip()
				ls = line.split('\t')
				genome, n50_value = ls[0], int(ls[1]) 
				n50s[genome] = n50_value

	# remove the directory with N50 stats after concatenating info into a single file.
	shutil.rmtree(n50_dir)
	
	return(n50s)

def writeN50sToFile(n50s, all_genomes_listing_file, concat_n50_result_file, logObject):
	"""
	Write N50 values to a file.
	********************************************************************************************************************
	Parameters:
	- n50s: Dictionary of N50 values for each genome.
	- all_genomes_listing_file: File containing the list of genomes.
	- concat_n50_result_file: File to save the concatenated N50 results.
	********************************************************************************************************************
	"""
	
	try:
		glmf_handle = open(all_genomes_listing_file, 'r')
	except:
		msg = 'Unable to open genome listing file %s.' % all_genomes_listing_file
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		sys.exit(1)

	n50_handle = open(concat_n50_result_file, 'w')
	for line in glmf_handle:
		line = line.strip()
		if line in n50s:
			n50_handle.write(line + '\t' + str(n50s[line]) + '\n')
	glmf_handle.close()
	n50_handle.close()

def filterMGEs(all_genomes_listing_file, outdir, proc_genomes_listing_file, logObject, threads=1, genomad_database=None, genomad_splits=8):
	mgecut_processed_dir = outdir + 'mgecut_processed_genomes/'
	mgecut_tmp_dir = outdir + 'mgecut_tmp/'
	setupDirectories([mgecut_processed_dir, mgecut_tmp_dir])
	
	parallel_jobs_4thread = max(math.floor(threads / 4), 1)
	multi_thread = 4
	if threads < 4:
		multi_thread = threads
		parallel_jobs_4thread = 1

	mge_proc_to_unproc_mapping = {}
	mge_unproc_to_proc_mapping = {}
	mgecut_cmds = []
	genomes = 0
	with open(all_genomes_listing_file) as oaglf:
		for line in oaglf:
			genome = line.strip()
			proc_genome = mgecut_processed_dir + genome.split('/')[-1]
			mge_proc_to_unproc_mapping[proc_genome] = genome
			mge_unproc_to_proc_mapping[genome] = proc_genome
			tmp_dir = mgecut_tmp_dir + genome.split('/')[-1]
			mgecut_cmd = ['mgecut', '-i', genome, '-o', proc_genome, '-d', tmp_dir]
			if genomad_database != None:
				mgecut_cmd += ['-m', 'genomad', '-gd', genomad_database, '-gs', str(genomad_splits), '-c', str(threads)]
			else:
				mgecut_cmd += ['-m', 'phispy', '-c', str(multi_thread)]
			mgecut_cmds.append(mgecut_cmd)
			genomes += 1
	try:
		if genomad_database == None:
			msg = 'Processing %d mgecut commands - using PhiSpy!' % len(mgecut_cmds)
			sys.stdout.write(msg + '\n')
			logObject.info(msg)
			p = multiprocessing.Pool(1)
			for _ in tqdm.tqdm(p.imap_unordered(multiProcess, mgecut_cmds), total=len(mgecut_cmds)):
				pass
			p.close()
		else:
			msg = 'Processing %d mgecut commands - using geNomad!' % len(mgecut_cmds)
			sys.stdout.write(msg + '\n')
			logObject.info(msg)
			p = multiprocessing.Pool(parallel_jobs_4thread)
			for _ in tqdm.tqdm(p.imap_unordered(multiProcess, mgecut_cmds), total=len(mgecut_cmds)):
				pass
			p.close()
	except Exception as e:
		msg = "An error occurred during multiprocessing: %s" % str(e)
		sys.stderr.write(msg + '\n')
		logObject.info(msg)

	glmf_handle = open(proc_genomes_listing_file, 'w')
	processed_genomes = 0
	for f in os.listdir(mgecut_processed_dir):
		mgecut_genome_file = mgecut_processed_dir + f
		glmf_handle.write(mgecut_genome_file + '\n')
		processed_genomes += 1
	glmf_handle.close()

	if processed_genomes != genomes:
		msg = 'Warning: not all genomes were processed by mgecut. %d of %d genomes were processed.' % (processed_genomes, genomes)
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		sys.exit(1)

	return [mge_proc_to_unproc_mapping, mge_unproc_to_proc_mapping]

def createLoggerObject(log_file):
	"""
	Description:
	This function creates a logging object.
	********************************************************************************************************************
	Parameters:
	- log_file: Path to file to which to write logging.
	********************************************************************************************************************
	Returns:
	- logger: A logging object.
	********************************************************************************************************************
	"""

	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger

def closeLoggerObject(logObject):
	"""
	Description:
	This function closes a logging object.
	********************************************************************************************************************
	Parameters:
	- logObject: A logging object.
	********************************************************************************************************************
	"""

	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)

def setupDirectories(directories, automate_flag=False):
	"""
	Description:
	This is a generalizable function to create directories.
	********************************************************************************************************************
	Parameters:
	- dictionaries: A list of paths to directories to create or recreate (after removing).
	********************************************************************************************************************
	"""
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				if automate_flag:
					os.system('rm -fr %s' % d)
				else:	
					response = input("The directory %s already exists, will delete and recreate, is this ok? (yes/no): " % d)
					if response.lower() == 'yes':
						os.system('rm -fr %s' % d)
					else:
						msg = 'Deletion not requested! Exiting ...'
						sys.stderr.write(msg + '\n')
						sys.exit(1)
			os.system('mkdir %s' % d)
	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def runCmd(cmd, logObject, check_files=[], check_directories=[], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	if logObject != None:
		logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr,
						executable='/bin/bash')
		for cf in check_files:
			assert (os.path.isfile(cf))
		for cd in check_directories:
			assert (os.path.isdir(cd))
		if logObject != None:
			logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		if logObject != None:
			logObject.error('Had an issue running: %s' % ' '.join(cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))


def is_fasta(fasta):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				SeqIO.parse(ogf, 'fasta')
		else:
			with open(fasta) as of:
				SeqIO.parse(of, 'fasta')
		return True
	except:
		return False

# Yield successive n-sized 
# chunks from l. 
def divide_chunks(l, n): 
	# function taken from: https://www.geeksforgeeks.org/break-list-chunks-size-n-python/
	# looping till length l 
	for i in range(0, len(l), n):
		yield l[i:i + n]

def compute_n50(inputs):
	input_fastas, output_file = inputs
	output_handle = open(output_file, 'w')
	for input_fasta in input_fastas:
		n50 = n50_calc(input_fasta)
		output_handle.write(input_fasta + '\t' + str(n50) + '\n')
	output_handle.close()

def n50_calc(genome_file):
	#Solution adapted from dinovski:
	#https://gist.github.com/dinovski/2bcdcc770d5388c6fcc8a656e5dbe53c
	lengths = []
	seq = ""	
	if genome_file.endswith('.gz'):
		with gzip.open(genome_file, 'rt') as fasta:
			for line in fasta:
				if line.startswith('>'):
					if seq != "":
						lengths.append(len(seq))
					seq = ""
				else:
					seq += line.strip()
		if seq != "":
			lengths.append(len(seq))
	else:	
		with open(genome_file) as fasta:
			for line in fasta:
				if line.startswith('>'):
					if seq != "":
						lengths.append(len(seq))
					seq = ""
				else:
					seq += line.strip()
		if seq != "":
			lengths.append(len(seq))

	## sort contigs longest>shortest
	all_len=sorted(lengths, reverse=True)
	csum=cumsum(all_len)

	n2=int(sum(lengths)/2)

	# get index for cumsum >= N/2
	csumn2=min(csum[csum >= n2])
	ind=where(csum == csumn2)
	n50 = all_len[int(ind[0])]
	return(n50)


def multiProcess(inputs):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list), with last item in list corresponding to a logging object handle for logging
	progress.
	"""
	input_cmd = inputs
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
	except Exception as e:
		sys.stderr.write(traceback.format_exc())

def filterGenomeUsingPhiSpy(phispy_annot_file, genome_file, filtered_genome_file):
	"""
	Filter genome using phage coordinate predictions from PhiSpy
	"""
	try:
		mge_coords = defaultdict(set)
		with open(phispy_annot_file) as opaf:
			for line in opaf:
				line = line.strip()
				ls = line.split('\t')
				scaffold, start, end = ls[1:4]
				start = int(start)
				end = int(end)
				for pos in range(start, end+1):
					mge_coords[scaffold].add(pos)

		outf = open(filtered_genome_file, 'w')
		with open(genome_file) as ogf:
			for rec in SeqIO.parse(ogf, 'fasta'):
				scaffold = rec.id
				if scaffold in mge_coords:
					nonmge_stretch = ''
					nonmge_stretch_id = 1
					for i, allele in enumerate(str(rec.seq)):
						pos = i + 1
						if pos in mge_coords[scaffold]:
							if nonmge_stretch != '':
								nsid = scaffold + '_' + str(nonmge_stretch_id)
								outf.write('>' + nsid + '\n' + str(nonmge_stretch) + '\n')
								nonmge_stretch_id += 1
							nonmge_stretch = ''
						else:
							nonmge_stretch += allele
					if nonmge_stretch != '':
						nsid = scaffold + '_' + str(nonmge_stretch_id)
						outf.write('>' + nsid + '\n' + str(nonmge_stretch) + '\n')
				else:
					outf.write('>' + rec.id + '_1\n' + str(rec.seq) + '\n')
		outf.close()
	except:
		sys.stderr.write('Issue with filtering genome based on MGE coords determined by PhiSpy.\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)


def filterGenomeUsingGeNomad(genomad_phage_annot_file, genomad_plasmid_annot_file, genome_file, filtered_genome_file):
	"""
	Filter genome using phage coordinate predictions from geNomad
	"""
	try:
		full_phage_or_plasmid_scaffs = set([])
		phage_coords = defaultdict(set)
		with open(genomad_phage_annot_file) as opaf:
			for i, line in enumerate(opaf):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				if ls[3] == 'NA':
					scaffold = ls[0]
					full_phage_or_plasmid_scaffs.add(scaffold)
				else:
					scaffold = '|'.join(ls[0].split('|')[:-1])
					start = int(ls[3].split('-')[0])
					end = int(ls[3].split('-')[1])
					for pos in range(start, end+1):
						phage_coords[scaffold].add(pos)

		with open(genomad_plasmid_annot_file) as opaf:
			for i, line in enumerate(opaf):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				scaffold = ls[0]
				full_phage_or_plasmid_scaffs.add(scaffold)

		outf = open(filtered_genome_file, 'w')
		with open(genome_file) as ogf:
			for rec in SeqIO.parse(ogf, 'fasta'):
				scaffold = rec.id
				if scaffold in full_phage_or_plasmid_scaffs:
					continue
				elif scaffold in phage_coords:
					nonmge_stretch = ''
					nonmge_stretch_id = 1
					for i, allele in enumerate(str(rec.seq)):
						pos = i + 1
						if pos in phage_coords[scaffold]:
							if nonmge_stretch != '':
								nsid = scaffold + '_' + str(nonmge_stretch_id)
								outf.write('>' + nsid + '\n' + str(nonmge_stretch) + '\n')
								nonmge_stretch_id += 1
							nonmge_stretch = ''
						else:
							nonmge_stretch += allele
					if nonmge_stretch != '':
						nsid = scaffold + '_' + str(nonmge_stretch_id)
						outf.write('>' + nsid + '\n' + str(nonmge_stretch) + '\n')
				else:
					outf.write('>' + rec.id + '_1\n' + str(rec.seq) + '\n')
		outf.close()
	except:
		sys.stderr.write('Issue with filtering genome based on MGE coords determined by geNomad.\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
