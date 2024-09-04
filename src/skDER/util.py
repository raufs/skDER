import os
import sys
import logging
import traceback
import subprocess
import numpy
from Bio import SeqIO
import gzip
from collections import defaultdict

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

def setupDirectories(directories):
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
				os.system('rm -rf %s' % d)
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

"""
# No longer used.
def compute_n50_pyfastx(inputs):
	# Uses pyfastx
	input_fastas, output_file, index_locally_flag = inputs
	output_handle = open(output_file, 'w')
	for input_fasta in input_fastas:
		genome_file = input_fasta
		if index_locally_flag:
			local_file = '/'.join(os.path.abspath(output_file).split('/')[:-1]) + '/' + '.'.join(input_fasta.split('/')[-1].split('.')[:-1]) + '.fasta'
			os.symlink(input_fasta, local_file)
			genome_file = local_file
		fa = pyfastx.Fasta(genome_file, build_index=True)
		n50, _ = fa.nl(50)
		try:
			index_file = genome_file + '.fxi'
			os.remove(index_file)
			if index_locally_flag:
				assert(input_fasta != genome_file)
				os.remove(genome_file)
		except:
			pass
		output_handle.write(input_fasta + '\t' + str(n50) + '\n')
	output_handle.close()
"""

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
	csum=numpy.cumsum(all_len)

	n2=int(sum(lengths)/2)

	# get index for cumsum >= N/2
	csumn2=min(csum[csum >= n2])
	ind=numpy.where(csum == csumn2)
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
