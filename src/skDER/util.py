import os
import sys
import logging
import traceback
import pyfastx
import subprocess
from Bio import SeqIO
import gzip

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
	logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr,
						executable='/bin/bash')
		for cf in check_files:
			assert (os.path.isfile(cf))
		for cd in check_directories:
			assert (os.path.isdir(cd))
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
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
	"""
	Uses pyfastx
	"""
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
		except:
			pass
		output_handle.write(input_fasta + '\t' + str(n50) + '\n')
	output_handle.close()

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
