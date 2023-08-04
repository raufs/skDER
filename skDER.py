#!/usr/bin/env python

### Program: skDER.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2023, Rauf Salamzade
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# setup is not adherent to pep recommendation for using toml at the moment - might update later
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import shutil
from time import sleep
import argparse
import gzip
import traceback
import logging
import subprocess
import pkg_resources
from Bio import SeqIO
import multiprocessing
import pyfastx


SKDER_DIR = '/'.join(os.path.abspath(__file__).split('/')[:-1]) + '/'
version = pkg_resources.require("skDER")[0].version

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: skder.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	skDER: efficient dynamic & high-resolution dereplication of microbial genomes to select representative genomes.

	skDER relies heavily on advances made by skani for fast ANI estimation while retaining accuracy - thus if you use skDER for your research it is essential to cite skani:
	- "Fast and robust metagenomic sequence comparison through sparse chaining with skani"

	Also please consider citing the lsaBGC manuscript - where a predecessor version of the dynamic dereplication stratedgy employed by skder was first described:
	- "Evolutionary investigations of the biosynthetic diversity in the skin microbiome using lsaBGC"
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-g', '--genomes', nargs='+', help='Genome assembly files in FASTA format (each file should end with either *.fasta, *.fa, or *.fna) [Optional].', required=False, default=[])
	parser.add_argument('-t', '--taxa_name', help='Genus or species identifier from GTDB (currently R214) for which to download genomes for [Optional].', required=False, default=None)
	parser.add_argument('-o', '--output_directory', help='Output directory.', required=True)
	parser.add_argument('-i', '--percent_identity_cutoff', type=float, help="ANI cutoff for dereplication [Default is 99.0].", required=False, default=99.0)
	parser.add_argument('-f', '--aligned_fraction_cutoff', type=float, help="Aligned cutoff threshold for dereplication - only needed by one genome [Default is 90.0].", required=False, default=90.0)
	parser.add_argument('-m', '--max_af_distance_cutoff', type=float, help="Maximum difference for aligned fraction between a pair to automatically disqualify the genome with a higher AF from being a representative.", required=False, default=10.0)
	parser.add_argument('-p', '--skani_triangle_parameters', help="Options for skani triangle. Note ANI and AF cutoffs\nare specified separately and the -E parameter is always\nrequested. [Default is \"\"].", default="", required=False)
	parser.add_argument('-l', '--symlink', action='store_true', help="Symlink representative genomes in results subdirectory instead of performing a copy of the files.", required=False, default=False)
	parser.add_argument('-c', '--cpus', type=int, help="Number of CPUs to use.", required=False, default=1)
	parser.add_argument('-v', '--version', action='store_true', help="Report version of skDER.", required=False, default=False)
	args = parser.parse_args()
	return args

def skder_main():
	# Parse arguments
	myargs = create_parser()

	genomes = myargs.genomes
	taxa_name = None
	if myargs.taxa_name:
		taxa_name = myargs.taxa_name.strip('"')
	outdir = os.path.abspath(myargs.output_directory) + '/'
	percent_identity_cutoff = myargs.percent_identity_cutoff
	aligned_fraction_cutoff = myargs.aligned_fraction_cutoff
	skani_triangle_parameters = myargs.skani_triangle_parameters
	max_af_distance_cutoff = myargs.max_af_distance_cutoff
	cpus = myargs.cpus
	symlink_flag = myargs.symlink
	version_flag = myargs.version

	if version_flag:
		sys.stderr.write('Version of skDER.py being used is: ' + str(version) + '\n')
		sys.exit(0)

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory already exists! Overwriting in 5 seconds, but only where needed will use checkpoint files to skip certain steps...\n")
		sleep(5)
	else:
		setupDirectories([outdir])

	# Create logging object
	log_file = outdir + 'Progress.log'
	logObject = createLoggerObject(log_file)
	parameters_file = outdir + 'Command_Issued.txt'
	sys.stdout.write('Running version %s\n' % version)
	sys.stdout.write("Appending command issued for future records to: %s\n" % parameters_file)
	sys.stdout.write("Logging more details at: %s\n" % log_file)
	logObject.info("\nNEW RUN!!!\n**************************************")
	logObject.info('Running version %s' % version)
	logObject.info("Appending command issued for future records to: %s" % parameters_file)

	parameters_handle = open(parameters_file, 'a+')
	parameters_handle.write(' '.join(sys.argv) + '\n')
	parameters_handle.close()

	all_genomes_listing_file = outdir + 'All_Genomes_Listing.txt'

	if taxa_name != "None" and taxa_name != None:
		# Step 0: parse GTDB information file, get list of Genbank accessions, and perform dry-run with
		# ncbi-genome-download if requested.
		gtdb_listing_file = SKDER_DIR + 'GTDB_R214_Information.txt.gz'

		if not os.path.isfile(gtdb_listing_file):
			sys.stdout.write("GTDB listing file not available, using wget to download it.\n")
			logObject.info("\nGTDB listing file not available, using wget to download it.")
			wget_cmd = ['wget', 'https://github.com/Kalan-Lab/lsaBGC/raw/main/db/GTDB_R214_Information.txt.gz', '-P', outdir]
			gtdb_listing_file = outdir + "GTDB_R214_Information.txt.gz"
			runCmd(wget_cmd, logObject, check_files=[gtdb_listing_file])

		genbank_accession_listing_file = outdir + 'NCBI_Genbank_Accession_Listing.txt'
		sys.stdout.write("--------------------\nStep 0\n--------------------\nBeginning by assessing which genomic assemblies are available for the taxa %s in GTDB and NCBI's Genbank db\n" % taxa_name)
		logObject.info("\n--------------------\nStep 0\n--------------------\nBeginning by assessing which genomic assemblies are available for the taxa %s in GTDB and NCBI's Genbank db" % taxa_name)

		if not os.path.isfile(genbank_accession_listing_file):
			genbank_accession_listing_handle = open(genbank_accession_listing_file, 'w')
			with gzip.open(gtdb_listing_file, 'rt') as ogtdb:
				for line in ogtdb:
					line = line.strip('\n')
					ls = line.split('\t')
					if ls[0] == 'none': continue
					if len(taxa_name.split()) == 1:
						if ls[1] == taxa_name:
							genbank_accession_listing_handle.write(ls[0] + '\n')
					elif len(taxa_name.split()) == 2:
						if ls[2] == taxa_name:
							genbank_accession_listing_handle.write(ls[0] + '\n')
			genbank_accession_listing_handle.close()

		ogalf = open(genbank_accession_listing_file)
		accession_count = len(ogalf.readlines())
		ogalf.close()

		if accession_count == 0:
			sys.stderr.write('Warning: no genomes found to belong the genus or species specified in GTDB.\n')
			logObject.info('Warning: no genomes found to belong the genus or species specified in GTDB.')
		else:
			genome_listing_file = outdir + 'NCBI_Genomes_from_Genbank_for_Taxa.txt'
			if not os.path.isfile(genome_listing_file):
				if accession_count != 0:
					ngd_dry_cmd = ['ncbi-genome-download', '--dry-run', '--section', 'genbank', '-A', genbank_accession_listing_file, 'bacteria', '>', genome_listing_file]
					runCmd(ngd_dry_cmd, logObject, check_files=[genome_listing_file])
				else:
					of = open(genome_listing_file, 'w')
					of.close()

			oglf = open(genome_listing_file)
			genome_count = len(oglf.readlines())
			oglf.close()
			if genome_count == 0:
				sys.stderr.write('Warning: no genomes could be downloaded with ncbi-genome-download.\n')
				logObject.info('Warning: no genomes could be downloaded with ncbi-genome-download.')
			else:
				# Download all genomes in FASTA format & prodigal gene calling
				genomes_directory = outdir + 'ngd_fasta/'
				if not os.path.isfile(all_genomes_listing_file):
					if genome_count != 0:
						ngd_real_cmd = ['ncbi-genome-download', '--formats', 'fasta', '--retries', '2', '--section',
								'genbank', '-A', genbank_accession_listing_file, '-o', genomes_directory,
								'--flat-output',  'bacteria']
						runCmd(ngd_real_cmd, logObject, check_directories=[genomes_directory])
						uncompress_cmds = []
						for f in os.listdir(genomes_directory):
							uncompress_cmds.append(['gunzip', genomes_directory + f])#, logObject])
						p = multiprocessing.Pool(cpus)
						p.map(multiProcess, uncompress_cmds)
						p.close()
						gf_listing_handle = open(all_genomes_listing_file, 'a+')
						for gf in os.listdir(genomes_directory):
							gfile = genomes_directory + gf  
							if not gfile.endswith('.fasta') and not gfile.endswith('.fa') and not gfile.endswith('.fna'): continue
							assert (is_fasta(gfile))
							gf_listing_handle.write(gfile + '\n')
						gf_listing_handle.close()

	if genomes:
		gf_listing_handle = open(all_genomes_listing_file, 'a+')
		for gf in genomes:
			if not gf.endswith('.fasta') and not gf.endswith('.fa') and not gf.endswith('.fna'): continue
			assert(is_fasta(gf))
			gf_listing_handle.write(gf + '\n')
		gf_listing_handle.close()

	number_of_genomes = None
	try:
		assert(os.path.isfile(all_genomes_listing_file))
		with open(all_genomes_listing_file) as oaglf:
			number_of_genomes= len(oaglf.readlines())
		assert(number_of_genomes >= 2)
	except Exception as e:
		logObject.error('Fewer than 2 genomes downloaded / provided. Exiting ...')
		sys.stderr.write('Fewer than 2 genomes downloaded / provided. Exiting ...\n')
		sys.exit(1)

	# calculate N50s
	n50_dir = outdir + 'Assembly_N50s/'
	os.system('mkdir ' + n50_dir)
	n50_inputs = []
	with open(all_genomes_listing_file) as oaglf:
		for line in oaglf:
			genome_path = line.strip()
			genome_file = genome_path.split('/')[-1]
			n50_inputs.append([genome_path, n50_dir + genome_file + '_n50.txt'])
	p = multiprocessing.Pool(cpus)
	p.map(compute_n50, n50_inputs)
	p.close()

	# concatenate N50 results into a single file
	concat_n50_result_file = outdir + 'Concatenated_N50.txt'
	os.system('time find %s -maxdepth 1 -type f | xargs cat >> %s' % (n50_dir, concat_n50_result_file))

	# run skani triangle
	skani_result_file = outdir + 'Skani_Triangle_Edge_Output.txt'
	skani_triangle_cmd = ['skani', 'triangle', '-l', all_genomes_listing_file, '-s', str(percent_identity_cutoff),
			      '--min-af', str(aligned_fraction_cutoff), '-E', skani_triangle_parameters,
			      '-o', skani_result_file]
	runCmd(skani_triangle_cmd, logObject, check_files=[skani_result_file])

	skder_result_file = outdir + 'skDER_Results.txt'
	
	skder_core_prog = SKDER_DIR + 'skDERcore'
	if not os.path.isfile(skder_core_prog):
		skder_core_prog = 'skDERcore'
	skder_core_cmd = [skder_core_prog, skani_result_file, concat_n50_result_file, str(max_af_distance_cutoff), '>', skder_result_file ]
	runCmd(skder_core_cmd, logObject, check_files=[skder_result_file])

	# copy over genomes which are non-redundant to a separate directory
	skani_drep_dir = outdir + 'Dereplicated_Representative_Genomes/'	
	if not os.path.isdir(skani_drep_dir):
		setupDirectories([skani_drep_dir])

	with open(skder_result_file) as osrf:
		for line in osrf:
			genome_path = line.strip()
			try:
				if symlink_flag:
					symlink_file = skani_drep_dir + genome_path.split('/')[-1]
					os.symlink(genome_path, symlink_file)
				else:
					shutil.copy2(genome_path, skani_drep_dir)
				
			except:
				sys.stderr.write('Warning: issues copying over representative genome %s to final dereplicated sub-directory.\n' % genome_path)
				logObject.warning('Issues copying over representative genome %s to final dereplicated sub-directory.' % genome_path)
				
	# close logging object and exit
	logObject.info('******************\nskDER finished!\n******************\nDirectory with dereplicated genomes can be found at: %s' % skani_drep_dir)
	sys.stdout.write('******************\nskDER finished!\n******************\nDirectory with dereplicated genomes can be found at: %s\n' % skani_drep_dir)
	closeLoggerObject(logObject)
	sys.exit(0)

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

def compute_n50(inputs):
	"""
	Uses pyfastx
	"""
	input_fasta, output_file = inputs
	fa = pyfastx.Fasta(input_fasta, build_index=True)
	n50, _ = fa.nl(50)
	output_handle = open(output_file, 'w')
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

if __name__ == '__main__':
	multiprocessing.set_start_method('fork')
	skder_main()
