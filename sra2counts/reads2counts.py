### Boas Pucker ###
### pucker@uni-bonn.de ###
### Based on previously published scripts: 10.1371/journal.pone.0280155 ###

__version__ = """v0.14	"""

__usage__ = """
					python3 reads2counts.py (""" +  __version__	+ """)
					--in <FILE_WITH_ONE_RUN_ID_PER_LINE>
					--cds <CDS_FILE>
					--out <OUTPUT_DIRECTORY>
					
					optional:
					--fastqdump <FULL_PATH_TO_FASTQ_DUMP>[fastq-dump]
					--kallisto <FULL_PATH_TO_KALLISTO>[kallisto]
					--cpus <NUMBER_OF_CPUS_TO_USE>[10]
					--min <MIN_PERCENT_EXPRESSION_ON_TOP100>[10]
					--max <MAX_PERCENT_EXPRESSION_ON_TOP100>[80]
					--mincounts <MIN_READ_NUMBER>[1000000]
					
					bug reports and feature requests: pucker@uni-bonn.de
					"""

import os, sys, glob, time, subprocess
try:
	import gzip
except ImportError:
	pass
try:
	import matplotlib.pyplot as plt
except ImportError:
	pass

# --- end of imports --- #


def get_data_for_jobs_to_run( folder, counttable_output_folder, index_file, tmp_cluster_folder ):
	"""! @brief collect all infos to run jobs """
	
	jobs_to_do = []
	ID = folder.split('/')[-2]
	status = True
	
	# --- get read file --- #
	PE_status = True
	SRA = False
	read_file1 = folder + ID + "_R1_001.fastq.gz"
	if not os.path.isfile( read_file1 ):
		#print "ERROR: file missing - " + read_file1
		PE_status = False
		if not PE_status:
			read_file1 = folder + ID + "_pass_1.fastq.gz"
			if os.path.isfile( read_file1 ):
				PE_status = True
				SRA = True
				read_file2 = folder + ID + "_pass_2.fastq.gz"
				if not os.path.isfile( read_file2 ):
					#print "ERROR: file missing - " + read_file2
					PE_status = False
			else:
				read_file1 = folder + ID + "_1.fastq.gz"
				if os.path.isfile( read_file1 ):
					PE_status = True
					SRA = True
					read_file2 = folder + ID + "_2.fastq.gz"
					if not os.path.isfile( read_file2 ):
						#print "ERROR: file missing - " + read_file2
						PE_status = False
				else:
					read_file1 = folder + ID + "_R1.fq.gz"
					if os.path.isfile( read_file1 ):
						PE_status = True
						SRA = True
						read_file2 = folder + ID + "_R2.fq.gz"
						if not os.path.isfile( read_file2 ):
							#print "ERROR: file missing - " + read_file2
							PE_status = False
					else:
						read_file1 = folder + ID + "_1.fq.gz"
						if os.path.isfile( read_file1 ):
							PE_status = True
							SRA = True
							read_file2 = folder + ID + "_2.fq.gz"
							if not os.path.isfile( read_file2 ):
								#print "ERROR: file missing - " + read_file2
								PE_status = False
						else:
							read_file1 = folder + ID + "_1.clean.fq.gz"
							if os.path.isfile( read_file1 ):
								PE_status = True
								SRA = True
								read_file2 = folder + ID + "_2.clean.fq.gz"
								if not os.path.isfile( read_file2 ):
									#print "ERROR: file missing - " + read_file2
									PE_status = False
							else:
								try:
									read_file1 = glob.glob( folder + "*_R1_001.fastq.gz" )[0]
									if os.path.isfile( read_file1 ):
										PE_status = True
										SRA = True
										read_file2 = glob.glob( folder + "*_R2_001.fastq.gz" )[0]
										if not os.path.isfile( read_file2 ):
											#print "ERROR: file missing - " + read_file2
											PE_status = False
								except:
									try:
										read_file1 = glob.glob( folder + "*_f1.fq.gz" )[0]
										if os.path.isfile( read_file1 ):
											PE_status = True
											SRA = True
											read_file2 = glob.glob( folder + "*_r2.fq.gz" )[0]
											if not os.path.isfile( read_file2 ):
												#print "ERROR: file missing - " + read_file2
												PE_status = False
									except:
										pass
	if not SRA:
		read_file2 = folder + ID + "_R2_001.fastq.gz"
		if not os.path.isfile( read_file2 ):
			#print "ERROR: file missing - " + read_file2
			PE_status = False
	if not PE_status:
		read_file1 = folder + ID + "_R1_001.fastq.gz"
		if not os.path.isfile( read_file1 ):
			read_file1 = folder + ID + ".fastq.gz"
			if not os.path.isfile( read_file1 ):
				read_file1 = folder + ID + "_1.fastq.gz"
				if not os.path.isfile( read_file1 ):
					read_file1 = folder + ID + "_R1.fq.gz"
					if not os.path.isfile( read_file1 ):
						read_file1 = folder + ID + "_pass.fastq.gz"
						if not os.path.isfile( read_file1 ):
							status = False
		read_file2 = False
	
	# --- get reference for quantification --- #
	output_dir = tmp_cluster_folder + ID + "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir ) 
	tmp_result_file = output_dir + "abundance.tsv"
	final_result_file = counttable_output_folder + ID + ".tsv"
	if os.path.isfile( final_result_file ):
		status = False
	if os.path.isfile( final_result_file + ".gz" ):
		status = False
	if status:
		jobs_to_do.append( { 'r1': read_file1, 'r2': read_file2, 'out': output_dir, 'index': index_file, 'tmp': tmp_result_file, 'fin': final_result_file, "ID": ID } )
	return jobs_to_do


def job_executer( jobs_to_run, kallisto, threads ):
	"""! @brief run all jobs in list """
	
	for idx, job in enumerate( jobs_to_run ):
		sys.stdout.write( "running job " + str( idx+1 ) + "/" + str( len( jobs_to_run ) ) + " - " + job["ID"] + "\n" )
		sys.stdout.flush()
		
		# --- run analysis --- #
		if job['r2']:
			cmd2 = " ".join( [ kallisto, "quant", "--index="+job['index'], "--output-dir="+job['out'], "--threads "+str( threads ), job['r1'], job['r2'] ] )
		else:
			cmd2 = " ".join( [ kallisto, "quant", "--index="+job['index'], "--single -l 200 -s 100", "--output-dir="+job['out'], "--threads "+str( threads ), job['r1'] ] )
		p = subprocess.Popen( args= cmd2, shell=True )
		p.communicate()
		
		# --- copy result file --- #
		p = subprocess.Popen( args= "cp " + job["tmp"] + " " + job["fin"] , shell=True )
		p.communicate()
		
		# --- compress result file (count table) --- #
		p = subprocess.Popen( args= "gzip " + job['fin'], shell=True )
		p.communicate()


def kallisto_quantification( cds_file, single_read_file_folders, kallisto, threads, counttable_output_folder, tmp_cluster_folder ):
	"""! @brief run everything """
	
	# --- prepare jobs to run --- #
	index_file = tmp_cluster_folder + "index"
	jobs_to_run = get_data_for_jobs_to_run( single_read_file_folders, counttable_output_folder, index_file, tmp_cluster_folder )
	sys.stdout.write( "number of jobs to run: " + str( len( jobs_to_run ) ) + "\n" )
	sys.stdout.flush()
	
	# --- run jobs --- #
	job_executer( jobs_to_run, kallisto, threads )


def load_counttable( counttable ):
	"""! @brief load data from counttable """
	
	counts = {}
	tpms = {}
	if counttable[-4:] == ".tsv":
		with open( counttable, "r" ) as f:
			f.readline()	#remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				counts.update( { parts[0]: float( parts[3] ) } )
				tpms.update( { parts[0]: float( parts[4] ) } )
				line = f.readline()
	else:
		with gzip.open( counttable, "rb" ) as f:
			f.readline()	#remove header
			line = f.readline().decode("utf-8")
			while line:
				parts = line.strip().split('\t')
				counts.update( { parts[0]: float( parts[3] ) } )
				tpms.update( { parts[0]: float( parts[4] ) } )
				line = f.readline().decode("utf-8")
	return counts, tpms


def generate_output_file( output_file, data ):
	"""! @brief generate output file for given data dictionary """
	
	samples = list( sorted( list( data.keys() ) ) )
	
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "GeneID" ] + samples ) + '\n' )
		if len( list(data.values()) ) > 0:
			for gene in list(sorted(list(data.values())[0].keys())):
				new_line = [ gene ]
				for sample in samples:
					new_line.append( data[ sample ][ gene ] )
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
		else:
			sys.stdout.write( "No sample files detected for generation of merged counttable.\n" )
			sys.stdout.flush()


def merge_kallisto_output( data_input_dir, counts_output_file, tpm_output_file ):
	"""! @brief run everything """
	
	counttables = glob.glob( data_input_dir + "*.tsv" ) + glob.glob( data_input_dir + "*.tsv.gz" )
	sys.stdout.write( "number of detected counttables: " + str( len( counttables ) ) + "\n" )
	sys.stdout.flush()

	count_data = {}
	tpm_data = {}
	for filename in counttables:
		ID = filename.split('/')[-1].split('.')[0]
		counts, tpms = load_counttable( filename )
		count_data.update( { ID: counts } )
		tpm_data.update( { ID: tpms } )
	
	if counts_output_file:
		generate_output_file( counts_output_file, count_data )
	if tpm_output_file:
		generate_output_file( tpm_output_file, tpm_data )


def load_all_TPMs( exp_file ):
	"""! @brief load all values from given TPM file """
	
	data = {}
	genes = []
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip()
		if "\t" in headers:
			headers = headers.split('\t')
		else:
			headers = [ headers ]
		if headers[0] == "gene":
			headers = headers[1:]
		elif headers[0] == exp_file.split('/')[-1][:10]:
			headers = headers[1:]
		for header in headers:
			data.update( { header: [] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes.append( parts[0] )
			for idx, val in enumerate( parts[1:] ):
				data[ headers[ idx ] ].append( float( val ) )
			line = f.readline()
	return data, genes


def  load_black_IDs( black_list_file ):
	"""! @brief load IDs from given black list """
	
	black_list = {}
	with open( black_list_file, "r" ) as f:
		line = f.readline()
		while line:
			if len( line ) > 3:
				black_list.update( { line.strip(): None } )
			line = f.readline()	
	return black_list


def filter_RNAseq( tpm_file, count_file, output_file, min_cutoff, max_cutoff, min_counts ):
	"""! @brief run everything """
	
	black_list = {}
	
	# --- run analysis of all data in folder/file --- #
	doc_file = output_file + ".doc"
	valid_samples = []
	with open( doc_file, "w" ) as out:
		out.write( "SampleName\tPercentageOfTop100\tPercentageOfTop500\tPercentageOfTop1000\n" )
		TPM_data, genes = load_all_TPMs( tpm_file )
		count_data, genes = load_all_TPMs( count_file )
		for key in sorted( list( TPM_data.keys() ) ):
			new_line = [ key ]
			selection = sorted( TPM_data[ key ] )
			counts = sum( count_data[ key ] )	#calculate counts per library
			if counts >= min_counts:	#check for sufficient library size
				try:	#check for ID presence on black list
					black_list[ key ]
					new_line.append( "ID on black list" )
					out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
				except KeyError:
					try:
						val = 100.0 * sum( selection[-100:] ) / sum( selection )
					except ZeroDivisionError:
						val = 0
					new_line.append( val )
					if min_cutoff < val < max_cutoff:
						valid_samples.append( key )
					if len( selection ) > 500 and val > 0:
						new_line.append( 100.0 * sum( selection[-500:] ) / sum( selection ) )
					else:
						new_line.append( "n/a" )
					if len( selection ) > 1000 and val > 0:
						new_line.append( 100.0 * sum( selection[-1000:] ) / sum( selection ) )
					else:
						new_line.append( "n/a" )
					out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
			else:
				new_line.append( "insufficient counts: " + str( counts ) )
				out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
	
	sys.stdout.write( "number of valid sample: " + str( len( valid_samples ) ) + "\n" )
	sys.stdout.write( "number of invalid sample: " + str( len( TPM_data.keys() ) - len( valid_samples ) ) + "\n" )
	sys.stdout.flush()
	
	# --- generate output file --- #
	if len( valid_samples ) > 0:
		with open( output_file, "w" ) as out:
			out.write( "gene\t" + "\t".join( valid_samples )+ "\n" )
			for idx, gene in enumerate( genes ):
				new_line = [ gene ]
				for sample in valid_samples:
					new_line.append( TPM_data[ sample ][ idx ] )
				out.write( "\t".join( list( map( str, new_line ) ) ) + "\n" )
	else:
		sys.stdout.write( "WARNING: no valid samples in data set!" )
		sys.stdout.flush()
	
	# --- generate figure --- #
	fig_file = output_file + ".pdf"
	values = []
	with open( doc_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				values.append( float( parts[1] ) )
			except ValueError:
				pass
			line = f.readline()
	
	values = [ x for x in values if str(x) != 'nan' ]
	
	try:
		fig, ax = plt.subplots()
		
		ax.hist( values, bins=100, color="green" )
		ax.set_xlabel( "Percentage of expression on top100 genes" )
		ax.set_ylabel( "Number of analyzed samples" )
		
		fig.savefig( fig_file )
	except:
		sys.stdout.write( "WARNING: figure generation error. Is matplotlib installed?\n" )
		sys.stdout.flush()


def main( arguments ):
	"""! @brief runs everything """
	
	information_input_file = arguments[ arguments.index('--in')+1 ]
	prefix = arguments[ arguments.index('--out')+1 ]
	cds_file = arguments[ arguments.index( '--cds' )+1 ]
	
	if '--cpus' in arguments:
		try:
			threads = int( arguments[ arguments.index( '--cpus' )+1 ] )
		except:
			threads = 10
	else:
		threads = 10
	
	if '--kallisto' in arguments:
		kallisto = arguments[ arguments.index( '--kallisto' )+1 ]
	else:
		kallisto = "kallisto"
	
	if "--fastqdump" in arguments:
		fastq_dump = arguments[ arguments.index( '--fastqdump' )+1 ]
	else:
		fastq_dump = "fastq-dump"
	
	if '--cache' in arguments:
		cache_dir = arguments[ arguments.index( '--cache' )+1 ]
	else:
		cache_dir = "/vol/tmp"
	
	if '--min' in arguments:
		min_cutoff = int( arguments[ arguments.index('--min')+1 ] )
	else:
		min_cutoff = 10	#value in percent (min read proportion assigned to top100)
	if '--max' in arguments:
		max_cutoff = int( arguments[ arguments.index('--max')+1 ] )
	else:
		max_cutoff = 80	#value in percent (max read proportion assigned to top100)
	
	if '--mincounts' in arguments:
		min_counts = int( arguments[ arguments.index('--mincounts')+1 ] )
	else:
		min_counts = 1000000	#min counts for RNAseq samples to be kept
	
	url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"
	
	if prefix[-1] != '/':
		prefix += "/"
	
	fastq_folder = prefix + "FASTQs/"
	counttable_output_folder = prefix + "output/"
	tmp_cluster_folder = prefix + "tmp/"
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	if fastq_folder[-1] != '/':
		fastq_folder += "/"
	if not os.path.exists( fastq_folder ):
		os.makedirs( fastq_folder )
	
	if counttable_output_folder[-1] != "/":
		counttable_output_folder += "/"
	if not os.path.exists( counttable_output_folder ):
		os.makedirs( counttable_output_folder )
	
	if tmp_cluster_folder[-1] != "/":
		tmp_cluster_folder += "/"
	if not os.path.exists( tmp_cluster_folder ):
		os.makedirs( tmp_cluster_folder )

	IDs = []
	with open( information_input_file, "r" ) as f:
		line = f.readline()
		while line:
			IDs.append( line.strip() )
			line = f.readline()
	sys.stdout.write( "number of IDs: " + str( len( IDs ) ) + "\n" )
	sys.stdout.flush()
	
	index_file = tmp_cluster_folder + "index"
	# --- generate index --- #
	if not os.path.isfile( index_file ):
		cmd1 = " ".join( [ kallisto, "index", "--index="+index_file, "--make-unique", cds_file ] )
		p = subprocess.Popen( args= cmd1, shell=True )
		p.communicate()
		sys.stdout.write( "Index file is ready. You may start more workers now.\n" )
		sys.stdout.flush()
	
	for idx, ID in enumerate( IDs ):
		sys.stdout.write( "processing " + str( idx+1 ) + "/" + str( len( IDs ) ) + "\t" + ID + "\n" )
		sys.stdout.flush()
		
		# --- check if sample was already done (count table) --- #
		count_table = counttable_output_folder + ID + ".tsv"
		compressed_count_table = counttable_output_folder + ID + ".tsv.gz"
		if os.path.exists( count_table ) + os.path.exists( compressed_count_table  ) < 1:	#check if a count table file exists
			output_directory = fastq_folder + ID + "/"
			if not os.path.exists( output_directory ):	#no folder for FASTQs of current SRA sample
				os.makedirs( output_directory )
				try:
					cmd = "".join( [ 	fastq_dump,
									" --split-files --outdir ",
									output_directory,
									" --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ",	#--readids 
									ID
								] )
					p = subprocess.Popen( args= cmd, shell=True )
					p.communicate()
					
					# --- run kallisto analysis --- #
					kallisto_quantification( cds_file, output_directory, kallisto, threads, counttable_output_folder, tmp_cluster_folder )
					
				except:
					sys.stdout.write("ERROR (unknown issue): " + ID + "\n" )
					sys.stdout.flush()
				
				# --- remove FASTQ folder --- #	
				try:
					p = subprocess.Popen( args= "rm -r " + output_directory, shell=True )
					p.communicate()
					try:
						p = subprocess.Popen( args= "rm -r " + tmp_cluster_folder + ID + "/", shell=True )
						p.communicate()
					except:
						sys.stdout.write("ERROR (TMP folder deletion failed): " + ID + "\n" )
						sys.stdout.flush()
				except:
					sys.stdout.write("ERROR (FASTQ folder deletion failed): " + ID + "\n" )
					sys.stdout.flush()
	
	# --- merge count tables --- #
	counts_output_file = prefix + "counts_table.txt"
	tpm_output_file = prefix + "tpm_table.txt"
	merge_kallisto_output( counttable_output_folder, counts_output_file, tpm_output_file )
	
	# --- filter RNAseq samples --- #
	clean_output_file = prefix + "tpm_table.clean.txt"
	filter_RNAseq( tpm_output_file, counts_output_file, clean_output_file, min_cutoff, max_cutoff, min_counts )


if '--in' in sys.argv and '--out' in sys.argv and '--cds' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
	

#--skip-technical prevents downloading of technical reads
#--readids will add .1 or .2 to read IDs to generate reads with different IDs (breaks BWA!)
#--read-filter pass will ensure that reads contain sequence and not just Ns
#--dumpbase ensures that only A, C, G, and T occur in the read and not color space
#--split-3 generates two files for paired reads and an additional file for singletons
#--clip removes adapters
