### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

# Functions from kallisto_pipeline3.py and merge_kallisto3.py have been re-used here

__usage__ = """
	python3 FASTQ_processor.py
	--in <FILE_WITH_ONE_RUN_ID_PER_LINE>
	--ref <REFERENCE_CDS_FILE>
	--out <OUTPUT_DIRECTORY>
	--tmp <TMP_DIRECTORY>
	
	optional:
	--kallisto <KALLISTO_PATH>[kallisto]
	--fastqdump <FASTQ-DUMP_PATH>[fastq-dump]
	--prefetch <PREFETCH_PATH>[prefetch]
	--cpus <NUMBER_OF_THREADS>[2]
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, glob, subprocess, gzip

# --- end of imports --- #


def get_data_for_job_to_run( folder, output_folder, index_file, tmp_dir ):
	"""! @brief collect all infos to run jobs """
	
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
	output_dir = tmp_dir + ID + "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir ) 
	tmp_result_file = output_dir + "abundance.tsv"
	final_result_file = output_folder + ID + ".tsv"
	if os.path.isfile( final_result_file ):
		status = False
	if os.path.isfile( final_result_file + ".gz" ):
		status = False
	if status:
		return [ { 'r1': read_file1, 'r2': read_file2, 'out': output_dir, 'index': index_file, 'tmp': tmp_result_file, 'fin': final_result_file, "ID": ID } ]
	else:
		print( read_file1 )
		print(PE_status  )
		print( read_file2 )
		return []


def run_kallisto( job, kallisto, threads ):
	"""! @brief run all jobs in list """
	
	sys.stdout.write( "running job " + job["ID"] + "\n" )
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


def convert_from_SRA( SRA_ID, tmp_dir, prefetch, fastq_dump ):
	"""! @brief starts and handles download of a single SRA entry """
	
	#run prefetch
	print( prefetch + " --output-directory " + tmp_dir + " " + SRA_ID )
	p = subprocess.Popen( args= prefetch + " --output-directory " + tmp_dir + " " + SRA_ID, shell=True )
	p.communicate()
	
	#run fastq-dump
	print(  fastq_dump + " --gzip --outdir " + tmp_dir + " " + tmp_dir + SRA_ID + "/" + SRA_ID+".sra" )
	p = subprocess.Popen( args= fastq_dump + " --gzip --outdir " + tmp_dir + SRA_ID + "/" + " " + tmp_dir + SRA_ID + "/" + SRA_ID+".sra", shell=True )
	p.communicate()


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


def generate_mapping_table( gff_file ):
	"""! @brief generate transcript to gene mapping table """
	
	transcript2gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA", "transcript" ]:
					try:
						ID = parts[-1].split(';')[0].split('=')[1]
						parent = parts[-1].split('arent=')[1]
						if ";" in parent:
							parent = parent.split(';')[0]
						transcript2gene.update( { ID: parent } )
					except:
						print(line)
			line = f.readline()
	return transcript2gene


def map_counts_to_genes( transcript2gene, counts ):
	"""! @brief map transcript counts to parent genes """
	
	error_collector = []
	gene_counts = {}
	for key in counts.keys():
		try:
			gene_counts[ transcript2gene[ key ] ] += counts[ key ]
		except KeyError:
			try:
				gene_counts.update( { transcript2gene[ key ]: counts[ key ] } )
			except KeyError:
				error_collector.append( key )
				gene_counts.update( { key: counts[ key ] } )
	if len( error_collector ) > 0:
		sys.stdout.write( "number of unmapped transcripts: " + str( len( error_collector ) ) + "\n" )
		sys.stdout.flush()
	return gene_counts


def generate_output_file( output_file, data ):
	"""! @brief generate output file for given data dictionary """
	
	samples = list( sorted( list( data.keys() ) ) )
	
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
		for gene in list(sorted(list(data.values())[0].keys())):
			new_line = [ gene ]
			for sample in samples:
				new_line.append( data[ sample ][ gene ] )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )


def merge_kallisto_results( count_table_folder, counts_output_file, tpm_output_file ):
	"""! @brief run everything """
	
	counttables = glob.glob( count_table_folder + "*.tsv" ) + glob.glob( count_table_folder + "*.tsv.gz" )
	sys.stdout.write( "number of detected counttables: " + str( len( counttables ) ) + "\n" )
	sys.stdout.flush()

	count_data = {}
	tpm_data = {}
	for filename in counttables:
		ID = filename.split('/')[-1].split('.')[0]
		counts, tpms = load_counttable( filename )
		#TPM are available and could be processed in the same way
		count_data.update( { ID: counts } )
		tpm_data.update( { ID: tpms } )
	
	generate_output_file( counts_output_file, count_data )
	generate_output_file( tpm_output_file, tpm_data )


def main( arguments ):
	"""! @brief runs everything """
	
	sra_accession_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	tmp_dir = arguments[ arguments.index('--tmp')+1 ]
	ref_file = arguments[ arguments.index('--ref')+1 ]
	
	if output_folder[-1] != '/':
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	count_table_folder = output_folder + "count_tables/"
	if not os.path.exists( count_table_folder ):
		os.makedirs( count_table_folder )
	
	if tmp_dir[-1] != '/':
		tmp_dir += "/"
	if not os.path.exists( tmp_dir ):
		os.makedirs( tmp_dir )
	
	if '--kallisto' in arguments:
		kallisto = arguments[ arguments.index('--kallisto')+1 ]
	else:
		kallisto = "kallisto"
	
	if '--cpus' in arguments:
		try:
			threads = int( arguments[ arguments.index( '--cpus' )+1 ] )
		except:
			threads = 2
	else:
		threads = 2	
	
	if '--prefetch' in arguments:
		prefetch = arguments[ arguments.index( '--prefetch' )+1 ]
	else:
		prefetch = "prefetch"
	
	if '--fastqdump' in arguments:
		fastq_dump = arguments[ arguments.index( '--fastqdump' )+1 ]
	else:
		fastq_dump = "fastq-dump"
	
	# ---- load SRA FASTQ accessions from input file --- #
	IDs = []
	with open( sra_accession_file, "r" ) as f:
		line = f.readline()
		while line:
			IDs.append( line.strip() )
			line = f.readline()
	sys.stdout.write( "number of IDs: " + str( len( IDs ) ) + "\n" )
	sys.stdout.flush()
	
	# --- prepare kallisto reference --- #
	index_file = tmp_dir + "index"
	if not os.path.isfile( index_file ):
		cmd1 = " ".join( [ kallisto, "index", "--index="+index_file, "--make-unique", ref_file ] )
		p = subprocess.Popen( args= cmd1, shell=True )
		p.communicate()
	
	# --- run downloads and analyze individual files --- #
	for idx, ID in enumerate( IDs ):
		sys.stdout.write( "processing " + str( idx+1 ) + "/" + str( len( IDs ) ) + "\t" + ID + "\n" )
		sys.stdout.flush()
		
		fastq_folder = tmp_dir + ID + "/"
		try:
			convert_from_SRA( ID, tmp_dir, prefetch, fastq_dump )
		except:
			sys.stdout.write( "ERROR: while downloading " + ID + "\n" )
			sys.stdout.flush()
		
		# --- run kallisto on FASTQ --- #
		job = get_data_for_job_to_run( fastq_folder, count_table_folder, index_file, tmp_dir )
		sys.stdout.write( str( job ) + "\n" )
		sys.stdout.flush()	
		if len( job ) == 1:
			run_kallisto( job[0], kallisto, threads )	#just run one job that is contained in list
		else:
			sys.stdout.write( "ERROR_WITH_JOB\n" )
			sys.stdout.flush()	
		
		p = subprocess.Popen( args= "rm -r " + fastq_folder, shell=True )
		p.communicate()
	
	# --- merge all count tables --- #
	counts_output_file = output_folder + "counts.txt"
	tpm_output_file = output_folder + "tpms.txt"
	merge_kallisto_results( count_table_folder, counts_output_file, tpm_output_file )


if '--in' in sys.argv and '--out' in sys.argv and '--tmp' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
	

	#--skip-technical prevents downloading of technical reads
	#--readids will add .1 or .2 to read IDs to generate reads with different IDs (breaks BWA!)
	#--read-filter pass will ensure that reads contain sequence and not just Ns
	#--dumpbase ensures that only A, C, G, and T occur in the read and not color space
	#--split-3 generates two files for paired reads and an additional file for singletons
	#--clip removes adapters
