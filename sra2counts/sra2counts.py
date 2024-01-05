### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.2 ###

__usage__ = """
					python3 sra2counts.py
					--sra <SRA_ID_FILE>
					--cds <REFERENCE_CDS_FILE>
					--out <FIGURE_OUTPUT_FILE>
					
					optional:
					--downloader <SRA_DOWNLOAD_SCRIPT>
					--pipeline <KALLISTO_PIPELINE_SCRIPT>
					--merger <KALLISTO_RESULT_MERGER_SCRIPT>
					--cleaner <COUNT_TABLE_FILTER_SCRIPT>
					"""

#REQUIREMENTS in $PATH: fastq-dump, kallisto

import os, sys, subprocess

# --- end of imports --- #

def load_IDs( input_sra_ID_file ):
	"""! @brief load SRA IDs from input file """
	
	IDs = []
	with open( input_sra_ID_file, "r" ) as f:
		content = f.read().strip()
		if "\n" in content:	#support for normal text files
			return content.split('\n')
		elif "\r" in content:	#support for text files created on Mac
			return content.split('\r')
		else:	#single ID i.e. no line break in file
			return [ content ]


def cds_check( cds_file, min_number=10000 ):
	"""! @brief run CDS FASTA file check """
	
	# --- load sequences --- #
	seqs = {}
	with open( cds_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					seqs.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		seqs.update( { header: "".join( seq ) } )	
	
	# --- perform checks --- #
	if len( list( seqs.keys() ) ) < min_number:
		return False
	first_seq = list( seqs.values() )[0].upper()
	if len( first_seq ) > ( first_seq.count( "A" ) + first_seq.count( "C" ) + first_seq.count( "G" ) + first_seq.count( "T" ) + first_seq.count( "N" )):
		return False
	
	return True


def main( arguments ):
	"""! @brief run generation of plots """
	
	input_sra_ID_file = arguments[ arguments.index('--sra')+1 ]
	cds_file = arguments[ arguments.index('--cds')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--downloader' in arguments:
		sra_downloader = arguments[ arguments.index('--downloader')+1 ]
	else:
		sra_downloader = "automatic_SRA_download.py"
	
	if '--pipeline' in arguments:
		kallisto_pipeline = arguments[ arguments.index('--pipeline')+1 ]
	else:
		kallisto_pipeline = "kallisto_pipeline3.py"
	
	if '--merger' in arguments:
		count_table_merger = arguments[ arguments.index('--merger')+1 ]
	else:
		count_table_merger = "merge_kallisto_output3.py"
	if '--cleaner' in arguments:
		count_table_cleaner = arguments[ arguments.index('--cleaner')+1 ]
	else:
		count_table_cleaner = "filter_RNAseq_samples.py"
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- check that only <=10 SRA IDs are included in the input file --- #
	IDs = load_IDs( input_sra_ID_file )
	if len( IDs ) > 10:
		IDs = IDs[:10]
		sys.stdout.write( "WARNING: too many SRA IDs submitted in this batch. Reduced to 10: " + ",".join( IDs ) + "\n" )
		sys.stdout.flush()
	sra_ID_file = output_folder + "SRA_IDs_to_do.txt"
	with open( sra_ID_file, "w" ) as out:
		out.write( "\n".join( IDs ) + "\n" )
	
	# --- check validity of cds file --- #
	cds_status = cds_check( cds_file )
	if not cds_status:
		sys.exit( "ERROR: CDS file is invalid." )
	
	# --- download RNA-seq samples --- #
	fastq_dir = output_folder + "fastq/"
	cache_dir = output_folder + "cache/"	#could be replaced by some kind of tmp folder
	if not os.path.exists( fastq_dir ):	#only run this part if the FASTQ folder does not exist yet
		os.makedirs( fastq_dir )
		os.makedirs( cache_dir )
		cmd = [ "python2", sra_downloader, "--in", sra_ID_file, "--out", fastq_dir, "--cache", cache_dir ]
		p = subprocess.Popen( args= " ".join( cmd ), shell=True )
		p.communicate()
	
	# --- run kallisto across all samples --- #
	kallisto_tmp_dir = output_folder + "kallisto_tmp/"
	kallisto_result_dir = output_folder + "kallisto_results/"
	if not os.path.exists( kallisto_result_dir ):
		os.makedirs( kallisto_tmp_dir )
		os.makedirs( kallisto_result_dir )
		cmd = [ "python3", kallisto_pipeline, "--cds", cds_file, "--reads", fastq_dir, "--tmp", kallisto_tmp_dir, "--out", kallisto_result_dir ]
		p = subprocess.Popen( args= " ".join( cmd ), shell=True )
		p.communicate()
	
	# --- merge kallisto result files --- #
	result_folder = output_folder + "results/"	#only download this folder
	tpm_file = result_folder + "TPMs.txt"
	count_file = result_folder + "counts.txt"
	if not os.path.exists( result_folder ):
		os.makedirs( result_folder )
		cmd = [ "python3", count_table_merger, "--in", kallisto_result_dir, "--tpms", tpm_file, "--counts", count_file ]
		p = subprocess.Popen( args= " ".join( cmd ), shell=True )
		p.communicate()
	
	# --- filter count table --- #
	clean_tpm_file = result_folder + "TPMs.clean.txt"
	if not os.path.isfile( clean_tpm_file ):
		cmd = [ "python3", count_table_cleaner, "--tpms", tpm_file, "--counts", count_file, "--out", clean_tpm_file ]
		p = subprocess.Popen( args= " ".join( cmd ), shell=True )
		p.communicate()
	
	
	#delete temporary folders to reduce download archive size !!!
	#os.rmdir( fastq_dir )
	#os.rmdir( cache_dir )
	#os.rmdir( kallisto_tmp_dir )
	#os.rmdir( kallisto_result_dir )


if '--sra' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
