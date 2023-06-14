### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python automatic_SRA_download.py
	--in <FILE_WITH_ONE_RUN_ID_PER_LINE>
	--out <OUTPUT_DIRECTORY>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #

def convert_from_SRA( SRA_ID, fastq_dump, output_directory ):
	"""! @brief starts and handles download of a single SRA entry """
	
	cmd = "".join( [ 	fastq_dump,
								" --split-files --outdir ",
								output_directory,
								" --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ",	#--readids 
								SRA_ID
							] )
	os.popen( cmd )


def main( arguments ):
	"""! @brief runs everything """
	
	information_input_file = arguments[ arguments.index('--in')+1 ]
	prefix = arguments[ arguments.index('--out')+1 ]
	
	if prefix[-1] != '/':
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	fastq_dump = "fastq-dump"
	cache_dir = "/vol/tmp"
	url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"
	
	IDs = []
	with open( information_input_file, "r" ) as f:
		line = f.readline()
		while line:
			IDs.append( line.strip() )
			line = f.readline()
	print "number of IDs: " + str( len( IDs ) )
	
	for idx, ID in enumerate( IDs ):
		print "processing " + str( idx+1 ) + "/" + str( len( IDs ) ) + "\t" + ID
		output_directory = prefix + ID + "/"
		if not os.path.exists( output_directory ):
			os.makedirs( output_directory )
			try:
				convert_from_SRA( ID, fastq_dump, output_directory )
			except:
				print "ERROR: " + ID


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
	

	#--skip-technical prevents downloading of technical reads
	#--readids will add .1 or .2 to read IDs to generate reads with different IDs (breaks BWA!)
	#--read-filter pass will ensure that reads contain sequence and not just Ns
	#--dumpbase ensures that only A, C, G, and T occur in the read and not color space
	#--split-3 generates two files for paired reads and an additional file for singletons
	#--clip removes adapters
