### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python fastq2fasta.py
					--fastq <FASTQ_FILE (INPUT)>
					--fasta <FASTA_FILE (OUTPUT)>
					
					WARNING: this script expects gzip compressed files
					"""

import gzip, os, sys

# --- end of imports --- #

def convert_fastq_to_fasta( fastq_file, fasta_file ):
	"""! @brief extraction of header and sequence from fastq file for fasta construction """
	
	with open( fasta_file, "w" ) as out:
		with gzip.open( fastq_file, "rb" ) as f:
			line = f.readline()
			while line:
				out.write(  '>' + line.replace(' ', '_') )
				out.write( f.readline() )
				f.readline()
				f.readline()
				line = f.readline()


def main( arguments ):
	"""! @brief run everything """
	
	fastq_file = arguments[ arguments.index('--fastq')+1 ]
	fasta_file = arguments[ arguments.index('--fasta')+1 ]
	
	convert_fastq_to_fasta( fastq_file, fasta_file )


if '--fastq' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

