### Boas Pucker ###
### b.pucker@tu-bs.de ###
## v0.11 ###


__usage__ = """
python alignment_input_mixer.py
--in <FULL_PATH_TO_INPUT_FILE>
--out <FULL_PATH_TO_OUTPUT_DIR>

bug reports and feature requests: b.pucker@tu-bs.de
					"""

import sys
from random import shuffle

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def main( arguments ):
	"""! @brief shuffles all sequences in a given FASTA file """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	seqs = load_sequences( input_file )
	seq_ids = seqs.keys()
	shuffle( seq_ids )
	with open( output_file, "w" ) as out:
		for ID in seq_ids:
			out.write( '>' + ID + '\n' + seqs[ ID ] + '\n' )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
	
