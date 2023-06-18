### Boas Pucker ###
### b.pucker@tu-bs.de ###

__cite__ = """Pucker, B., 2023. xxx"""

__version__ = "v0.1"

__usage__ = """
					Extract single copy BUSCO peptide sequences (""" + __version__ + """
					python3 get_single_copy_BUSCOs.py
					--busco <BUSCO_RESULT_TSV_FILE>
					--pepin <PEPTIDE_SEQUENCE_INPUT_FILE>
					--pepout <PEPTIDE_SEQUENCE_OUTPUT_FILE>
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""

import os, sys

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	busco_result_file = arguments[ arguments.index('--busco')+1 ]
	pep_input_file = arguments[ arguments.index('--pepin')+1 ]
	pep_output_file = arguments[ arguments.index('--pepout')+1 ]
	
	peps = load_sequences( pep_input_file )
	single_copy_ids = []
	with open( busco_result_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[1] == "Complete":
					single_copy_ids.append( parts[2] )
			line = f.readline()

	print( "number of BUSCOs: " + str(len( single_copy_ids ) ) )

	with open( pep_output_file, "w" ) as out:
		for ID in single_copy_ids:
			out.write( '>' + ID + '\n' + peps[ ID ] + "\n" )


if '--busco' in sys.argv and '--pepin' in sys.argv and '--pepout' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
