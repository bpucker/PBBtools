### Boas Pucker ###
### b.pucker@tu-bs.de ###

__version__ = "v0.1"

__usage__ = """
					python3 rename_ncbi_collection.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					"""

import re, sys, os


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]

	seqs = load_sequences( input_file )

	with open( output_file, "w" ) as out:
		for key in list( seqs.keys() ):
			try:
				ID = re.findall( "[A-Z0-9_]+", key )[0]
				spec = re.findall( "\[[a-zA-Z\.\ ]+\]", key )[-1]
				out.write( '>' + spec[1:-1].replace(" ", "_").replace( ".", "_" ).replace( "__", "_" ) + "@" + ID + "\n" + seqs[ key ] + "\n" )
			except IndexError:
				print( "ERROR: " + key )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
