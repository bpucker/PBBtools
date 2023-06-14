### Boas Pucker ###
### b.pucker@tu-bs.de ###

__version__ = "v0.1"

__usage__ = """
					Collect sequences identified by collect_best_BLAST_hits.py (""" + __version__ + """)
					python3 collect_and_rename_seqs.py
					--in <INPUT_FOLDER>
					--out <OUTPUT_FILE>
					"""

import os, sys, glob

# --- end of imports --- #

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
		if len( "".join( seq ) ) > 0:
			sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	"""! @brief run everything """

	input_folder = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if input_folder[-1] != "/":
		input_folder += "/"

	relevant_files = glob.glob( input_folder + "*/best_hit_sequences.fasta" )

	seq_collection = {}
	for filename in relevant_files:
		spec = filename.split('/')[-2]
		seqs = load_sequences( filename )
		if len( list( seqs.keys() ) ) > 0:
			for key in list( seqs.keys() ):
				seq_collection.update( { spec + "@" + key: seqs[ key ] } )
		else:
			print( "WARNING: no sequences detected for species " + spec )

	with open( output_file, "w" ) as out:
		for key in list( seq_collection.keys() ):
			out.write( '>' + key + "\n" + seq_collection[ key ] + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
