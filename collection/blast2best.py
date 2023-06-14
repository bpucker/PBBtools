### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python blast2best.py
					--in <FULL_PATH_TO_INPUT_FILE>
					--out  <FULL_PATH_TO_OUTPUT_FILE>
					
					optional:
					--num <NUMBER_OF_BEST_HITS>
					"""

import os, sys
from operator import itemgetter

# --- end of imports --- #

def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				best_hits[ parts[0] ].append( { 'score': float( parts[-1] ), 'line': line } )
			except:
				best_hits.update( { parts[0]: [ { 'score': float( parts[-1] ), 'line': line } ] } )
			line = f.readline()
	
	final_sorted_hits = {}
	for key in sorted( best_hits.keys() ):
		data = sorted( best_hits[ key ], key=itemgetter('score') )[::-1]
		final_sorted_hits.update( { key: data } )
	return final_sorted_hits


def main( arguments ):
	
	input_file = arguments[ arguments.index( '--in' ) + 1 ]
	output_file = arguments[ arguments.index( '--out' ) + 1 ]
	
	if '--num' in arguments:
		number_of_hits = int( arguments[ arguments.index( '--num' ) + 1 ] )
	else:
		number_of_hits = 1
	
	hits = load_best_blast_hit( input_file )
	with open( output_file, "w" ) as out:
		for key in sorted( hits.keys() ):
			if len( hits[ key ] ) > number_of_hits:
				for each in hits[ key ][ : number_of_hits ]:
					out.write( each[ 'line' ] )
			else:
				for each in hits[ key ]:
					out.write( each[ 'line' ] )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
