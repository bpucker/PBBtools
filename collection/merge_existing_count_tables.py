### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.2 ###

__usage__ = """
					python3 merge_existing_count_tables.py
					--in <COMMA_SEPARATED_LIST_OF_COUNT_TABLES>
					--out <OUTPUT_FILE>
					"""

import sys, os

# --- end of imports --- #

def load_exp( exp_file ):
	"""! @brief load expression values """
	
	exp = {}
	genes = []
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip().split('\t')
		if headers[0] == "gene":
			headers = headers[1:]
		for header in headers:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes.append( parts[0] )
			for idx, val in enumerate( parts[1:] ):
				exp[ headers[ idx ] ].update( { parts[0]: val } )
			line = f.readline()
	return exp, headers, genes


def main( arguments ):
	"""! @brief run generation of plots """
	
	input_data = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	# --- collect count table input files --- #
	input_files = []
	for each in input_data.split(','):
		if os.path.isfile( each ):
			input_files.append( each )
	
	# --- load all data from individual count tables --- #
	data = []
	for filename in input_files:
		e, h, g = load_exp( filename )
		data.append( { 'exp': e, 'samples': h, 'genes': g, 'id': filename } )
	
	# --- check for consistency --- #
	gene_number = len( data[0]['genes'] )
	for entry in data:
		if len( entry['genes'] ) != gene_number:
			print( "WARNING: GENE NUMBERS DO NOT MATCH " + entry['id'] )
	
	# --- merge everything --- #
	all_exp = {}
	for entry in data:
		all_exp.update( entry['exp'] )
	
	genes = data[0]['genes']
	samples = sorted( list( all_exp.keys() ) )
	
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + "\n" )
		for gene in genes:
			new_line = [ gene ]
			for sample in samples:
				new_line.append( all_exp[ sample ][ gene ] )
			out.write( "\t".join( new_line ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
