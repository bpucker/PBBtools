### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.25 ###

__usage__ = """
					python generate_subset.py
					--in <FULL_PATH_TO_EXPRESSION_INPUT_FILE>
					--out <FULL_PATH_TO_EXPRESSION_OUTPUT_FILE>
					--accessions <FULL_PATH_TO_FILE_WITH_ACCESSION_NUMBERS>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os

# --- end of imports --- #


def load_expression_values( data_file ):
	"""! @brief load expression values """
	
	exp_data = {}
	with open( data_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		for header in headers:
			exp_data.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, value in enumerate( parts[1:] ):
				exp_data[ headers[ idx ] ].update( { parts[0]:  float( value ) } )
			line = f.readline()
	return exp_data


def load_accessions( sample_file ):
	"""! @brief load all accession numbers from given file """
	
	with open( sample_file, "r" ) as f:
		content = f.read().strip()
	
	if ',' in content:
		parts = content.split( ',' )
	elif ";" in content:
		parts = content.split( ';' )
	elif "\n" in content:
		 parts = content.split( '\n' )
	elif "\t" in content:
		 parts = content.split( '\n' )
	elif " " in content:
		 parts = content.split( ' ' )
	else:
		sys.exit( "ERROR: no accession number detected in file (uncommon separator?)" )
	
	accessions = []
	for part in parts:
		if len( part.strip() ) > 5:
			accessions.append( part.strip() )
	if len( accessions ) > 0:
		return accessions
	else:
		sys.exit( "ERROR: no accession number detected in file" )


def generate_output_file( accessions, exp_data, output_file ):
	"""! @brief generate output file """
	
	valid_accessions = []
	for accession in accessions:
		try:
			exp_data[ accession ]
			valid_accessions.append( accession )
		except KeyError:
			pass
	
	sys.stdout.write( "number of valid accessions: " + str( len( valid_accessions ) ) + '\n' )
	sys.stdout.flush()
	
	if len( valid_accessions ) < 1:
		sys.exit( "ERROR: no provided accession IDs matched to existing data set." )
	
	genes = exp_data.values()[0].keys()
	
	with open( output_file, "w" ) as out:
		out.write( "genes\t" + "\t".join( valid_accessions ) + '\n' )
		for gene in genes:
			new_line = [ gene ]
			for accession in valid_accessions:
				new_line.append( exp_data[ accession ][ gene ] )
			out.write( "\t".join( map( str, new_line ) )+"\n" )


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	sample_file = arguments[ arguments.index('--accessions')+1 ]
	
	if output_file[0] == "/":	#only generate new folder if absolute path is given
		output_dir = "/".join( output_file.split('/')[:-1] )
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir )
	
	accessions = load_accessions( sample_file )
	
	sys.stdout.write( "number of detected accessions: " + str( len( accessions ) ) + '\n' )
	sys.stdout.flush()
	
	exp_data = load_expression_values( input_file )
	
	generate_output_file( accessions, exp_data, output_file )


if '--in' in sys.argv and '--out' in sys.argv and '--accessions' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
