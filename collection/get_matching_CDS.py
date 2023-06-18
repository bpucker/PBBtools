### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python3 get_matching_CDS.py
					--in <INPUT_FASTA_FILE>
					--out <OUTPUT_FASTA_FILE>
					--data <DATA_FOLDERS_COMMA_SEPARATED>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, glob, re, os

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
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def main( arguments ):
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	data = arguments[ arguments.index('--data')+1 ]
	
	if "," in data:
		data_folders = data.split(',')
	else:
		data_folders = [ data ]
	
	data_file_names = []
	for folder in data_folders:
		data_file_names += glob.glob( folder + "*.fasta" ) + glob.glob( folder + "*.fa" )
	
	seqs_of_interest = load_sequences( input_file )
	print( "number of peptide sequences: " + str( len( seqs_of_interest.keys() ) ) )
	
	cds_seqs = {}
	for filename in data_file_names:
		cds_seqs.update( load_sequences( filename ) )
	
	print( "number of coding sequences: " + str( len( cds_seqs.keys() ) ) )
	
	with open( output_file, "w" ) as out:
		for key in seqs_of_interest.keys():
			try:
				out.write( '>' + key + '\n' + cds_seqs[ key ] + "\n" )
			except KeyError:
				try:
					if "_" in key:
						alt_key = "_".join( key.split('_')[1:] )
						out.write( '>' + alt_key + '\n' + cds_seqs[ alt_key ] + "\n" )
					else:
						print( key )
				except KeyError:
					print( key )


if '--in' in sys.argv and '--out' in sys.argv and '--data' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
