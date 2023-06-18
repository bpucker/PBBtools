### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

__usage__ = """
					python rename_seqs_in_tree.py
					--in <INPUT_TREE_FILE>
					--out <OUTPUT_TREE_FILE>
					--taxon <TAXON_FILE>
					
					optional:
					--lineage <INCLUDES_A_LINEAGE_PREFIX>
					"""


import re, os, sys

# --- end of imports --- #

def load_spec_name_mapping_table( taxon_file, origin ):
	"""! @brief load species name mapping table (include origin in name if specified) """
	
	mapping_table = {}
	with open( taxon_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if origin:
				mapping_table.update( { parts[1]: parts[3] + "_" + parts[2] } )
			else:
				mapping_table.update( { parts[1]: parts[2] } )
			line = f.readline()
	return mapping_table


def main( arguments ):
	"""! @brief replace all sequence ID in tree file with sequence names """

	tree_input_file = arguments[ arguments.index( '--in' )+1 ]
	tree_output_file = arguments[ arguments.index( '--out' )+1 ]
	taxon_file = arguments[ arguments.index( '--taxon' )+1 ]
	
	if '--lineage' in arguments:
		origin = True
	else:
		origin = False

	mapping_table = load_spec_name_mapping_table( taxon_file, origin )
	
	with open( tree_input_file, "r" ) as f:
		tree = f.read()
	tree_parts = tree.split('(')
	
	for idx, part in enumerate( tree_parts ):
		subparts = part.split(',')
		for idx2, subpart in enumerate( subparts ):
			IDs = re.findall( "[a-zA-Z]+[a-zA-Z0-9_\.\-]+", subpart )
			if len( IDs ) > 0:
				try:
					subpart = subpart.replace( IDs[0], mapping_table[ IDs[0] ] + "||" + IDs[0] )
				except KeyError:
					pass	#print IDs[0]
			subparts[ idx2 ] = subpart
		tree_parts[ idx ] = ",".join( subparts )
	
	with open( tree_output_file, "w" ) as out:
		out.write( "(".join( tree_parts ) )


if '--in' in sys.argv and '--out' in sys.argv and '--taxon' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
