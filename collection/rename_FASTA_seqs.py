### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

__usage__ = """
					python rename_FASTA_seqs.py
					--in <INPUT_FASTA_FILE>
					--out <OUTPUT_FASTA_FILE>
					--taxon <TAXON_FILE>
					
					optional:
					--lineage <ADDS_LINEAGE_PREFIX>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, glob, re, os

# --- end of imports --- #

def load_taxon_mapping( taxon_file, lineage ):
	
	mapping = {}
	with open( taxon_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if lineage:
				mapping.update( { parts[1]: parts[3] + "_" + parts[2] } )
			else:
				mapping.update( { parts[1]: parts[2] } )
			line = f.readline()
	return mapping



def main( arguments ):
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	taxon_file = arguments[ arguments.index('--taxon')+1 ]
	
	if '--lineage' in arguments:
		lineage = True
	else:
		lineage = False
	
	mapping_table = load_taxon_mapping( taxon_file, lineage )
	
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == ">":
					if "@" in line:
						#print line[1:].split('@')[0]
						try:
							out.write( '>' + mapping_table[ line[1:].split('@')[0] ] + "@" + line.split('@')[1] )
						except KeyError:
							#print "ERROR: " + line.strip()
							out.write( line )
					else:
						out.write( line )
				else:
					out.write( line )
				line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv and '--taxon' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
