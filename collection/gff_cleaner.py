### Boas Pucker ###
### b.pucker@tu-bs.de ###


__version__ = "0.11"


__usage__ = """
					python3 gff_cleaner.py
					--in <INPUT_GFF_FILE>
					--out <OUTPUT_GFF_FILE>
					
					optional:
					--feature <mRNA_FEATURE_NAME>[mRNA]
					--ID <ID_TAG_IN_LAST_COLUMN>[ID]
					
					Bug reports and feature requests: b.pucker@tu-bs.de
					"""

import os, sys

# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	if '--feature' in arguments:
		feature = arguments[ arguments.index('--feature')+1 ]
	else:
		feature = "mRNA"
	if '--ID' in arguments:
		ID = arguments[ arguments.index('--ID')+1 ]
	else:
		ID = "ID"

	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len( parts ) > 2:
						if parts[2] == feature:
							start, end = int( parts[3] ), int( parts[4] )
							if start < end:
								subparts = parts[-1].split(';')
								for each in subparts:
									try:
										if each.split('=')[0] == ID:
											out.write( "\t".join( parts[:-1] + [ each ] ) + "\n" )
									except IndexError:
										print( "ERROR: no correct feature ID detected in line - " + line )
				line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
