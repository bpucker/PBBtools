### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python3 avg_cov_per_contig.py
					--in <COV_FILE>
					--out <OUTPUT_FILE>
					feature requests and bug reports: b.pucker@tu-bs.de
					"""

import os, sys

# --- end of imports --- #


def calculate_values( cov_file ):
	"""! @brief calculate the average coverage per contig """
	
	avg_cov_per_contig = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		prev_contig = line.split('\t')[0]
		values = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_contig:
				if len( values ) > 0:
					avg_cov_per_contig.update( { prev_contig: sum( values ) / len( values) } )
				else:
					avg_cov_per_contig.update( { prev_contig: 0 } )
				prev_contig = parts[0] + ""
				values = []
			values.append( float( parts[-1] ) )
			line = f.readline()
		if len( values ) > 0:
			avg_cov_per_contig.update( { prev_contig: sum( values ) / len( values) } )
		else:
			avg_cov_per_contig.update( { prev_contig: 0 } )
	return avg_cov_per_contig


def main( arguments ):
	
	cov_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	avg_cov_per_contig = calculate_values( cov_file )
	
	with open( output_file, "w" ) as out:
		for key in list( avg_cov_per_contig.keys() ):
			out.write( key + "\t" + str( avg_cov_per_contig[ key ] ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )


