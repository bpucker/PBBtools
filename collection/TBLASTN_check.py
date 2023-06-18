### Boas Pucker ###
### v0.11 ###
### b.pucker@tu-bs.de ###

__usage__ = """ 
				python3 TBLASTN_check.py
				--out <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS>
				--in <BAIIT_FASTA>
				--ref <REF_GENOME_FASTA>
				--gff <GFF3_FILE>
				
				feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import re, os, sys, subprocess
from operator import itemgetter


# --- end of imports --- #

def load_gene_positions_per_chr( gff_file ):
	"""! @brief load gene positions """
	
	gene_pos_per_chr = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if len( parts ) > 2:
					if parts[2] == "gene":
						try:
							gene_pos_per_chr[ parts[0] ].append( { 'start': int( parts[3] ), 'end': int( parts[4] ), 'id': parts[-1] } )
						except KeyError:
							gene_pos_per_chr.update( { parts[0]: [ { 'start': int( parts[3] ), 'end': int( parts[4] ), 'id': parts[-1] } ] } )
			line = f.readline()
	return gene_pos_per_chr


def load_blast_hits( tblastn_result_file ):
	"""! @brief load best blast hit per query """
	
	blast_hits_per_query = {}
	with open( tblastn_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				blast_hits_per_query[ parts[0] ].append( { 'chr': parts[1], 'start': int( parts[8] ), 'end': int( parts[9] ) } )
			except KeyError:
				blast_hits_per_query.update( { parts[0]: [ { 'chr': parts[1], 'start': int( parts[8] ), 'end': int( parts[9] ) } ] } )
			line = f.readline()
	return blast_hits_per_query


def compare_blast_hits_to_gene_positions( hits, positions, output_file ):
	"""! @brief compare BLAST hits to annotated gene positions """
	
	black_list = []
	with open( output_file, "w" ) as out:
		for hit in hits:
			match = []
			try:
				for pos in positions[ hit['chr'] ]:
					if hit['start'] < pos['end']:
						if hit['end'] > pos['start']:
							match.append( pos['id'] )
			except KeyError:
				pass
			if len( match ) < 1:
				out.write( hit['chr'] + "\t" + str( hit['start'] ) + "\t" + str( hit["end"] ) + "\n" )
			else:
				if match[0] not in black_list:
					out.write( match[0] + "\n" )
					black_list.append( match[0] )


def main( arguments ):
	"""! @brief run everything """
	
	bait_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	ref_file = arguments[ arguments.index('--ref')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	
	tblastn = "tblastn"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	tblastn_result_file = output_folder + "tblastn.txt"
	if not os.path.isfile( tblastn_result_file ):
		p = subprocess.Popen( args= tblastn + " -query " + bait_file + " -subject " + ref_file + " -out " + tblastn_result_file + " -outfmt 6 -evalue 0.00001", shell=True )
		p.communicate()
	
	gene_positions_per_chr = load_gene_positions_per_chr( gff_file )
	
	blast_hits = load_blast_hits( tblastn_result_file )
	
	for each in list( blast_hits.keys() ):
		output_file = output_folder + each + ".txt"
		compare_blast_hits_to_gene_positions( blast_hits[ each ], gene_positions_per_chr, output_file )
	


if '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv and '--gff' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
