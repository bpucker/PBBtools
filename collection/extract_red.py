### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.11 ###

__usage__ = """
					python extract_red.py
					--tree <TREE_FILE>
					--seq <SEQUENCE_FILE>
					--out <OUTPUT_FILE>
					
					optional:
					--taxon <taxon_table>
					--color <color_to_extract>
					
					"""


import re, sys, os

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().replace( ".", "_" ).replace( " ", "" )
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update( { header: "".join( seq ) } )
				header = line.strip()[1:].replace( ".", "_" ).replace( " ", "" )
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_taxon_mapping( taxon_file ):
	"""! @brief load taxon mappping """
	
	taxon_mapping = {}
	with open( taxon_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				taxon_mapping[ parts[3] + "_" + parts[2] ].append( parts[1] )
			except KeyError:
				taxon_mapping.update( { parts[3] + "_" + parts[2]: [ parts[1] ] } )
			line = f.readline()
	return taxon_mapping


def main( arguments ):
	"""! @brief run everything """

	tree_file = arguments[ arguments.index('--tree')+1 ]
	seq_file = arguments[ arguments.index('--seq')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	if '--color' in arguments:
		color = arguments[ arguments.index('--color')+1 ]
	else:
		color = "#ff0000"
	
	if '--taxon' in arguments:
		taxon_file = arguments[ arguments.index('--taxon')+1 ]
		taxon_mapping = load_taxon_mapping( taxon_file )
		print len( taxon_mapping.keys() )
	else:
		taxon_mapping = {}

	with open( tree_file, "r" ) as f:
		content = f.read()

	IDs = re.findall( "'{0,1}[a-zA-Z0-9_\-\.@|\[\]\?]+'{0,1}\[&!color=" + color + "\]", content )
	print "number of selected sequences: " + str( len( IDs ) )
	seqs = load_sequences( seq_file )

	clean_IDs = []
	for ID in IDs:
		if "'" in ID:
			clean_IDs.append( ID[1:-18] )
		else:
			clean_IDs.append( ID[:-17] )
	
	counter = 0
	with open( output_file, "w" ) as out:
		for ID in clean_IDs:
			try:
				if "||" in ID:
					out.write( '>' + ID.split('||')[-1] + '\n' + seqs[ ID.split('||')[-1] ] + "\n" )
					counter += 1
				else:
					out.write( '>' + ID + '\n' + seqs[ ID ] + "\n" )
					counter += 1
			except KeyError:
				try:
					status = False
					for each in taxon_mapping[ ID.split('@')[0] ]:
						altID = each + "@" + ID.split('@')[1]
						try:
							out.write( '>' + altID + '\n' + seqs[ altID ] + "\n" )	#" " + ID +
							status = True
						except KeyError:
							pass
					if not status:
						print ID
				except KeyError:
					pass	#print ID
	print "number of successfully extracted seqs: " + str( counter )


if '--tree' in sys.argv and '--seq' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
