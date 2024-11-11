### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.3 ###


__usage__ = """
	python3 aa2cds.py
	--aln <ALIGNMENT_INPUT_FILE>
	--cds <CDS_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	
	bug reports and feature requests: b.pucker@tu-bs.de
	"""

import sys, os

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load sequences of given FASTA file into dictionary with sequence IDs as keys and sequences as values """
	
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


def generate_aligned_cds( aln_seq, cds_seq ):
	"""! @brief convert aligned peptide sequence into aligned corresponding coding sequence """
	
	cds_chunks = list(cds_seq[ 0+i:3+i ] for i in range( 0, len( cds_seq ), 3 ))
	aa_length = len(aln_seq)-aln_seq.count('-')
	cds_length = len( cds_seq )
	if aa_length*3 != cds_length:
		sys.exit( "WARNING: sequence length mismatch\n" + aln_seq + "\n" )
	
	new_seq = []
	for idx, aa in enumerate( aln_seq ):
		if aa == "-":
			new_seq.append( "---" )
		else:
			new_seq.append( cds_chunks[ idx - aln_seq[:idx].count('-') ] )
	return "".join( new_seq )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	out_file = arguments[ arguments.index('--out')+1 ]
	aln_file = arguments[ arguments.index('--aln')+1 ]
	cds_file = arguments[ arguments.index('--cds')+1 ]
	
	cds = load_sequences( cds_file )
	aln = load_sequences( aln_file )
	cds_aln = {}
	for key in list( aln.keys() ):
		aln_seq = aln[ key ]
		try:
			cds_seq = cds[ key ]
			#check if divisble by 3
			if len( cds_seq ) % 3 != 0:
				sys.stdout.write( "invalid CDS sequence (length not multiple of 3): " + key + "(length=" + str( len( cds_seq ) ) + ")\n" )
				sys.stdout.flush()
			else:
				cds_aln.update( { key: generate_aligned_cds( aln_seq, cds_seq ) } )
		except KeyError:
			sys.stdout.write( "missing sequence in CDS: " + key + "\n" )
			sys.stdout.flush()
	with open( out_file, "w" ) as out:
		for key in list( cds_aln.keys() ):
			out.write( '>' + key + '\n' + cds_aln[ key ] + "\n" )


if '--aln' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
