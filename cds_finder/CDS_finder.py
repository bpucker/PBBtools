### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator and KIPEs ###

__version__ = "v0.2"

__usage__ = """
					python3 CDS_finder.py
					--in <INPUT_FILE>
					--out <OUTPUT_FILE>
					
					optional:
					--len <MIN_CDS_LENGTH>[100]
					--atg <CDS_MUST_START_WITH_ATG>[on]
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""

import os, re, sys, subprocess
from operator import itemgetter

# --- end of imports --- #

def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update( { header: "".join( seq ).upper() } )
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ).upper() } )
	return sequences


def identify_best_cds_candidate( cds_candidates ):
	"""! @brief find the longest CDS candidate """
	
	if len( cds_candidates ) > 0:
		data = []
		for each in cds_candidates:
			data.append( { 'len': len( each ), 'seq': each } )
		return sorted( data, key=itemgetter('len') )[-1]['seq']
	else:
		return ""


def get_cds_candidates_per_frame( seq, atg_status ):
	"""! @brief identify all possible CDS in given reading frame """
	
	candidates = []
	codons = [ seq[i:i+3] for i in range( 0, len( seq ), 3 ) ]
	current_seq = []
	for codon in codons:
		if codon in [ "TAA", "TGA", "TAG" ]:
			current_seq.append( codon )
			candidates.append( "".join( current_seq ) )
			current_seq = []
		else:
			if len( current_seq ) > 0:
				current_seq.append( codon )
			elif atg_status == True and codon == "ATG":
				current_seq.append( codon )
	if len( current_seq ) > 0:
		candidates.append( "".join( current_seq ) )
	return candidates


def find_longest_ORF( assembly, len_cutoff, atg_status ):
	"""! @brief find the longest ORF in contig """
	
	best_cds_candidates = {}
	for key in list( assembly.keys() ):
		cds_candidates = []
		seq = assembly[ key ]
		cds_candidates += get_cds_candidates_per_frame( seq, atg_status )	#frame1
		cds_candidates += get_cds_candidates_per_frame( seq[1:], atg_status )	#frame2
		cds_candidates += get_cds_candidates_per_frame( seq[2:], atg_status )	#frame3
		revseq = revcomp( seq )
		cds_candidates += get_cds_candidates_per_frame( revseq, atg_status )	#frame4
		cds_candidates += get_cds_candidates_per_frame( revseq[1:], atg_status )	#frame5
		cds_candidates += get_cds_candidates_per_frame( revseq[2:], atg_status )	#frame6
		best_cds_candidate = identify_best_cds_candidate( cds_candidates )
		if len( best_cds_candidate ) > len_cutoff:
			best_cds_candidates.update( { key: best_cds_candidate } )
	return best_cds_candidates


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--len' in arguments:
		try:
			min_CDS_len = int( arguments[ arguments.index('--len')+1 ] )
		except ValueError:
			min_CDS_len = 100
	else:
		min_CDS_len = 100
	
	if '--atg' in arguments:
		atg_status = arguments[ arguments.index('--atg')+1 ].upper()
		if atg_status == "OFF":
			atg_status = False
		else:
			atg_status = True
	else:
		atg_status = True
	
	assembly = load_sequences( input_file )
	
	CDS_collection = find_longest_ORF( assembly, min_CDS_len, atg_status )
	
	with open( output_file, "w" ) as out:
		for key in list( CDS_collection.keys() ):
			out.write( '>' + key + "\n" + CDS_collection[key] + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
