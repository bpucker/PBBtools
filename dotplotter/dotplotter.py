### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.1 ###

__usage__ = """
					python dotplotter.py
					option1:
					--seq1 <SEQ1>
					--seq2 <SEQ2>
					--out <FIGURE_OUTPUT_FILE>
					
					option2:
					--in1 <FASTA_FILENAME1>
					--in2 <FASTA_FILENAME2>
					--out <FIGURE_OUTPUT_FILE>
					
					optional:
					--kmer <KMER_SIZE>[21]
					--name1 <NAME_OF_SEQ1>[seq1]
					--name2 <NAME_OF_SEQ2>[seq2]
					"""

import os, sys, re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# --- end of imports --- #

def generate_kmers( seq1, kmer ):
	"""! @brief generate kmers """
	
	kmers = []
	if len( seq1 ) >= kmer:
		for i in range( len( seq1 )-kmer ):
			kmers.append( seq1[ i:i+kmer ] )
	else:
		sys.exit( "ERROR: k-mer size exceeds sequence length. Use smaller k-mer size." )
	return kmers


def get_all_kmer_positions( kmers, seq1, seq2 ):
	"""! @brief calculate x and y positions for plotting based on k-mer positions """
	
	xpositions, ypositions = [], []
	for each in kmers:
		xpos = [ x.start() for x in re.finditer( each, seq1 ) ]
		ypos = [ y.start() for y in re.finditer( each, seq2 ) ]
		for x in xpos:
			for y in ypos:
				xpositions.append( x )
				ypositions.append( y )
	return xpositions, ypositions


def generate_dotplot( x_kmer_positions, y_kmer_positions, figfile, name1, name2 ):
	"""! @brief generate dotplots """
	
	fig, ax = plt.subplots()
	
	ax.plot( x_kmer_positions, y_kmer_positions, markersize=1, marker=".", linestyle=" ", color="green" )
	
	ax.set_xlabel( name1 )
	ax.set_ylabel( name2 )
	
	fig.tight_layout()
	fig.savefig( figfile, dpi=300 )


def load_seq_from_file( seqfile ):
	"""! @brief load sequence from given FASTA file """
	
	with open( seqfile, "r" ) as f:
		header = f.readline().strip()[1:]
		seq = []
		line = f.readline()
		while line:
			if line[0] == ">":
				return "".join( seq ).upper(), header
			seq.append( line.strip() )
			line = f.readline()
		return "".join( seq ).upper(), header


def main( arguments ):
	"""! @brief run generation of plots """
	
	if '--seq1' in arguments:
		seq1 = arguments[ arguments.index('--seq1')+1 ].upper()
		seq2 = arguments[ arguments.index('--seq2')+1 ].upper()
	elif '--in1' in arguments:
		input_file1 = arguments[ arguments.index('--in1')+1 ]
		input_file2 = arguments[ arguments.index('--in2')+1 ]
		seq1, fasta1_seq_name = load_seq_from_file( input_file1 )
		seq2, fasta2_seq_name = load_seq_from_file( input_file2 )		
	
	figfile = arguments[ arguments.index('--out')+1 ]
	
	if '--kmer' in arguments:
		kmer = int( arguments[ arguments.index('--kmer')+1 ] )
	else:
		kmer = 21
	
	if '--name1' in arguments:
		name1 = arguments[ arguments.index('--name1')+1 ]
	else:
		try:
			name1 = fasta1_seq_name + ""
		except:
			name1 = "seq1"
	
	if '--name2' in arguments:
		name2 = arguments[ arguments.index('--name2')+1 ]
	else:
		try:
			name2 =  fasta2_seq_name + ""
		except:
			name2 = "seq2"
	
	kmers = generate_kmers( seq1, kmer )
	x_kmer_positions, y_kmer_positions = get_all_kmer_positions( kmers, seq1, seq2 )
	
	generate_dotplot( x_kmer_positions, y_kmer_positions, figfile, name1, name2 )


if '--seq1' in sys.argv and '--seq2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
elif '--in1' in sys.argv and '--in2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
