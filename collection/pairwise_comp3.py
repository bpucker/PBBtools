### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.15 ###

__usage__ = """
					python3 pairwise_comp3.py
					--in <FULL_PATH_TO_INPUT_FASTA>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""


import sys, os, subprocess

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
					sequences.update( { header: "".join( seq ).upper() } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ).upper() } )	
	return sequences


def get_identity( name1, name2, seqs, tmp_folder ):
	"""! @brief get identity between two sequences """
	
	seq_file = tmp_folder + name1 + "_vs_" + name2 + ".fasta"
	with open( seq_file, "w" ) as seqout:
		seqout.write( '>' + name1 + "\n" + seqs[ name1 ] + "\n>" + name2 + "\n" + seqs[ name2 ] )
	mafft_file = tmp_folder + name1 + "_vs_" + name2 + ".fasta.aln"
	if not os.path.isfile( mafft_file ):
		p = subprocess.Popen( args= "mafft " + seq_file + " > " + mafft_file, shell=True )
		p.communicate()
	
	seq1, seq2 = load_sequences( mafft_file ).values()
	counter = 0
	for idx, aa in enumerate( seq1 ):
		if aa == seq2[ idx ]:
			counter += 1
	value = round( 100.0 * counter / len( seq1 ), 2 )
	return str( value ) + "%"


def main( arguments ):
	"""! @brief run everything """
	
	input_fasta = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	tmp_folder = output_folder + "tmp/"
	if not os.path.exists( tmp_folder ):
		os.makedirs( tmp_folder )
	
	seqs = load_sequences( input_fasta )
	names = sorted( list( seqs.keys() ) )
	
	output_file = output_folder + "summary.txt"
	with open( output_file, "w" ) as out:
		out.write( "\t" + "\t".join( names ) + "\n" )
		for idx1, name1 in enumerate( names ):
			new_line = [ name1 ]
			for idx2, name2 in enumerate( names ):
				if idx2 < idx1:
					new_line.append( "-" )
				elif idx2 == idx1:
					new_line.append( "100%" )
				else:
					new_line.append( get_identity( name1, name2, seqs, tmp_folder ) )
			out.write( "\t".join( map( str, new_line ) ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
