### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
						python generate_gff3_from_peptides.py
						--tblastn <TBLASTN_RESULT_FILE>
						--pep <PEPTIDE_FILE>
						--assembly <TRANSCRIPTOME_ASSEMBLY_FILE>
						--out <OUTPUT_FOLDER>
						"""

from operator import itemgetter
import matplotlib.pyplot as plt
import time, os, sys

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq.replace( "X", "" ) } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq.replace( "X", "" ) } )
	return sequences


def load_blast_results( tblastn_result_file ):
	"""! @brief load all tBLASTn results from file """
	
	hits = {}
	
	with open( tblastn_result_file ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] == parts[1] and parts[2] == "100.000":	#take only perfect hits
				try:
					entry = hits[ parts[0] ]
					if float( parts[-1] ) > entry['score']:
						del hits[ parts[0] ]
						if int( parts[6] ) < int( parts[7] ):
							if int( parts[8] ) < int( parts[9] ):
								hits.update( { parts[0]: { 'id': parts[0], 'q_start': int( parts[6] ), 'q_end': int( parts[7] ), 's_start': int( parts[8] ), 's_end': int( parts[9] ), 'orientation': True, 'score': float( parts[-1] )} } )
							else:
								hits.update( { parts[0]: { 'id': parts[0], 'q_start': int( parts[6] ), 'q_end': int( parts[7] ), 's_start': int( parts[8] ), 's_end': int( parts[9] ), 'orientation': False, 'score': float( parts[-1] ) } } )
						else:
							print "ERROR: wrong orientation of query!"
				except KeyError:
					if int( parts[6] ) < int( parts[7] ):
						if int( parts[8] ) < int( parts[9] ):
							hits.update( { parts[0]: { 'id': parts[0], 'q_start': int( parts[6] ), 'q_end': int( parts[7] ), 's_start': int( parts[8] ), 's_end': int( parts[9] ), 'orientation': True, 'score': float( parts[-1] ) } } )
						else:
							hits.update( { parts[0]: { 'id': parts[0], 'q_start': int( parts[6] ), 'q_end': int( parts[7] ), 's_start': int( parts[8] ), 's_end': int( parts[9] ), 'orientation': False, 'score': float( parts[-1] ) } } )
					else:
						print "ERROR: wrong orientation of query!"
			line = f.readline()
	print "number of hits: " + str( len( hits.keys() ) )
	return hits


def revcomp( seq ):
	"""! @brief construct reverse complement for given sequence """
	
	seq = seq.lower()
	rev_comp_dict = { 'a': 't', 't':'a', 'g':'c', 'c':'g', 'n':'n' }
	new_seq = []
	for nt in seq:
		try:
			new_seq.append( rev_comp_dict[ nt ] )
		except KeyError:
			new_seq.append( 'n' )
	return "".join( new_seq[::-1] ).upper()


def load_genetic_code( genetic_code_file ):
	"""! @brief load genetic code from file """
	
	genetic_code = {}
	
	with open( genetic_code_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			codons = parts[1].split(', ')
			for codon in codons:
				genetic_code.update( { codon: parts[0] } )
			line = f.readline()
	return genetic_code


def translate( seq, genetic_code ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	seq = seq.upper()
	
	peptide = []
	
	for i in range( int( len( seq ) / 3.0 ) ):
		codon = seq[i*3:i*3+3]
		try:
			peptide.append( genetic_code[ codon ] )
		except:
			peptide.append( "*" )
	return "".join( peptide )


def get_start_and_end( transcript, peptide, genetic_code ):
	"""! @brief get position of peptide within translated transcript """
	
	transcript = transcript.upper()
	peptide = peptide.upper()
	
	pep = translate( transcript, genetic_code )
	if peptide in pep:
		return 3*(pep.index( peptide ))+1, 3*(pep.index( peptide ))+3*len( peptide )
	else:
		pep = translate( transcript[1:], genetic_code )
		if peptide in pep:
			return 3*(pep.index( peptide ))+2 , 3*(pep.index( peptide ))+3*len( peptide )+1
		else:
			pep = translate( transcript[2:], genetic_code )
			if peptide in pep:
				return 3*(pep.index( peptide )) + 3, 3*(pep.index( peptide ))+3*len( peptide )+2
			
			else:
				pep = translate( revcomp( transcript ), genetic_code )
				if peptide in pep:
					return "ERROR", 3*pep.index( peptide ) + 1, 3*pep.index( peptide )+3*len( peptide )+2
				else:
					pep = translate( revcomp( transcript )[1:] , genetic_code )
					if peptide in pep: 
						return "ERROR", 3*pep.index( peptide ) + 2, 3*pep.index( peptide )+3*len( peptide )+2
					else:
						pep = translate( revcomp( transcript )[2:], genetic_code )
						if peptide in pep:
							return "ERROR", 3*pep.index( peptide ) + 3, 3*pep.index( peptide )+3*len( peptide )+2


def construct_gff3_file_based_on_tblastn_results( blast_results, peptides, transcripts, gff3_file, genetic_code ):
	"""! @brief construct GFF3 file based on information from tBLASTn of peptide sequences vs. encoding contigs """
	
	with open( gff3_file, "w" ) as out:
		for key in sorted( peptides.keys() ):
			ok_status = True
			new_line = [ key, ".", "CDS" ]
			try:
				entry = blast_results[ key ]
				
				#feature on forward strand
				if entry['orientation']:
					new_line = [ entry['id'], ".", "CDS" ]
					if entry['q_start'] == 1 and entry['q_end'] == len( peptides[ key ] ):
						new_line.append( str( entry['s_start'] ) )
						new_line.append( str( min( [ entry['s_end']+3, len( transcripts[ key ] ) ] ) ) )
					else:
						try:
							start, end = get_start_and_end( transcripts[ key ], peptides[ key ], genetic_code )
						except:
							error, start, end = get_start_and_end( transcripts[ key ], peptides[ key ], genetic_code )
							print "ERROR" + entry['id']
							ok_status = False
						new_line.append( str( start ) )
						new_line.append( str( min( [ end+3,  len( transcripts[ key ] ) ] ) ) )
					new_line.append( "." )
					new_line.append( "+" )
				
				#feature on reverse strand
				else:
					if entry['q_start'] == 1 and entry['q_end'] == len( peptides[ key ] ):
						new_line.append( str( max( [ 1, entry['s_end']-3 ] ) ) )
						new_line.append( str( entry['s_start'] ) )
					else:
						try:
							start, end = get_start_and_end( revcomp( transcripts[ key ] ), peptides[ key ], genetic_code )
						except:
							error, start, end = get_start_and_end( revcomp( transcripts[ key ] ), peptides[ key ], genetic_code )
							print "ERROR" + entry['id']
							ok_status = False
						new_start = max( [ 1, len( transcripts[ key ]  ) - end + 1 - 3 ] )
						new_end = len( transcripts[ key ]  ) - start + 1
						new_line.append( str( new_start ) )
						new_line.append( str( new_end ) )
					new_line.append( "." )
					new_line.append( "-" )
			except KeyError:
				transcript = transcripts[ key ]
				peptide = peptides[ key ]
				fw_results = get_start_and_end( transcript, peptide, genetic_code )
				if len( fw_results ) == 2:
					start, end = fw_results
					new_line.append( str( start ) )
					new_line.append( str( end ) )
					new_line.append( "." )
					new_line.append( "+" )
				else:
					start, end = get_start_and_end( revcomp( transcript ), peptide, genetic_code )
					new_start = max( [ 1, len( transcripts[ key ]  ) - end + 1 - 3 ] )
					new_end = len( transcripts[ key ]  ) - start + 1
					new_line.append( str( new_start ) )
					new_line.append( str( new_end ) )
					new_line.append( "." )
					new_line.append( "-" )
			new_line.append( "." )
			new_line.append( "ID=" + key )
			if ok_status:
				out.write( "\t".join( new_line ) + '\n' )


def construct_CDS_file( transcripts, gff3_file, CDS_file ):
	"""! @brief construct CDS file based on transcript sequences and information from GFF3 file """
	
	CDS = {}
	
	with open( CDS_file, "w" ) as out:
		with open( gff3_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if parts[6] == "+":
						seq =  transcripts[ parts[0] ][ int( parts[3] )-1: int( parts[4] ) ] 
						out.write( '>' + parts[0] + '\n' + seq + '\n' )
						CDS.update( { parts[0]: seq } )
					else:
						seq = revcomp( transcripts[ parts[0] ][ int( parts[3] )-1: int( parts[4] ) ] )
						out.write( '>' + parts[0] + '\n' + seq + '\n' )
						CDS.update( { parts[0]: seq } )
				line = f.readline()
	return CDS


def main( arguments ):
	
	tblastn_result_file = arguments[ arguments.index('--tblastn')+1 ]
	peptide_file = arguments[ arguments.index('--pep')+1 ]
	assembly_file = arguments[ arguments.index('--assembly')+1 ]
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	if '--name' in arguments:
		name = arguments[ arguments.index('--name')+1 ]
	else:
		name = "x"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	gff3_file = output_folder + name + ".cds.gff3"
	CDS_file = output_folder + name + ".cds.fasta"
	
	transcripts = load_sequences( assembly_file )
	print "number of transcripts: " + str( len( transcripts.keys() ) )
	peptides = load_sequences( peptide_file )
	print "number of peptides: " + str( len( peptides.keys() ) )
	blast_results = load_blast_results( tblastn_result_file )
	genetic_code = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'}
	
	construct_gff3_file_based_on_tblastn_results( blast_results, peptides, transcripts, gff3_file, genetic_code )
	CDS = construct_CDS_file( transcripts, gff3_file, CDS_file )

	
if '--tblastn' in sys.argv and '--pep' in sys.argv and '--assembly' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
