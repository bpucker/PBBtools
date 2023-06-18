### Boas Pucker ###
### boas.pucker@rub.de ###
### v0.15 ###

__usage__ = """
					python get_DEGs_from_DESeq2_output.py
					--in <DESEQ2_RESULT_FILE>
					--out <OUTPUT_FILE>
					
					optional:
					--l2fc <LOG2_FOLD_CHANGE>[1]
					--p <ADJUSTED_P_VALUE>[0.05]
					"""

import os, sys
import scipy.stats as s
import matplotlib.pyplot as plt
from operator import itemgetter

# --- end of imports --- #

def load_genes_of_interest( input_file, l2fc_cutoff, adj_pvalue_cutoff ):
	"""! @brief get genes of interest from given count table """
	
	gois = []
	with open( input_file, "r" ) as f:
		header = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split(' ')
			if parts[-1] != "NA":
				if float( parts[-1] ) < adj_pvalue_cutoff:
					if abs( float( parts[2] ) ) > l2fc_cutoff:
						gois.append( { 'ID': parts[0].replace( '"', '' ), 'baseMean': float( parts[1] ), 'l2fc': float( parts[2] ), 'padj': float( parts[-1] ) } )
			line = f.readline()
	return gois


def load_annotation( annotation_file ):
	"""! @brief load all content from annotation file """
	
	mapping_table = {}
	with open( annotation_file, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 1:
				mapping_table.update( { parts[0]: ";".join( parts[1:] ) } )
			line = f.readline()
	return mapping_table


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index( '--anno' )+1 ]
		anno = load_annotation( anno_file )
	else:
		anno = {}
	
	if '--l2fc' in arguments:
		try:
			l2fc_cutoff = float( arguments[ arguments.index( '--l2fc' )+1 ] )
		except:
			print "ERROR: l2fc input not recognized ... falling back to default."
			l2fc_cutoff = 1
	else:
		l2fc_cutoff = 1
		
	if '--p' in arguments:
		try:
			adj_pvalue_cutoff = float( arguments[ arguments.index( '--p' )+1 ] )
		except:
			print "ERROR: p input not recognized ... falling back to default."
			adj_pvalue_cutoff = 0.05
	else:
		adj_pvalue_cutoff = 0.05

	
	print "adjusted p-value cutoff: " + str( adj_pvalue_cutoff )
	print "log2() fold change cutoff: " + str( l2fc_cutoff )
	
	goi = load_genes_of_interest( input_file, l2fc_cutoff, adj_pvalue_cutoff )
	with open( output_file, "w" ) as out:
		out.write( "ID\tbaseMean\tl2fc\tpadj\tannotation\n" )
		for gene in sorted( goi, key=itemgetter('padj') ):
			try:
				out.write( "\t".join( map( str, [ gene['ID'], gene['baseMean'], gene['l2fc'], gene['padj'], anno[ gene['ID'] ] ] ) ) + '\n' )
			except KeyError:
				out.write( "\t".join( map( str, [ gene['ID'], gene['baseMean'], gene['l2fc'], gene['padj'] ] ) ) + '\tn/a\n' )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
