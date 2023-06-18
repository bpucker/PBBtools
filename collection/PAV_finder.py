### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.41 ###

__usage__ = """
					python3 PAV_finder3.py
					--cov1 <COVERAGE_FILE1>
					--cov2 <COVERAGE_FILE2>
					--out <OUTPUT_FOLDER>
					
					optional:
					--mode <MODE_OF_PAV_DETECTION(gene|genomic|zcr)>[genomic]
					--gff <GFF_ANNOTATION_FILE>
					--anno <FUNCTIONAL_ANNOTATION_FILE>
					--mincov <MINIMAL_COMBINED_COVERAGE_OF_BOTH_SAMPLES_PER_GENE>
					--blocksize <SIZE_FOR_GENOMIC_PAV_OR_ZCR_DETECTION>[3000]
					--maxrelcov <RELATIVE_COVERAGE_CUTOFF_OF_ZCR_DETECTION>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, glob, sys, gzip, re, math
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

# --- end of imports --- #

def get_mean( values ):
	"""! @brief calculate mean of given list of values """
	
	if values != []:
		return sum( values ) / float( len( values ) )
	else:
		return 0


def load_coverage( cov_file ):
	"""! @brief load coverage from given file """
	
	if cov_file[-3:].lower() == "cov" or cov_file[-3:].lower() == "txt":	# uncompressed coverage file
		coverage = {}
		with open( cov_file, "r" ) as f:
			line = f.readline()
			prev_chr = line.split('\t')[0]
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0]
					cov = []
				try:
					cov.append( float( parts[2] ) )
				except IndexError:
					print( line )
				line = f.readline()
			coverage.update( { prev_chr: cov } )
		return coverage
	
	else:	# compressed coverage file
		coverage = {}
		with gzip.open( cov_file, "rb" ) as f:
			line = f.readline()
			prev_chr = line.split('\t')[0]
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0]
					cov = []
				cov.append( float( parts[2] ) )
				line = f.readline()
			coverage.update( { prev_chr: cov } )
		return coverage


def load_genes( gff_file ):
	"""! @brief load gene positions from GFF3 file """
	
	genes = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "mRNA":	
					if ";" in parts[-1]:
						ID = parts[-1].split(';')[0].split('=')[1]
					else:
						ID = parts[-1].split('=')[1]
					genes.update( { ID: { 'chr': parts[0], 'start': int(parts[3]), 'end': int( parts[4] ) } } )
			line = f.readline()
	return genes


def generate_genomic_blocks( cov, blocksize ):
	"""! @brief get genomic blocks to replace genes """
	
	genes = {}
	for chromosome in cov.keys():
		start = 0
		end = 0 + blocksize
		while end < len( cov[ chromosome ] ):
			if min( [ end, len( cov[ chromosome ] ) ] ) - start >= blocksize:
				genes.update( { 	chromosome + "_%_" + str( start ) + "_%_" + str( min( [ end, len( cov[ chromosome ] ) ] ) ):
											{ 'chr': chromosome, 'start': start, 'end': min( [ end, len( cov[ chromosome ] ) ] ) } 
										} )
			start += blocksize
			end += blocksize
	return genes


def calculate_cov_per_gene( cov, genes ):
	"""! @brief calculate coverage per gene """
	
	cov_per_gene = {}
	for key in genes.keys():
		values = cov[ genes[ key ]['chr'] ][ genes[ key ]['start']: genes[ key ]['end'] ]
		cov_per_gene.update( { key: get_mean( values ) } )
	return cov_per_gene


def generate_figure( values_to_plot, fig_file ):
	"""! @brief generate figure with coverage values per gene """
	
	max_value = 5 * get_mean( values_to_plot )
	
	fig, ax = plt.subplots()
	
	ax.hist( values_to_plot, bins=int( max_value ), range=( 0, max_value ) )
	ax.set_xlabel( "average coverage per gene" )
	ax.set_ylabel( "number of genes" )
	
	ax.set_xlim( 0, max_value )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	fig.savefig( fig_file, dpi=600 )
	plt.close( "all" )


def load_cov_per_gene( doc_file ):
	"""! @brief load coverage per gene """
	
	coverages = {}
	with open( doc_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			coverages.update( { parts[0]: float( parts[1] ) } )
			line = f.readline()
	return coverages


def PAV_detection( cov_per_gene1, cov_per_gene2, PAV_file, anno, cov_cutoff ):
	"""! @brief detection of PAVs """
	
	avg_cov1 = float( np.median( list( cov_per_gene1.values() ) ) )
	avg_cov2 = float( np.median( list( cov_per_gene2.values() ) ) )
	
	# --- calculate ratios --- #
	
	ratios = {}
	norm_cov_per_gene1 = {}
	norm_cov_per_gene2 = {}
	for key in list( cov_per_gene1.keys() ):
		cov1 = cov_per_gene1[ key ]
		cov2 = cov_per_gene2[ key ]
		if cov1+cov2 > cov_cutoff:
			if cov2 == 0:
				if cov1 == 0:
					ratios.update( { key: 0 } )
				else:
					cov1 = cov1 / avg_cov1
					cov2 = 0.01 / avg_cov2
					ratios.update( { key: np.log2( cov1 / cov2 ) } )	
					norm_cov_per_gene1.update( { key: cov1 } )
					norm_cov_per_gene2.update( { key: cov2 } )
			else:
				if cov1 == 0:
					cov1 = 0.01 / avg_cov1
					cov2 = cov2 / avg_cov2
				else:
					cov1 = cov1 / avg_cov1
					cov2 = cov2 / avg_cov2
				ratios.update( { key: np.log2( cov1 / cov2 ) } )
				norm_cov_per_gene1.update( { key: cov1 } )
				norm_cov_per_gene2.update( { key: cov2 } )
	
	mean = np.mean( list( ratios.values() ) )
	sd = np.std( list( ratios.values() ) )
	
	# --- calculate modified z-scores --- #
	
	mad = []
	ratio_med = np.median( list( ratios.values() ) )
	for gene in list( ratios.keys() ):
		if ratios[ gene ] > ratio_med:
			mad.append( ratios[ gene ] - ratio_med )
		else:
			mad.append( ratio_med - ratios[ gene ] )
	mad = np.median( mad )
	mod_z_scores = {}
	for gene in list( ratios.keys() ):
		mod_z_scores.update( { gene: ( 0.6745 * ( ratios[ gene ] - ratio_med ) ) / mad } )
	
	# --- sort values and generate output file --- #
	
	data_for_sorting = []
	for key in list( ratios.keys() ):
		data_for_sorting.append( { 'id': key, 'ratio': ratios[ key ], 'z': mod_z_scores[ key ], 'abs_z': abs( mod_z_scores[ key ] ) } )
	data_for_sorting = sorted( data_for_sorting, key=itemgetter('abs_z') )[::-1]
	
	with open( PAV_file, "w" ) as out:
		out.write( "GeneID\tCovS1\tCovS2\tlog2(NormCov1/NormCov2)\tMod_Z-score\tAnnotation\n" )
		for entry in data_for_sorting:
			try:
				try:
					out.write( "\t".join( map( str, [ 	entry['id'],
																		round( norm_cov_per_gene1[ entry['id'] ], 3 ),
																		round( norm_cov_per_gene2[ entry['id'] ], 3 ),
																		round( entry['ratio'], 3 ),
																		round( entry['z'], 3 ),
																		anno[ entry['id'] 
																	] ] ) ) + '\n' )
				except KeyError:
					out.write( "\t".join( map( str, [ 	entry['id'],
																		round( norm_cov_per_gene1[ entry['id'] ], 3 ),
																		round( norm_cov_per_gene2[ entry['id'] ], 3 ),
																		round( entry['ratio'], 3 ),
																		round( entry['z'], 3 ),
																		"n/a" 
																	] ) ) + '\n' )
			except KeyError:
				pass


def ZCR_detection( cov_per_gene1, cov_per_gene2, ZCR_file, cov_cutoff, maxrelcov ):
	"""! @brief detection of PAVs """
	
	avg_cov1 = float( np.median( list( cov_per_gene1.values() ) ) )
	avg_cov2 = float( np.median( list( cov_per_gene2.values() ) ) )
	
	# --- calculate ratios --- #
	
	ZCRs = []
	for key in list( sorted( list( cov_per_gene1.keys() ) ) ):
		cov1 = cov_per_gene1[ key ]
		cov2 = cov_per_gene2[ key ]
		if cov1<= cov_cutoff and cov2 <= cov_cutoff:
			norm_cov1 = cov1 / avg_cov1
			norm_cov2 = cov2 / avg_cov2
			if norm_cov1 <= maxrelcov and norm_cov2 <= maxrelcov:
				ZCRs.append( key )
	with open( ZCR_file, "w" ) as out:
		out.write( "\n".join( ZCRs ) + '\n' )
	return ZCRs


def load_annotation( anno_file ):
	"""! @brief load annotation from given file """
	
	anno = {}
	
	with open( anno_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno.update( { parts[0]: parts[1] } )
			line = f.readline()
	
	return anno


def main( arguments ):
	"""! @brief runs everything """
	
	cov_file1 = arguments[ arguments.index( '--cov1' )+1 ]
	cov_file2 = arguments[ arguments.index( '--cov2' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	
	if '--mode' in arguments:
		mode = arguments[ arguments.index( '--mode' )+1 ].lower()
	else:
		mode = "genomic"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index( '--anno' )+1 ]
		anno = load_annotation( anno_file )
	else:
		anno = {}
	
	if '--mincov' in arguments:
		cov_cutoff = float( arguments[ arguments.index( '--mincov' )+1 ] )
	else:
		cov_cutoff = -1
	
	if '--blocksize' in arguments:
		blocksize = int( arguments[ arguments.index( '--blocksize' )+1 ] )
	else:
		blocksize = 3000
	
	if '--maxrelcov' in arguments:
		maxrelcov = float( arguments[ arguments.index( '--maxrelcov' )+1 ] )
	else:
		maxrelcov = 0.1
	
	# --- stop execution if results are impossible --- #
	
	if mode == "zcr":
		if cov_cutoff < 0:
			sys.exit( "ERROR: coverage cannot be a negative value (--mincov needs to be 0 or higher)" )
		if maxrelcov < 0:
			sys.exit( "ERROR: coverage cannot be a negative value (--maxrelcov needs to be 0 or higher)" )
	
	
	doc_file1 = output_dir + "cov_per_gene1.txt"
	doc_file2 = output_dir + "cov_per_gene2.txt"
	
	# --- detection of PAVs --- #
	
	if mode in [ "gene", "genomic" ]:
		if mode == "gene":
			gff_file = arguments[ arguments.index( '--gff' )+1 ]
			genes = load_genes( gff_file )
			print( "number of genes to check: " + str( len( list( genes.keys() ) ) ) )
		
		
		# --- handle coverage 1 --- #
		
		genomic_block_state = False
		if not os.path.isfile( doc_file1 ):
			cov = load_coverage( cov_file1 )
			if mode == "genomic":
				genes = generate_genomic_blocks( cov, blocksize )
				print( "number of genomic blocks to check: " + str( len( list( genes.keys() ) ) ) )
				genomic_block_state = True
			cov_per_gene1 = calculate_cov_per_gene( cov, genes )
			with open( doc_file1, "w" ) as out:
				for gene in sorted( list( cov_per_gene1.keys() ) ):
					out.write( gene + '\t' + str( cov_per_gene1[gene] ) + '\n' )
		else:
			cov_per_gene1 = load_cov_per_gene( doc_file1 )
		
		# --- handle coverage 2 --- #
		
		if not os.path.isfile( doc_file2 ):
			cov = load_coverage( cov_file2 )
			if mode == "genomic" and not genomic_block_state:
				genes = generate_genomic_blocks( cov, blocksize )
				print( "number of genomic blocks to check: " + str( len( list( genes.keys() ) ) ) )
			cov_per_gene2 = calculate_cov_per_gene( cov, genes )
			with open( doc_file2, "w" ) as out:
				for gene in sorted( cov_per_gene2.keys() ):
					out.write( gene + '\t' + str( cov_per_gene2[gene] ) + '\n' )
		else:
			cov_per_gene2 = load_cov_per_gene( doc_file2 )
	
	# --- detection of ZCRs --- #
	
	elif mode == "zcr":
		# --- handle coverage 1 --- #
		genomic_block_state = False
		if not os.path.isfile( doc_file1 ):
			cov = load_coverage( cov_file1 )
			genes = generate_genomic_blocks( cov, blocksize )
			print( "number of genoimc blocks to check: " + str( len( list( genes.keys() ) ) ) )
			genomic_block_state = True
			cov_per_gene1 = calculate_cov_per_gene( cov, genes )
			with open( doc_file1, "w" ) as out:
				for gene in list( sorted( list( cov_per_gene1.keys() ) ) ):
					out.write( gene + '\t' + str( cov_per_gene1[gene] ) + '\n' )
		else:
			cov_per_gene1 = load_cov_per_gene( doc_file1 )
		
		# --- handle coverage 2 --- #
		if not os.path.isfile( doc_file2 ):
			cov = load_coverage( cov_file2 )
			if not genomic_block_state:
				genes = generate_genomic_blocks( cov, blocksize )
				print( "number of genoimc blocks to check: " + str( len( list( genes.keys() ) ) ) )
			cov_per_gene2 = calculate_cov_per_gene( cov, genes )
			with open( doc_file2, "w" ) as out:
				for gene in list( sorted( list( cov_per_gene2.keys() ) ) ):
					out.write( gene + '\t' + str( cov_per_gene2[gene] ) + '\n' )
		else:
			cov_per_gene2 = load_cov_per_gene( doc_file2 )
	
	
	fig_file1 = output_dir + "coverage_per_gene1.pdf"
	fig_file2 = output_dir + "coverage_per_gene2.pdf"
	generate_figure( cov_per_gene1.values(), fig_file1 )
	generate_figure( cov_per_gene2.values(), fig_file2 )
	
	if mode in [ "gene", "genomic" ]:
		PAV_file = output_dir + "PAVs.txt"
		PAV_detection( cov_per_gene1, cov_per_gene2, PAV_file, anno, cov_cutoff )
	elif mode == "zcr":
		ZCR_file = output_dir + "ZCRs.txt"
		
		ZCRs = ZCR_detection( cov_per_gene1, cov_per_gene2, ZCR_file, cov_cutoff, maxrelcov )
		print( "number of detected ZCRs: " + str( len( ZCRs ) ) )


if '--cov1' in sys.argv and '--cov2' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
