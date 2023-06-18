### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python3 exp_plots.py
					--genes <GENES_FILE>
					--exp <EXPRESSION_FILE>
					--out <FIGURE_OUTPUT_FILE>
					
					optional:
					--cutfac <NUMBER_OF_IQRs_TO_DEFINE_OUTLIERS>
					--logscale <THIS_ACTIVATES_LOGSCALE>[off]
					--filteroff <SWITCHES_OUTLIER_FILTER_OFF>
					
					reference:  Apiaceae FNS I originated from F3H through tandem gene duplication. Boas Pucker, Massimo Iorizzo. bioRxiv 2022.02.16.480750; doi: https://doi.org/10.1101/2022.02.16.480750
					"""

import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
from scipy import stats
from pandas import DataFrame

# --- end of imports --- #


def load_genes( gene_file ):
	"""! @brief load genes """
	
	genes = {}
	gene_order = []
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			if "\t" in line:
				parts = line.strip().split('\t')
				genes.update( { parts[0]: parts[1] } )
				gene_order.append( parts[0] )
			else:
				genes.update( { line.strip(): line.strip() } )
				gene_order.append( line.strip() )
			line = f.readline()
	return genes, gene_order


def load_exp( exp_file ):
	"""! @brief load expression values """
	
	exp = {}
	with open( exp_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		for header in headers:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			values = list( map( float, parts[1:] ) )
			for idx, val in enumerate( values ):
				exp[ headers[ idx ] ].update( { parts[0]: val } )
			line = f.readline()
	return exp, headers


def generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_names, outlier_cutoff_factor, logscale, filter_status ):
	"""! @brief generate figure """
	
	universal_fontsize = 14
	sample_size_fontsize = 8	
	
	datamatrix = []
	max_val_per_gene = {}	#max value for each gene
	sample_size = {}	#number of valid values
	for idx, gene in enumerate( gene_order ):		
		# --- remove outliers --- #
		tmp = []
		for sample in sample_names:
			try:
				tmp.append( exp[ sample ][ gene ] )
			except KeyError:
				pass
		if len( tmp ) > 0:
			median = np.median( tmp )
			iqr = stats.iqr( tmp )
			valid = []
			for k, each in enumerate( tmp ):
				if filter_status:
					if abs( each - median ) < outlier_cutoff_factor * iqr:
						if logscale:
							datamatrix.append( [ "$\it{" + gene.replace( "-", "'").replace( "_", "\_" ) +"}$" , "x"+str( k ), math.log( each+1, 2) ] )
							valid.append( math.log( each+1, 2) )
						else:
							datamatrix.append( [ "$\it{" + gene.replace( "-", "'").replace( "_", "\_" ) +"}$" , "x"+str( k ), each ] )
							valid.append( each )
				else:
					if logscale:
						datamatrix.append( [ "$\it{" + gene.replace( "-", "'").replace( "_", "\_" ) +"}$" , "x"+str( k ), math.log( each+1, 2) ] )
						valid.append( math.log( each+1, 2) )
					else:
						datamatrix.append( [ "$\it{" + gene.replace( "-", "'").replace( "_", "\_" ) +"}$" , "x"+str( k ), each ] )
						valid.append( each )
	
		try:
			max_val_per_gene.update( { gene: max( valid ) } )
		except ValueError:
			max_val_per_gene.update( { gene: 0 } )
		sample_size.update( { gene: len( valid ) } )

	#gene, pigment_state, spec, value
	# datamatrix = [ 	[ "CHS", "x1", 10 ],
								# [ "CHS", "x2", 8 ],
								# [ "CHI", "x1", 5 ],
								# [ "CHI", "x2", 7 ]
							# ]
	try:
		print( sample_size )
	except:
		pass
	
	df = DataFrame( datamatrix, columns=[ "gene", "sample", "gene expression" ])
	
	fig, ax = plt.subplots( )
	x = sns.violinplot( 	x="gene",
									y="gene expression",
									data = df,	#DataFrame: column=variable, row=observation
									palette = [ "#87CEFA", "#FF8000" ],
									scale = "count",
									split = True,
									cut = 0
								)
	
	ax.set_xlabel( "" )
	ax.set_ylabel( "gene expression [TPM]" )
	
	plt.subplots_adjust( left=0.1, right=0.999, top=0.99, bottom=0.13, hspace=0.03  )
	
	fig.savefig( figfile, dpi=300 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief run generation of plots """
	
	gene_file = arguments[ arguments.index('--genes')+1 ]
	exp_file = arguments[ arguments.index('--exp')+1 ]
	figfile = arguments[ arguments.index('--out')+1 ]
	if '--cutfac' in arguments:
		try:
			outlier_cutoff_factor = float( arguments[ arguments.index('--cutfac')+1 ] )
		except:
			outlier_cutoff_factor = 3
	else:
		outlier_cutoff_factor = 3
	if "--logscale" in arguments:
		logscale = True
	else:
		logscale = False
	
	if '--filteroff' in arguments:
		filter_status = False
	else:
		filter_status = True
	
	genes, gene_order = load_genes( gene_file )

	exp, sample_names = load_exp( exp_file )

	generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_names, outlier_cutoff_factor, logscale, filter_status )


if '--genes' in sys.argv and '--exp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
