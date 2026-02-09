### Boas Pucker ###
### pucker@uni-bonn.de ###
### v0.31 ###

__usage__ = """
					python3 heatmap_plotter.py
					--genes <GENES_FILE>
					--exp <EXPRESSION_FILE>
					--out <FIGURE_OUTPUT_FILE>
					
					optional:
					--norm (activates normalization per gene)
					--samples <SAMPLE_METADATA>
					"""

import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import pandas as pd

# --- end of imports --- #


def load_genes( gene_file ):
	"""! @brief load genes """
	
	genes = {}
	gene_order = []
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip()
			if "\t" in parts:
				parts = parts.split('\t')
				genes.update( { parts[0]: parts[1] } )
				gene_order.append( parts[0] )
			else:
				genes.update( { parts: parts } )
				gene_order.append( parts )
			line = f.readline()
	return genes, gene_order


def load_exp( gene_file ):
	"""! @brief load expression values """
	
	exp = {}
	with open( gene_file, "r" ) as f:
		sample_order = f.readline().strip().split('\t')
		if sample_order[0] in [ "genes", "gene", "GeneID", "GeneIDs" ]:
			sample_order = sample_order[1:]
		for header in sample_order:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			values = map( float, parts[1:] )
			for idx, val in enumerate( values ):
				exp[ sample_order[ idx ] ].update( { parts[0]: val } )
			line = f.readline()
	return exp, sample_order


def generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_order, norm, samples ):
	"""! @brief generate figure """
	
	exp_dataset = []	#[ [ column1 ], [column2], [column3] ]
	gene_lables = []
	for gene in gene_order:
		values = []
		for sample in sample_order:
			values.append( exp[ sample ][ gene ] )
		exp_dataset.append( values )
		gene_lables.append( genes[ gene ] )
	sample_name_display = []
	for each in sample_order:
		try:
			sample_name_display.append( samples[ each ] )
		except KeyError:
			sample_name_display.append( each )
	data = pd.DataFrame( exp_dataset, columns=sample_name_display )
	data.index = gene_lables + []
	
	if norm:
		data = data.apply( lambda x: (x-x.mean())/x.std(), axis = 1)
	
	fig, ax = plt.subplots()
	
	h = sns.heatmap( data, annot=True, fmt=".1f", cmap="crest", annot_kws={"size": 3}, xticklabels = 1, yticklabels=1 )
	
	ax.set(xlabel="", ylabel="")
	ax.xaxis.tick_top()
	h.set_yticklabels( h.get_yticklabels(), rotation=0, size = 6 )	#turn y-axis labels to horizontal
	h.set_xticklabels( h.get_xticklabels(), rotation=90, size = 6 )	#turn y-axis labels to horizontal
	
	fig.tight_layout()
	
	#plt.subplots_adjust( left=0.2, right=0.999, top=0.9, bottom=0.01, hspace=0.03  )
	
	fig.savefig( figfile, dpi=300 )
	plt.close( "all" )


def load_sample_names( sample_file ):
	"""! @brief load sample names from input file """
	
	samples = {}
	sample_order = []
	with open( sample_file, "r" ) as f:
		line = ''.join([ i if ord(i) < 128 else ' ' for i in f.readline() ])
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 1:
				samples.update( { parts[0]: parts[1] } )
				sample_order.append( parts[0] )
			line = ''.join([ i if ord(i) < 128 else ' ' for i in f.readline() ])
	return samples, sample_order


def main( arguments ):
	"""! @brief run generation of plots """
	
	gene_file = arguments[ arguments.index('--genes')+1 ]
	exp_file = arguments[ arguments.index('--exp')+1 ]
	figfile = arguments[ arguments.index('--out')+1 ]
	
	if '--norm' in arguments:
		norm = True
	else:
		norm = False
	
	genes, gene_order = load_genes( gene_file )

	exp, sample_order = load_exp( exp_file )
	
	if '--samples' in arguments:
		sample_file = arguments[ arguments.index('--samples')+1 ]
		samples, sample_order = load_sample_names( sample_file )
	else:
		samples = {}

	generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_order, norm, samples )


if '--genes' in sys.argv and '--exp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
