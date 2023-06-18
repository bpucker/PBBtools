### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

### some functions taken from MYB_annotator.py and KIPEs ###

__version__ = "v0.1"

__usage__ = """
					python3 region_check.py
					--config <CONFIG_FILE>
					--out <OUTPUT_FOLDER>
					"""

import os, re, sys, subprocess, math, glob
from operator import itemgetter

# --- end of imports --- #

def load_config_file( config_file ):
	"""! @brief load config file """
	info = []
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				info.append( { 'id': parts[0], 'gff': parts[1], 'cds': parts[2], 'spec': parts[3], 'feature_type': parts[4], 'id_type': parts[5] } )
			line = f.readline()
	return info


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
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_flanking_genes( gff_file, target_gene_ID, feature_type, id_type, num_to_take ):
	"""! @brief load flanking genes from given GFF3 file """
	
	genes_per_seq = {}
	target_seq = False
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == feature_type:	#gene row is reached
					gene_ID = parts[-1].split( id_type+"=" )[1]
					if ";" in gene_ID:
						gene_ID = gene_ID.split(';')[0]
						if gene_ID == target_gene_ID:
							target_seq = parts[0]
						try:
							genes_per_seq[ parts[0] ].append( { 'id': gene_ID, 'start': int( parts[3] ) } )
						except KeyError:
							genes_per_seq.update( { parts[0]: [ { 'id': gene_ID, 'start': int( parts[3] ) } ] } )
			line = f.readline()
	if target_seq:
		genes_on_target_seq = list( sorted( genes_per_seq[ target_seq ], key=itemgetter("start") ) )
		target_position = False
		for idx, each in enumerate( genes_on_target_seq ):
			if each['id'] == target_gene_ID:
				target_position = idx + 0
		selection = genes_on_target_seq[ max( [ 0, target_position-num_to_take ] ): min( [ target_position+num_to_take, len( genes_on_target_seq )-1 ] ) ]
		selected_IDs = []
		for each in selection:
			selected_IDs.append( each['id'] )
		return selected_IDs
	else:
		return []


def load_best_blast_hit_per_query( blast_result_file ):
	"""! @brief load best BLAST hit per query sequence """
	
	best_her_per_query = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[-1] ) > best_her_per_query[ parts[0] ]['score']:
					best_her_per_query[ parts[0] ] = { 'subject': parts[1], 'score': float( parts[-1] ) }
			except KeyError:
				best_her_per_query.update( { parts[0]: { 'subject': parts[1], 'score': float( parts[-1] ) } } )
			line = f.readline()
	
	best_IDs = []
	for each in list( best_her_per_query.values() ):
		best_IDs.append( each['subject'] )
	return list( set( best_IDs ) )


def main( arguments ):
	"""! @brief run everything """
	
	config_file = arguments[ arguments.index('--config')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	num_to_take = 100
	cpu = 4
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	config_data = load_config_file( config_file )
	summary_across_species = {}
	flanking_genes_per_spec = {}
	all_species = []
	for entry in config_data:
		cds_seqs = load_sequences( entry['cds'] )
		all_species.append( entry['spec'] )
		
		# --- identify flanking genes --- #
		flanking_genes = load_flanking_genes( entry['gff'], entry['id'], entry['feature_type'], entry['id_type'], num_to_take )	#load feature and ID type from config file
		print ( entry['spec'] + ": " + str( len( flanking_genes ) ) )
		flanking_genes_per_spec.update( { entry['spec']:  flanking_genes} )
		
		# --- construct bait sequence file --- #
		bait_file = output_folder + entry['spec'] + ".baits.fasta"
		if not os.path.isfile( bait_file ):	#construct bait sequence file only if it does not exist
			with open( bait_file, "w" ) as out:
				for gene in flanking_genes:
					out.write( '>' + gene + "\n" + cds_seqs[ gene ] + "\n" )
		
		# --- run BLAST search vs. all others --- #
		result_summary = {}
		for entry2 in config_data:
			if entry2['spec'] != entry['spec']:	#do not run BLAST search vs. self
				blast_result_file = output_folder + entry['spec'] + "_vs_" + entry2['spec'] + ".results.txt"
				if not os.path.isfile( blast_result_file ):	#run BLAST search only if result files are not present yet
					blast_db = output_folder + entry2['spec'] + "_blastdb"
					p = subprocess.Popen( args= "makeblastdb -in " + entry2[ 'cds' ] + " -out " + blast_db + " -dbtype nucl", shell=True )
					p.communicate()
					p = subprocess.Popen( args= "tblastx -query " + bait_file+ " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( cpu ), shell=True )
					p.communicate()
				best_blast_hit_per_query = load_best_blast_hit_per_query( blast_result_file )	#returns only a list of the best hits
				result_summary.update( { entry2['spec']: best_blast_hit_per_query } )
		summary_across_species.update( { entry['spec']:  result_summary} )
	
	# --- compare best hits between species to find matches --- #
	#summary_across_species = { 'spec1': { 'spec2': [ g1, g2, g3 ], 'spec3': [ g1, g2, g3], ... },  'spec2': { 'spec1': [ g1, g2, g3 ], 'spec3': [ g1, g2, g3], ... }, ... }
	for idx, spec in enumerate( all_species ):
		for idx2, spec2 in enumerate( all_species ):
			if idx2 > idx:
				self_genes = flanking_genes_per_spec[ spec ]
				hit_self_genes = summary_across_species[ spec2 ][ spec ]
				matches = []
				for gene in self_genes:
					if gene in hit_self_genes:
						matches.append( gene )
				print ( spec + " vs. " + spec2 + ": " + str( len( matches ) ) )
				if len( matches ) > 0:
					print ( matches )


if '--config' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
