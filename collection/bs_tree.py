### Boas Pucker ###
### b.pucker@tu-bs.de ###
__version__ ="v0.21"

__usage__ = """
	python bs_tree.py
	--pep <FULL_PATH_TO_PEPTIDE_INPUT_FILE>
	--cds <FULL_PATH_TO_CDS_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIR>
	--taxon <TAXON_FILE>
	
	
	optional:
	--tree <fast|raxml>[raxml]
	--occ <FLOAT, occupancy required per alignment column>[0.1]
	--name <STRING, prefix for final alignment file>
	--mafft <FULL_PATH_TO_MAFFT>
	--pxaa2cdn <PATH_TO_pxaa2cdn>
	--pxclsq <PATH_TO_pxclsq>
	--fasttree <PATH_TO_FastTree>
	--raxml <PATH_TO_RAXML-NG>
	--renametree <PATH_TO_RENAME_TREE_SEQS.PY>
	--renmaefasta <PATH_TO_RENAME_FASTA_SEQS.PY>
					"""

import os, sys

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().replace( " ", "" ).replace('(', '_').replace(')', '_').replace( "	", "" ).replace('.', '_').replace(':', '_').replace(';', '_')
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
				if len( seq ) > 0:
					sequences.update( { header: seq } )
				header = line.strip()[1:].replace( " ", "" ).replace('(', '_').replace(')', '_').replace( "	", "" ).replace('.', '_').replace(':', '_').replace(';', '_')
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		if len( seq ) > 0:
			sequences.update( { header: seq } )
	return sequences


def main( arguments ):
	"""! @brief handle everything """
	
	pep_input_file = arguments[ arguments.index('--pep') + 1 ]
	cds_input_file = arguments[ arguments.index('--cds') + 1 ]
	output_dir = arguments[ arguments.index('--out') + 1 ]
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--mafft') + 1 ]
	else:
		mafft = "mafft"
	
	if '--pxaa2cdn' in arguments:
		pxaa2cdn = arguments[ arguments.index('--pxaa2cdn') + 1 ]
	else:
		pxaa2cdn = "pxaa2cdn"
	
	if '--pxclsq' in arguments:
		pxclsq = arguments[ arguments.index('--pxclsq') + 1 ]
	else:
		pxclsq = "pxclsq"
	
	if '--fasttree' in arguments:
		fasttree = arguments[ arguments.index('--fasttree') + 1 ]
	else:
		fasttree = "/grp/pbb/tools/FastTree"
	
	if '--raxml' in arguments:
		raxml = arguments[ arguments.index('--raxml') + 1 ]
	else:
		raxml = "/grp/pbb/tools/RAxML9/raxml-ng"
	
	if '--renametree' in arguments:
		renametree = arguments[ arguments.index('--renametree') + 1 ]
	else:
		renametree = "/grp/pbb/scripts/rename_seqs_in_trees.py"
	
	if '--renamefasta' in arguments:
		renamefasta = arguments[ arguments.index('--renamefasta') + 1 ]
	else:
		renamefasta = "/grp/pbb/scripts/rename_FASTA_seqs.py"
	
	if '--occ' in arguments:
		occupancy = float( arguments[ arguments.index('--occ') + 1 ] )
	else:
		occupancy = 0.1
	
	if '--name' in arguments:
		name = arguments[ arguments.index('--name') + 1 ]
	else:
		name = "XXX"
	
	if "--tree" in arguments:
		tree_tool = arguments[ arguments.index('--tree') + 1 ]
		if tree_tool not in [ "fast", "raxml" ]:
			tree_tool = "fast"
	else:
		tree_tool = "raxml"
		
	
	if output_dir[-1] != '/':
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	# --- modify FASTA files to ensure proper and compatible names --- #
	mod_pep_FASTA = output_dir + "names_modified.pep.fasta"
	mod_cds_FASTA = output_dir + "names_modified.cds.fasta"
	peps = load_sequences( pep_input_file )
	cdss = load_sequences( cds_input_file )
	
	with open( mod_pep_FASTA, "w" ) as out:
		for key in peps.keys():
			out.write( '>' + key + "\n" + peps[ key ] + "\n" )
	with open( mod_cds_FASTA, "w" ) as out:
		for key in cdss.keys():
			out.write( '>' + key + "\n" + cdss[ key ] + "\n" )
	
	# --- construct CDS alignment --- #
	pep_alignment_file = mod_pep_FASTA + ".aln"
	os.popen( " ".join( [ mafft, mod_pep_FASTA, ">", pep_alignment_file ] ) )
	
	# --- convert pep to cds --- #
	cds_alignment_file = mod_cds_FASTA + ".aln"
	os.popen( " ".join( [ pxaa2cdn, "-a", pep_alignment_file, "-n", mod_cds_FASTA, "-o", cds_alignment_file ] ) )
	
	# --- alignment cleaning --- #
	clean_cds_alignment_file = cds_alignment_file + ".cln"
	os.popen( " ".join( [ pxclsq, "-s", cds_alignment_file, "-o", clean_cds_alignment_file, "-p", str( occupancy ) ] ) )
	
	
	if tree_tool == "fast":
		# --- construct fasttree tree --- #
		tree_file = clean_cds_alignment_file + ".tre"
		os.popen( " ".join( [ fasttree, "-gtr -nt <", clean_cds_alignment_file, ">", tree_file ] ) )
	elif tree_tool == "raxml":
		# --- contructing RAxML-NG tree --- #
		prefix = clean_cds_alignment_file.split('.cds.fa')[0]
		tree_file = prefix + ".raxml.bestTree"
		os.popen( " ".join( [ raxml, "--all --threads 1 --model GTR+F+G --bs-trees 200 --msa", clean_cds_alignment_file, "--prefix", prefix ] ) )
		#Felsenstein slow bootstrap
		
	# --- rename sequences in tree and FASTA files --- #
	if '--taxon' in arguments:
		taxon_file = arguments[ arguments.index('--taxon') + 1 ]
		output_tree_file = output_dir + name + ".tre"
		os.popen( " ".join( [ "python", renametree, "--in", tree_file, "--out", output_tree_file, "--taxon", taxon_file ] ) )
		
		renamed_aln_file = output_dir + name + ".aln.cln"
		os.popen( " ".join( [ "python", renamefasta, "--in", clean_cds_alignment_file, "--out", renamed_aln_file, "--taxon", taxon_file ] ) )


if '--pep' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
