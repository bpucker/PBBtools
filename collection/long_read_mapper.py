### Boas Pucker ###
### boas.pucker@rub.de ###
### v0.15 ###

__usage__ = """
					python long_read_mapper.py
					--reads <READ_FILE>
					--ref <REFERENCE_SEQUENCE_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--tech <ONT|PB>[ONT]
					--cov (activates generation of COVERAGE file)
					--cpus <INT, number of threads for mapping and sorting>[20]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #


def main( arguments ):
	"""! @brief run everthing """
	
	read_file = arguments[ arguments.index('--reads')+1 ]
	ref_file = arguments[ arguments.index('--ref')+1 ]
	output_directory = arguments[ arguments.index('--out')+1 ]
	
	if '--tech' in arguments:
		mode = arguments[ arguments.index('--tech')+1 ]
	else:
		mode = "ONT"
	
	if '--cpus' in arguments:
		cpus = int( arguments[ arguments.index('--cpus')+1 ] )
	else:
		cpus = 20
	
	
	minimap2_path = "minimap2"
	samtools_path = "samtools"
	cov_file_script = "/grp/pbb/scripts/construct_cov_file.py"
	
	if output_directory[-1] != '/':
		output_directory += "/"
	
	if not os.path.exists( output_directory ):
		os.makedirs( output_directory )
	
	# --- mapping step ---- #
	mapping_file = output_directory + "mapping.sam"
	if mode == "ONT":
		cmd = " ".join( [ minimap2_path, "-ax map-ont --secondary=no -t", str( cpus ), ref_file, read_file, ">", mapping_file ] )
	else:
		cmd = " ".join( [ minimap2_path, "-ax map-pb --secondary=no -t", str( cpus ), ref_file, read_file, ">", mapping_file ] )
	os.popen( cmd )
	
	# --- converting SAM to BAM --- #
	mapping_bam_file = output_directory + "mapping.bam"
	cmd = " ".join( [ samtools_path, "view -Sb", mapping_file, ">", mapping_bam_file ] )
	os.popen( cmd )
	
	# --- sort BAM file --- #
	sorted_mapping_bam_file = output_directory + "mapping_sorted.bam"
	cmd = " ".join( [ samtools_path, "sort -@ "+ str( cpus ), "-o " + sorted_mapping_bam_file, mapping_bam_file ] )
	os.popen( cmd )
	
	# --- generation of coverage (COV) file --- #
	if '--cov' in arguments:
		cov_file = output_directory + "mapping_sorted.cov"
		cmd = " ".join( [ "python", cov_file_script, "--in", sorted_mapping_bam_file, "--out", cov_file, "--bam_is_sorted" ] )
		os.popen( cmd )


if '--reads' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
