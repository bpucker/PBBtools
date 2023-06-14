### Boas Pucker ###
### b.pucker@tu-bs.de ###

__version__ = "v0.1"

__usage__ = """
					python BAM_splitter.py """ + __version__ + """
					--bamin <PATH_TO_INPUT_BAM_FILE>
					--ref <PATH_TO_REFERENCE_FASTA_FILE>
					--bamout <PATH_TO_BAM_OUTPUT_FILE>
					--header <HEADER>
					--samtools <PATH_TO_SAMTOOLS>
					--piccard_tools <PATH_TO_PICCARD_TOOLS>
					--java <PATH_TO_JAVA>
					"""

import os, sys

# --- end of imports --- #

def index_bam_file( sorted_bam_file, piccard_tools, java ):
	"""! @brief index given BAM file """
	
	index_file = sorted_bam_file + ".bai"
	
	cmd = [ java + " -Xmx8g -jar ",
			piccard_tools,
			" BuildBamIndex I=",
			sorted_bam_file,
			" O=",
			index_file
		 ]
	cmd = "".join( cmd )
	os.popen( cmd )
	return index_file


def add_read_groups( input_bam_file, output_bam_file, piccard_tools, java ):
	
	cmd = "".join( [ java + " -Xmx8g -jar ",
					piccard_tools,
					" AddOrReplaceReadGroups ",
					" I= ",
					input_bam_file,
					" O=",
					output_bam_file,
					" RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
					])
	os.popen( cmd )


def main( arguments ):
		
		# --- required inputs --- #
		if '--original_bam_file'  in arguments:
			original_bam_file = arguments[ arguments.index( '--original_bam_file' )+1 ]
		else:
			original_bam_file = arguments[ arguments.index( '--bamin' )+1 ]
		if '--ref_file' in arguments:
			ref_file = arguments[ arguments.index( '--ref_file' )+1 ]
		else:
			ref_file = arguments[ arguments.index( '--ref' )+1 ]
		if  '--bam_file' in arguments:
			bam_file = arguments[ arguments.index( '--bam_file' )+1 ]
		else:
			bam_file = arguments[ arguments.index( '--bamout' )+1 ]
		
		header = arguments[ arguments.index( '--header' )+1 ]
		samtools = arguments[ arguments.index( '--samtools' )+1 ]
		piccard_tools = arguments[ arguments.index( '--piccard_tools' )+1 ]
		if '--java' in arguments:
			java = arguments[ arguments.index( '--java' )+1 ]
		else:
			java = "/vol/java-8/bin/java"
		
		# --- construction of bam file --- #
		cmd = samtools + " view " + original_bam_file + " " + header + " -b > " + bam_file
		os.popen( cmd )
		
		
		# --- indexing BAM file --- #
		index_bam_file( bam_file, piccard_tool, java )
		
		
		# --- adding read group --- #
		rg_bam_file = bam_file + "_RG_added.bam"
		add_read_groups( bam_file, rg_bam_file, piccard_tools, java )
		
		
		# --- indexing BAM file --- #
		index_bam_file( rg_bam_file, piccard_tools, java )


if '--original_bam_file' in sys.argv and '--ref_file' in sys.argv and '--bam_file' in sys.argv and '--header' in sys.argv and '--samtools' in sys.argv and '--piccard_tools' in sys.argv:
	main( sys.argv )
elif '--bamin' in sys.argv and '--ref' in sys.argv and '--bamout' in sys.argv and '--header' in sys.argv and '--samtools' in sys.argv and '--piccard_tools' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
