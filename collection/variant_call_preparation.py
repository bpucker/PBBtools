### Boas Pucker ###
### boas.pucker@rub.de ###
### v0.2 ###

### WARNING: DEPENDING ON THE CLUSTER ENVIRONMENT THE JAVA PATH NEEDS AN UPDATE ###

import os, sys, subprocess

# --- end of imports --- #

def index_bam_file( sorted_bam_file, piccard_tools, java_path ):
	"""! @brief index given BAM file """
	
	index_file = sorted_bam_file + ".bai"
	
	cmd = [ java_path,
			" -Xmx16g -XX:ParallelGCThreads=1 -jar ",
			piccard_tools,
			" BuildBamIndex I=",
			sorted_bam_file,
			" O=",
			index_file
		 ]
	p = subprocess.Popen( args="".join( cmd ), shell=True )
	p.communicate()
	return index_file


def add_read_groups( input_bam_file, output_bam_file, piccard_tools, java_path ):
	
	cmd = "".join( [ java_path,
					" -Xmx16g -XX:ParallelGCThreads=1 -jar ",
					piccard_tools,
					" AddOrReplaceReadGroups ",
					" I= ",
					input_bam_file,
					" O=",
					output_bam_file,
					" RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
					])
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()


def main( arguments ):
		
		# --- required inputs --- #
		original_bam_file = arguments[ arguments.index( '--original_bam_file' )+1 ]
		ref_file = arguments[ arguments.index( '--ref_file' )+1 ]
		bam_file = arguments[ arguments.index( '--bam_file' )+1 ]
		header = arguments[ arguments.index( '--header' )+1 ]
		samtools = arguments[ arguments.index( '--samtools' )+1 ]
		GATK = arguments[ arguments.index( '--GATK' )+1 ]
		piccard_tools = arguments[ arguments.index( '--piccard_tools' )+1 ]
		
		if '--java' in arguments:
			java_path = arguments[ arguments.index( '--java' )+1 ]
		else:
			java_path = "/vol/java-8/bin/java"
		
		# --- construction of bam file --- #
		#chr:start-end
		cmd = samtools + " view " + original_bam_file + " " + header + " -b > " + bam_file
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
		
		
		# --- indexing BAM file --- #
		index_bam_file( bam_file, piccard_tools, java_path )
		
		
		# --- adding read group --- #
		rg_bam_file = bam_file + "_RG_added.bam"
		add_read_groups( bam_file, rg_bam_file, piccard_tools, java_path )
		
		
		# --- indexing BAM file --- #
		index_bam_file( rg_bam_file, piccard_tools, java_path )


main( sys.argv )
