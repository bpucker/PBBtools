# Boas Pucker #
# b.pucker@tu-bs.de #
# v0.21 #


import os, sys, re, glob, datetime, subprocess, time
from operator import itemgetter

# --- end of imports --- #

__usage__ = """ python GATK4_wrapper_slurm.py
							--in <PATH_TO_BAM_FILE>
							--ref <PATH_TO_REFERENCE_FILE>
							--out <PATH_TO_DIRECTORY>
							--user <USER_ID>
							
							optional:
							--bam_is_sorted (prevents sorting of bam file)
							--block <BLOCK_SIZE_IN_BP>[5Mbp]
							--para_jobs <NUMBER_OF_JOBS_TO_RUN_IN_PARALLEL_ON_CLUSTER>[100]
							--script <PATH_TO_AUXILIARY_SCRIPT>[variant_call_preparation.py]
							--samtools <SAMTOOLS_PATH>[samtools]
							--picard <PICARD_PATH>[picard.jar]
							--gatk <PATH_TO_GATK>[gatk-package-4.1.9.0-local.jar]
							
							bug reports and feature requests: b.pucker@tu-bs.de
						"""


def load_sequences( filename ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( filename ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def get_seq_lengths( fasta ):
	"""! brief get lengths of all sequences """
	
	seqs = load_sequences( fasta )
	lengths = {}
	for key in seqs.keys():
		lengths.update( { key: len( seqs[ key ] ) } )
	return lengths


def sort_bam_file( prefix, input_bam_file, samtools_path, piccard_tools ):
	"""! @brief sort given bam file by position """
	
	# --- sort bam file --- #
	
	tmp_file = prefix + input_bam_file.split('/')[-1] + "_tmp.bam"
	
	cmd = samtools_path + " sort -@ 8 -m 16G -o " +  tmp_file + " " + input_bam_file
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()
	
	# --- adding read group IDs --- #
	output_file = prefix + input_bam_file.split('/')[-1] + "_sorted.bam"
	cmd = [ "java -Xmx32g -jar ",
			piccard_tools,
			" AddOrReplaceReadGroups I=",
			tmp_file,
			" O=",
			output_file,
			" RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
		]
	p = subprocess.Popen( args= "".join( cmd ), shell=True )
	p.communicate()
	
	return output_file


def index_bam_file( sorted_bam_file, piccard_tools ):
	"""! @brief index given BAM file """
	
	index_file = sorted_bam_file + ".bai"
	
	cmd = [ "java -Xmx64g -XX:ParallelGCThreads=4 -jar ",
			piccard_tools,
			" BuildBamIndex I=",
			sorted_bam_file,
			" O=",
			index_file
		 ]
	cmd = "".join( cmd )
	p = subprocess.Popen( args= "".join( cmd ), shell=True )
	p.communicate()
	return index_file


def index_ref_file( ref_fasta_file, piccard_tools, samtools ):
	"""! @brief construct indexed reference file """
	
	dictionary_file = '.'.join( ref_fasta_file.split('.')[:-1] ) + ".dict"
	cmd1 = [ 	"java -Xmx64g -XX:ParallelGCThreads=4 -jar ",
				piccard_tools,
				" CreateSequenceDictionary R=",
				ref_fasta_file,
				" O=",
				dictionary_file
			]
	p = subprocess.Popen( args="".join( cmd1 ), shell=True )
	p.communicate()
	time.sleep(30)
	sys.stdout.write( "reference dictionary file created\n" )
	sys.stdout.flush()
	
	cmd2 = [ 	samtools,
				" faidx ",
				ref_fasta_file,
			]	
	p = subprocess.Popen( args="".join( cmd2 ), shell=True )
	p.communicate()
	time.sleep(30)
	sys.stdout.write( "reference index file created\n" )
	sys.stdout.flush()


def prepare_files_for_variant_calling_on_cluster( prefix, input_bam_file, headers, raw_bam_files, input_ref_file, samtools, piccard_tools, GATK, para_jobs, script_name, userID ):
	"""! @brief do all preparation in one file per chromosome on cluster """
	
	IDs_to_check = []
	files_to_delete = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, header in enumerate( headers ):
		
		ID = "PP" + batch_ID + '_' + str( idx ).zfill(3)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		cmd = ''.join( [ 	"python3 " + script_name,
									" --header " + header,
									' --original_bam_file ' + input_bam_file,
									" --bam_file ",
									raw_bam_files[ idx ],
									" --ref_file ",
									input_ref_file,
									" --samtools " + samtools,
									" --GATK " + GATK,
									" --piccard_tools " + piccard_tools
								] )
		
		with open( sh_file, "w" ) as out:
			out.write( "#!/bin/bash\n")
			out.write( "#SBATCH --job-name=" + str(ID)  + "\n" )
			out.write( "#SBATCH --cpus-per-task=1\n" )
			out.write( "#SBATCH --mem=30G\n" )
			out.write( "#SBATCH --ntasks=1\n" )
			out.write( "#SBATCH --output=" + out_file + "\n")
			out.write( "#SBATCH --error=" + err_file + "\n" )
			out.write( cmd + "\n" )
		
		result=subprocess.run("chmod +x " + sh_file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		# Check if the command completed successfully
		if result.returncode!=0:
			sys.stdout.write("An error occurred while executing the command:" + result.stderr.decode())
			sys.stdout.flush()

		result=subprocess.run("sbatch " + sh_file, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		sys.stdout.write("JOB SUBMITTED (" + ID + ")\n")
		sys.stdout.flush()

		# Check if the command completed successfully
		if result.returncode!=0:
			sys.stdout.write("An error occurred while executing the command:"+result.stderr.decode())
			sys.stdout.flush()

		time.sleep( 1 )
		
		waiting_status = True
		while waiting_status:
			squeue = subprocess.run([ "squeue -u " + userID + ' --format="%30j"' ], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
			content = squeue.stdout
			qstat_IDs = re.findall("PP_" + batch_ID + "_\d{3}", content)
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		squeue = subprocess.run([ "squeue -u " + userID + ' --format="%30j"' ], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		content = squeue.stdout
		qstat_IDs = re.findall("PP_" + batch_ID + "_\d{3}", content)
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )


def run_haplotype_caller_on_cluster( prefix, bam_files, input_ref_file, GATK, para_jobs, userID ):
	"""! @brief run variant detection itself """
	
	IDs_to_check = []
	files_to_delete = []
	resulting_vcf_files = []
	
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, bam_file in enumerate( bam_files ):
		output_file = bam_file + "_raw_variants.vcf"
		resulting_vcf_files.append( output_file )
		
		ID = "HC" + batch_ID + '_' + str( idx ).zfill(3)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		files_to_delete.append( sh_file )
		out_file = prefix + ID + '.out'
		files_to_delete.append( out_file )
		err_file = prefix + ID + '.report'
		
		cmd = "java -Xmx8g -XX:ParallelGCThreads=4 -jar " + GATK + " HaplotypeCaller --reference " + input_ref_file + " --input " + bam_file + " --output " + output_file

		with open( sh_file, "w" ) as out:
			out.write( "#!/bin/bash\n")
			out.write( "#SBATCH --job-name=" + str(ID)  + "\n" )
			out.write( "#SBATCH --cpus-per-task=1\n" )
			out.write( "#SBATCH --mem=30G\n" )
			out.write( "#SBATCH --ntasks=1\n" )
			out.write( "#SBATCH --output=" + out_file + "\n")
			out.write( "#SBATCH --error=" + err_file + "\n" )
			out.write( "".join( cmd ) + "\n" )
		
		result=subprocess.run("chmod +x " + sh_file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		# Check if the command completed successfully
		if result.returncode!=0:
			sys.stdout.write("An error occurred while executing the command:" + result.stderr.decode())
			sys.stdout.flush()

		result=subprocess.run("sbatch " + sh_file, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		sys.stdout.write("JOB SUBMITTED (" + ID + ")\n")
		sys.stdout.flush()

		# Check if the command completed successfully
		if result.returncode!=0:
			sys.stdout.write("An error occurred while executing the command:"+result.stderr.decode())
			sys.stdout.flush()

		time.sleep( 1 )
		
		waiting_status = True
		while waiting_status:
			squeue = subprocess.run([ "squeue -u " + userID + ' --format="%30j"' ], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
			content = squeue.stdout
			qstat_IDs = re.findall("HC_" + batch_ID + "_\d{3}", content)
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		squeue = subprocess.run([ "squeue -u " + userID + ' --format="%30j"' ], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		content = squeue.stdout
		qstat_IDs = re.findall("HC_" + batch_ID + "_\d{3}", content)
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
		
	return resulting_vcf_files


def combine_vcf_files( prefix, vcf_files ):
	"""! @brief combines VCF files of all chromosomes for final filtering """
	
	combined_raw_vcf_file = prefix + "combined_raw_vcf_file.vcf"
	
	# --- collect header lines of all VCF files --- #
	fileformat_line = ""
	filter_lines = []
	format_lines = []
	GATK_command_lines = []
	info_lines = []
	contig_lines = []
	header = ""
	
	with open( vcf_files[0], "r" ) as f:
		line = f.readline()
		while line:
			if '##fileformat=' in line:
				fileformat_line = line
			#ALT
			elif '##FILTER=' in line:
				filter_lines.append( line )
			elif '##FORMAT=' in line:
				format_lines.append( line )
			elif '##GATKCommandLine.HaplotypeCaller=' in line:
				GATK_command_lines.append( line )
			##GVCFBlock
			elif '##INFO=' in line:
				info_lines.append( line )
			elif '##contig=' in line:
				contig_lines.append( line )
			##source
			elif '#CHROM' in line:
				header = '\t'.join( line.split('\t')[:-1] ) + '\tDATA_SET\n'
			
			if line[0] != '#':
				break
			line = f.readline()
	
	for vcf_file in vcf_files[1:]:
		with open( vcf_file, "r" ) as f:
			line = f.readline()
			while line:
				if '##GATKCommandLine.HaplotypeCaller=' in line:
					GATK_command_lines.append( line )
				
				if line[0] != '#':
					break
				line = f.readline()
	
	# --- write header lines into output VCF file --- #
	with open( combined_raw_vcf_file, "w" ) as out:
		out.write( fileformat_line )
		out.write( ''.join( filter_lines ) )
		out.write( ''.join( format_lines ) )
		out.write( ''.join( GATK_command_lines ) )
		out.write( ''.join( info_lines ) )
		out.write( ''.join( contig_lines ) )
		out.write( header )
		
		# --- write variant content into combined VCF file --- #
		number_of_raw_variants = 0
		for vcf_file in vcf_files:
			with open( vcf_file, "r" ) as f:
				line = f.readline()
				while line:
					if line[0] != '#':
						out.write( line )
						number_of_raw_variants += 1
					line = f.readline()
	sys.stdout.write( "number of raw variants: " + str( number_of_raw_variants ) + "\n" )
	return combined_raw_vcf_file


def load_vcf_data( vcf_file ):
	"""! @brief load all data from VCF file """
	
	data = []
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				data.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'line': line } )
			line = f.readline()
	return data


def construct_sorted_output_file( data, outputfile ):
	"""! @brief write all data in sorted way (chr, pos) into output file """
	
	sorted_data = sorted( data, key=itemgetter( 'chr', 'pos', 'line' ) )
	with open( outputfile, "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDATA_SET\n" )
		for each in sorted_data:
			out.write( each['line'] )


def separate_snps_and_indels( prefix, combined_raw_vcf_file, input_ref_file, GATK ):
	"""! @brief separation of SNPs and InDels and separate filtering of both data sets """
	
	snps_only_file = prefix + "raw_snps_only_file.vcf"
	indels_only_file = prefix + "raw_indels_only_file.vcf"
	
	clean_snps_file = prefix + "clean_snps_file.vcf"
	clean_indels_file = prefix + "clean_indels_file.vcf"
	
	clean_all_variants_file = prefix + "clean_all_variants_file.vcf"
	
	
	# --- extraction of SNPs --- #
	cmd = ''.join( [ 	"java -Xmx32g -XX:ParallelGCThreads=1 -jar ",
						GATK,
						" SelectVariants -R ",
						input_ref_file,
						" -V ",
						combined_raw_vcf_file,
						" --select-type-to-include SNP --output ",
						snps_only_file
					] )
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()
	
	# --- filtering of SNPs --- #
	cmd = ''.join( [ 	"java -Xmx32g -XX:ParallelGCThreads=1 -jar ",
						GATK,
						" VariantFiltration -R ",
						input_ref_file,
						" -V ",
						snps_only_file,
						' --filter-expression "QD < 2.0" --filter-name "QD_filter"',
						' --filter-expression "FS > 60.0" --filter-name "FS_filter"',
						' --filter-expression "MQ < 40.0" --filter-name "MQ_filter"',
						" --output ",
						clean_snps_file
					] )
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()
	
	# --- extraction of InDels --- #
	cmd = ''.join( [ 	"java -Xmx32g -XX:ParallelGCThreads=1 -jar ",
						GATK,
						" SelectVariants -R ",
						input_ref_file,
						" -V ",
						combined_raw_vcf_file,
						" --select-type-to-include INDEL --output ",
						indels_only_file
					])
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()
	
	# --- filtering of InDels --- #
	cmd = ''.join( [ 	"java -Xmx32g -XX:ParallelGCThreads=1 -jar ",
						GATK,
						" VariantFiltration -R ",
						input_ref_file,
						" -V ",
						indels_only_file,
						' --filter-expression "QD < 2.0" --filter-name "QD_filter"',
						' --filter-expression "FS > 200.0" --filter-name "FS_filter"',
						#' --filterExpression "DP > 300" --filterName "high_DP_filter"',
						#' --filterExpression "DP < 30" --filterName "low_DP_filter"',
						' --output ',
						clean_indels_file
					] )
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()
	
	# --- combine all clean variants in one vcf file --- #
	snps = load_vcf_data( clean_snps_file )
	indels = load_vcf_data( clean_indels_file )
	construct_sorted_output_file( snps+indels, clean_all_variants_file )
	
	
	# --- print final report --- #
	sys.stdout.write( "\n\nFINAL REPORT\n\ncleaned SNPs are located here: " + clean_snps_file + "\n" )
	sys.stdout.write(  "\ncleaned InDels are located here: " + clean_indels_file + "\n" )
	sys.stdout.write(  "\nall cleaned variants are located here: " + clean_all_variants_file + "\n" )
	sys.stdout.flush()


def get_header_strings( header, length, block_size ):
	"""! @brief calculate header strings based on given chromosome name and lengths """
	
	strings = []
	start = 1
	end = 1 + block_size
	while end < length:
		strings.append( header + ":" + str( int( start ) ) + "-" + str( end ) )
		start += ( block_size * 0.99 )
		end += block_size
	strings.append( header + ":" + str( int( start ) ) + "-" + str( length ) )
	return strings



def load_fasta_headers( input_fasta ):
	"""! @brief load fasta headers as they occur within file """
	
	headers = []
	
	with open( input_fasta, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == '>':
				headers.append( line.strip()[1:].split(' ')[0] )
			line = f.readline()
	return headers


def load_vcf_file( input_vcf ):
	"""! @brief load VCF file into list of dictionaries """
	
	comment_lines = []
	variants = []
	
	with open( input_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == '#':
				comment_lines.append( line )
			else:
				parts = line.strip().split('\t')
				try:
					variants.append( { 'chr': parts[0], 'pos': int( parts[1] ), 'line': line } )
				except:
					pass
			line = f.readline()
	return sorted( variants, key=itemgetter( 'chr', 'pos' ) ), comment_lines


def construct_sorted_vcf_file( output_vcf, variants,  headers, comment_lines ):
	"""! @brief construct sorted VCF file """
	
	with open( output_vcf, "w" ) as out:
		for line in comment_lines:
			out.write( line )
		for header in headers:
			for variant in variants:
				if header == variant['chr']:
					out.write( variant['line'] )


def sort_vcf_by_fasta( input_fasta, input_vcf,  output_vcf ):
	"""! @brief calls all functions for sorting VCF files """
	
	variants, comment_lines = load_vcf_file( input_vcf )
	headers = load_fasta_headers( input_fasta )
	construct_sorted_vcf_file( output_vcf, variants,  headers, comment_lines )


def main( arguments):
	"""! @brief calls all functions involved in the GATK variant detection """
	
	input_bam_file = arguments[ arguments.index( '--in' )+1 ]
	input_ref_file = arguments[ arguments.index( '--ref' )+1 ]
	userID = arguments[ arguments.index( '--user' )+1 ]
	
	prefix = arguments[ arguments.index( '--out' )+1 ]
	if prefix[-1] != '/':
		prefix += "/"
	
	if '--para_jobs' in arguments:
		para_jobs = int( arguments[ arguments.index( '--para_jobs' )+1 ] )
	else:
		para_jobs=300
	
	if '--block' in arguments:
		block_size = int( arguments[ arguments.index( '--block' )+1 ] )
		if block_size < 1000000:
			block_size = 1000000
	else:
		block_size = 5000000	#5Mbp default block size; min block size 1Mbp
	
	
	if '--script' in arguments:
		script_name = arguments[ arguments.index( '--script' )+1 ]
	else:
		script_name = "variant_call_preparation.py"
	
	if '--samtools' in parameters:
		samtools = parameters[ parameters.index( '--samtools' )+1 ]
	else:
		samtools = "samtools"	#v2.5.0
	
	if '--picard' in parameters:
		piccard_tools = parameters[ parameters.index( '--picard' )+1 ]
	else:
		piccard_tools = "picard.jar"	#v4.1.9.0
	
	if '--gatk' in parameters:
		GATK = parameters[ parameters.index( '--gatk' )+1 ]
	else:
		GATK = "gatk-package-4.1.9.0-local.jar"
	
	
	# --- construct directory for all files --- #
	if os.path.exists( prefix ):
		sys.stdout.write( "directory already present\n" )
		sys.stdout.flush()
	else:
		os.makedirs( prefix )
	
	# -- transfering and indexing reference file --- #
	cmd = "cp " + input_ref_file + " " + prefix
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()
	input_ref_file = prefix + input_ref_file.split('/')[-1]
	if sum( [ os.path.isfile( input_ref_file + ".fai" ), os.path.isfile( input_ref_file.replace( ".fasta", ".dict" ) ) ] ) < 2:
		index_ref_file( input_ref_file, piccard_tools, samtools )
	
	
	# --- sorting and indexing BAM file prior to splitting --- #
	if not '--bam_is_sorted' in arguments:
		cmd = "cp " + input_bam_file + " " + prefix
		os.popen( cmd )
		sorted_bam_file = prefix + input_bam_file.split('/')[-1] + "_sorted.bam"
		sort_bam_file( prefix, input_bam_file, samtools, piccard_tools )
	else:
		sorted_bam_file = input_bam_file
	

	# --- delete all existing bai files to avoid outdated index files --- #
	bai_files = glob.glob( prefix + '*.bai' )
	for each in bai_files:
		os.popen( "rm " + each )
	sys.stdout.write( "number of deleted .bai files: " + str( len( bai_files ) ) + "\n" )
	sys.stdout.flush()
	time.sleep( 10 )
	
	# --- index input BAM file --- #
	if not os.path.isfile( sorted_bam_file + ".bai" ):
		index_bam_file( sorted_bam_file, piccard_tools )
	
	# --- preparation of splitting and processing separate BAM files on cluster --- #
	tmp_headers = load_sequences( input_ref_file ).keys()
	raw_bam_files = []
	headers = []
	
	vcf_files = glob.glob( prefix + "*.vcf" )
	black_list = {}
	for vcf in vcf_files:
		black_list.update( { vcf.split('/')[-1].split('.bam')[0]: None } )
	
	seq_lengths = get_seq_lengths( input_ref_file )	#get lengths of all sequences in FASTA file
	
	for header in tmp_headers:	#everything is done for a single chromosome/block
		try:
			black_list[ header ]
		except KeyError:
			header_strings = get_header_strings( header, seq_lengths[ header ], block_size )
			for each in header_strings:
				try:
					black_list[ each.replace(':', 'x_x_x') ]
				except KeyError:
					headers.append( each )
					bam_file = prefix + each.replace(':', 'x_x_x') + ".bam"
					raw_bam_files.append( bam_file )
	
	sys.stdout.write( "number of preparation jobs to run: " + str( len( headers ) ) + "\n" )
	sys.stdout.flush()
	
	# --- prepare files --- #
	if len( headers ) > 0:
		prepare_files_for_variant_calling_on_cluster( prefix, sorted_bam_file, headers, raw_bam_files, input_ref_file, samtools, piccard_tools, GATK, para_jobs, script_name, userID )
	
	
	# --- identify varaints via HaplotypeCaller --- #
	resulting_bam_files = glob.glob( prefix + "*_RG_added.bam" )
	bams_to_do = []
	for bam in resulting_bam_files:
		ID = bam.split('/')[-1].split('.bam')[0]
		try:
			black_list[ ID ]
		except KeyError:
			bams_to_do.append( bam )
	
	x =  run_haplotype_caller_on_cluster( prefix, bams_to_do, input_ref_file, GATK, para_jobs, userID )
	
	
	# --- processing of VCFs --- #
	combined_raw_vcf_file = prefix + "combined_raw_vcf_file.vcf"
	if not os.path.exists( combined_raw_vcf_file ):
		tmp_raw_vcf_files = sorted( glob.glob( prefix + "*.vcf" ) )
		raw_vcf_files = []
		for each in tmp_raw_vcf_files:
			if each.split('/')[-1].split('.')[0] not in [ "clean_all_variants_file", "clean_indels_file", "clean_snps_file", "combined_raw_vcf_file", "raw_indels_only_file", "raw_snps_only_file" ]:
				raw_vcf_files.append( each )
		sys.stdout.write( "number of generated VCFs in this round: " + str( len( x ) ) + "\n" )
		sys.stdout.write( "total number of present VCFs: " + str( len( raw_vcf_files ) ) + "\n" )
		sys.stdout.flush()
		combined_raw_vcf_file = combine_vcf_files( prefix, raw_vcf_files )
	
	
	# --- sort variants in combined vcf file --- #
	sorted_combined_raw_vcf_file = combined_raw_vcf_file + "_sorted.vcf"
	if not os.path.isfile( sorted_combined_raw_vcf_file ):
		sort_vcf_by_fasta( input_ref_file, combined_raw_vcf_file,  sorted_combined_raw_vcf_file )
	
	
	# --- separate and filter variants --- #
	separate_snps_and_indels( prefix, sorted_combined_raw_vcf_file, input_ref_file, GATK )
	
	
	# --- analysing number of completed jobs --- #
	sys.stdout.write( "Process completed:\n" )
	tmp_raw_vcf_files = sorted( glob.glob( prefix + "*.vcf" ) )
	raw_vcf_files = []
	for each in tmp_raw_vcf_files:
		if each.split('/')[-1].split('.')[0] not in [ "clean_all_variants_file", "clean_indels_file", "clean_snps_file", "combined_raw_vcf_file", "raw_indels_only_file", "raw_snps_only_file" ]:
			raw_vcf_files.append( each )
	sys.stdout.write( "number of generated VCFs in this round: " + str( len( x ) ) + "\n" )
	sys.stdout.write( "total number of present VCFs: " + str( len( raw_vcf_files ) ) + "\n" )
	sys.stdout.flush()


if '--in' in sys.argv and '--ref' in sys.argv and '--out' in sys.argv and '--user' in sys.argv:
	main( sys.argv)
else:
	sys.exit( __usage__ )
