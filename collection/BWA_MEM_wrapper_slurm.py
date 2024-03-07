### Boas Pucker ###
### v0.9 ###
### b.pucker@tu-bs.de ###

__usage__ = """ python3 BWA_MEM_wrapper_slurm.py
				MANDATORY:
				--output_dir <FULL_PATH_TO_DIRECTORY>
				--reference <FULL_PATH_TO_FILE>
				--read1 <FASTQ_FILES_OF_FW_READS_SEPARATED_BY_COMMA>
				--user <USER_ID>
				
				OPTIONAL:
				--read2 <FASTQ_FILES_OF_RV_READS_SEPARATED_BY_COMMA>
				--reference_is_indexed [PREVENTS_CONSTRUCTION_OF_NEW_REFERENCE_INDEX]
				--bwa <BWA_MEM_PATH>[bwa]
				--samtools <SAMTOOLS_PATH>[samtools]
				--picard <PICARD_PATH>[picard.jar]
				--jobs <NUMBER_OF_PARALLEL_JOBS>[20]
				
				Paired-end mapping is performed if two files are provided. Otherwise single end mapping is performed with data provided via --read1.
				
				feature requests and bug reports: b.pucker@tu-bs.de
			"""


import os, re, sys, glob, time, datetime, subprocess 

# --- end of imports --- #


def submit_BWA_MEM( BWA_path, reference, output_directory, read_file_pairs, para_jobs, userID ):
	"""! @brief submit BWA MEM mapping jobs to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-4:]
	for idx, read_file_pair in enumerate( sorted( read_file_pairs ) ):	#iterate over all provided data files
		time.sleep( 1 )
		
		ID = "W_" + batch_ID + '_' + str( idx ).zfill(3)
		IDs_to_check.append( ID )
		
		# --- setting file names --- #
		err_file = output_directory + ID + ".err"
		out_file = output_directory + ID + ".out"
		sh_file = output_directory + ID + ".sh"
		
		output_file = output_directory + ID + ".sam.gz"
		
		if len( read_file_pair ) > 1:
			BWA_cmd = [ 	BWA_path,
							"mem -M -t 11 -c 1000",
							reference,
							read_file_pair[0],
							read_file_pair[1],
							"| gzip -3 >",
							output_file
						]
		else:
			BWA_cmd = [ 	BWA_path,
							"mem -M -t 11 -c 1000",
							reference,
							read_file_pair[0],
							"| gzip -3 >",
							output_file
						]
		BWA_cmd = " ".join( BWA_cmd )
		
		with open( sh_file, "w" ) as out:
			out.write( "#!/bin/bash\n")
			out.write( "#SBATCH --job-name=" + str(ID)  + "\n" )
			out.write( "#SBATCH --cpus-per-task=1\n" )
			out.write( "#SBATCH --mem=30G\n" )
			out.write( "#SBATCH --ntasks=1\n" )
			out.write( "#SBATCH --output=" + out_file + "\n")
			out.write( "#SBATCH --error=" + err_file + "\n" )
			out.write( BWA_cmd + "\n" )
		
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
			qstat_IDs = re.findall("W_" + batch_ID + "_\d{3}", content)
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
		qstat_IDs = re.findall("W_" + batch_ID + "_\d{3}", content)
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
	
	sys.stdout.write( "BWA MEM read mapping done...\n" )
	sys.stdout.flush()


def convert_sam_to_bam( bwa_result_files, samtools_path, output_directory, para_jobs, userID ):
	"""! @brief submits jobs to cluster to convert sam to bam """
	
	#print "submitting BWA MEM jobs to cluster ... "
	time.sleep( 10 )
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-4:]
	for idx, filename in enumerate( sorted( bwa_result_files ) ):	#iterate over all provided data files
		time.sleep( 1 )
		
		ID = "C_" + batch_ID + '_' + str( idx ).zfill(3)
		IDs_to_check.append( ID )
		
		# --- setting file names --- #
		err_file = output_directory + ID + ".err"
		out_file = output_directory + ID + ".out"
		
		sh_file = output_directory + ID + ".sh"
		
		output_file = output_directory + ID + ".sam.gz"
		
		samtools_cmd = [ 	samtools_path,
							" view -Sb -@ 8 -m 20G ",
							filename,
							" > ",
							filename.split('.sam.gz')[0] + ".bam.gz"
						]
		samtools_cmd = "".join( samtools_cmd )
		
		with open( sh_file, "w" ) as out:
			out.write( "#!/bin/bash\n")
			out.write( "#SBATCH --job-name=" + str(ID)  + "\n" )
			out.write( "#SBATCH --cpus-per-task=1\n" )
			out.write( "#SBATCH --mem=30G\n" )
			out.write( "#SBATCH --ntasks=1\n" )
			out.write( "#SBATCH --output=" + out_file + "\n")
			out.write( "#SBATCH --error=" + err_file + "\n" )
			out.write( samtools_cmd + "\n" )
		
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
			qstat_IDs = re.findall("C_" + batch_ID + "_\d{3}", content)
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
		qstat_IDs = re.findall("C_" + batch_ID + "_\d{3}", content)
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 10 )
	
	sys.stdout.write( "conversion from SAM to BAM done...\n" )
	sys.stdout.flush()


def main( parameters ):
	"""! @brief run everything """
	
	prefix = parameters[ parameters.index( '--out' )+1 ]
	if prefix[-1] != '/':
		prefix += '/'
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	input_reference_file = parameters[ parameters.index( '--ref' )+1 ]
	userID = parameters[ parameters.index( '--user' )+1 ]
	
	if '--jobs' in parameters:
		try:
			para_jobs = int( parameters[ parameters.index( '--jobs' )+1 ] )
		except:
			para_jobs = 20
	else:
		para_jobs = 20
	
	if '--bwa' in parameters:
		BWA_path = parameters[ parameters.index( '--bwa' )+1 ]
	else:
		BWA_path = "bwa"
	
	if '--samtools' in parameters:
		samtools_path = parameters[ parameters.index( '--samtools' )+1 ]
	else:
		samtools_path = "samtools"
	
	if '--picard' in parameters:
		picard_tools_path = parameters[ parameters.index( '--picard' )+1 ]
	else:
		picard_tools_path = "picard.jar"
	
	read1_file = parameters[ parameters.index( '--read1' )+1 ]
	if "," in read1_file:
		read_files1 = read1_file.split( ',' )
	else:
		read_files1 = [ read1_file ]
	try:
		read2_file = parameters[ parameters.index( '--read2' )+1 ]
		if "," in read2_file:
			read_files2 = read2_file.split( ',' )
		else:
			read_files2 = [ read2_file ]
	except:
		read2_file = False
	
	
	# --- producing output directory --- #
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	
	# -- continue with original input file(s) --- #
	tmp_dir = prefix + str( datetime.datetime.now() )[-4:] + "_tmp_dir/"
	splitted_read_pair_files = []
	RG_IDs = []
	
	for i in range( len( read_files1 ) ):
		if read2_file:
			if os.path.isfile( read_files1[i] ) and os.path.isfile( read_files2[i] ):
				splitted_read_pair_files += [ ( read_files1[i], read_files2[i] ) ]
				RG_IDs.append( read_files1[i].split('/')[-1].split('.')[0] )
		else:
			if os.path.isfile( read_files1[i] ):
				splitted_read_pair_files += [ [ read_files1[i] ] ]
				RG_IDs.append( read_files1[i].split('/')[-1].split('.')[0] )

	# --- construct reference if necessary --- #
	if not '--reference_is_indexed' in parameters:
		sys.stdout.write( "constructing reference ... \n" )
		sys.stdout.flush()
		reference_file = prefix + input_reference_file.split('/')[-1]
		cmd = "cp " + input_reference_file + " " + reference_file
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
		
		cmd = BWA_path + " index " + reference_file
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	else:
		reference_file = input_reference_file
	
	
	# --- submit mapping jobs to cluster --- #
	sys.stdout.write( "running BWA ... \n" )
	sys.stdout.flush()
	bwa_mapping_result_path = prefix + "bwa_mapping_result_path/"
	if not os.path.exists( bwa_mapping_result_path ):
		os.makedirs( bwa_mapping_result_path )
	bwa_result_files = glob.glob( bwa_mapping_result_path + "*.sam.gz" )
	if len( bwa_result_files ) == 0:
		submit_BWA_MEM( BWA_path, reference_file, bwa_mapping_result_path, splitted_read_pair_files, para_jobs, userID )
	
	# --- #convert to bam --- #
	bwa_result_files = glob.glob( bwa_mapping_result_path + "*.sam.gz" )
	sys.stdout.write( "number of identified BWA result files: " + str( len( bwa_result_files ) ) + "\n" )
	sys.stdout.flush()
	bam_files = glob.glob( bwa_mapping_result_path + "*.bam.gz" )
	if len( bam_files ) == 0:
		convert_sam_to_bam( bwa_result_files, samtools_path, bwa_mapping_result_path, para_jobs, userID )
	else:
		if os.path.getsize( bam_files[0] ) < 1000:
			convert_sam_to_bam( bwa_result_files, samtools_path, bwa_mapping_result_path, para_jobs, userID )
	
	
	# --- merge results --- #
	sys.stdout.write( "merging BAM files ... \n" )
	sys.stdout.flush()
	final_bam_file = prefix + "final_bam_file.bam.gz"
	cmd = samtools_path + " merge " + final_bam_file + " " + bwa_mapping_result_path + "*.bam.gz"
	if not os.path.isfile( final_bam_file ):
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	elif os.path.getsize( final_bam_file ) < 1000:
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	
	# --- sort by position --- #
	sys.stdout.write( "sorting BAM file ... \n" )
	sys.stdout.flush()
	sorted_final_bam_file = final_bam_file.replace( ".bam.gz", "_sorted.bam.gz" )
	cmd = samtools_path + " sort -@ 4 -m 3G -o " +  sorted_final_bam_file + " " + final_bam_file
	if not os.path.isfile( sorted_final_bam_file ):
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	elif os.path.getsize( sorted_final_bam_file ) < 1000:
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	
	# --- remove duplicates --- #
	sys.stdout.write( "removing duplicate reads ... \n" )
	sys.stdout.flush()
	removed_duplicates_file = sorted_final_bam_file.replace( ".bam.gz", "_duplicates_removed.bam.gz" )
	metrics_file = prefix + "remove_duplicates_metrics.txt"
	cmd = [ "java -Xmx11g -jar ",
			picard_tools_path,
			" MarkDuplicates I=",
			sorted_final_bam_file,
			" O=" + removed_duplicates_file,
			" METRICS_FILE=" + metrics_file,
			" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=5000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000"
		]
	if not os.path.isfile( removed_duplicates_file ):
		p = subprocess.Popen( args="".join( cmd ), shell=True )
		p.communicate()
	elif os.path.getsize( removed_duplicates_file ) < 1000:
		p = subprocess.Popen( args="".join( cmd ), shell=True )
		p.communicate()

	# --- sort by position --- #
	sys.stdout.write( "sorting BAM file ... \n" )
	sys.stdout.flush()
	final_removed_duplicates_file = removed_duplicates_file.replace( ".bam.gz", "_sorted.bam.gz" )
	cmd = samtools_path + " sort -@ 4 -m 3G -o " +  final_removed_duplicates_file + " " + removed_duplicates_file
	if not os.path.isfile( final_removed_duplicates_file ):
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	elif os.path.getsize( final_removed_duplicates_file ) < 1000:
		p = subprocess.Popen( args=cmd, shell=True )
		p.communicate()
	
	
	# --- adding read group IDs --- #
	sys.stdout.write( "adding read group ... \n" )
	sys.stdout.flush()
	rg_final_file = final_removed_duplicates_file.replace( ".bam.gz", "_rg.bam.gz" )
	cmd = [ "java -Xmx8g -jar ",
			picard_tools_path,
			" AddOrReplaceReadGroups I=",
			final_removed_duplicates_file,
			" O=",
			rg_final_file,
			" RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"
		]
	if not os.path.isfile( rg_final_file ):
		p = subprocess.Popen( args="".join( cmd ), shell=True )
		p.communicate()
	elif os.path.getsize( rg_final_file ) < 1000:
		p = subprocess.Popen( args="".join( cmd ), shell=True )
		p.communicate()


if '--out' in sys.argv and '--ref' in sys.argv and '--read1' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )	
