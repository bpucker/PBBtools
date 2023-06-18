### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.2 ###

### Pucker & Brockington, 2018. doi: https://doi.org/10.1186/s12864-018-5360-z ###

__usage__ = """
					python run_BUSCO_across_species.py
					--in <PEPTIDE_FILE_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					
					optional:
					--augustus <AUGUSTUS_PATH>
					--augustusseqex <AUGUSTUS_SEQ_EXTRACTION_SEQUENCES>
					--busco <BUSCO>
					--augustusbinpath <AUGUSTUS_BIN_PATH>
					--augustusconfigpath <AUGUSTUS_CONFIG_PATH>
					"""

import sys, os, re, glob
from operator import itemgetter

# --- end of imports --- #

def run_BUSCO( input_file, prefix, busco_path, augustus_path, augustus_config_path ):
	"""! @brief run BUSCO in genome mode on given assembly file """
	
	os.chdir( prefix )
	os.environ["PATH"] = augustus_path + ":" + os.environ["PATH"]
	print os.environ["PATH"]
	os.environ["AUGUSTUS_CONFIG_PATH"] = augustus_config_path
	print os.environ["AUGUSTUS_CONFIG_PATH"]
	cmd = "python " + busco_path + " --in " + input_file + " --cpu 10 -m proteins --out busco_run > " + prefix +"log.txt"
	os.popen( cmd )


def main( arguments ):
	"""! @brief runts everything """
	
	input_dir = arguments[ arguments.index('--in')+1 ]
	output_dir = arguments[ arguments.index('--out')+1 ]
	
	if output_dir[-1] != "/":
		output_dir += "/"
	
	output_file = output_dir + "SUMMARY.txt"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	if '--augustus' in arguments:
		augustus = arguments[ arguments.index('--augustus')+1 ]
	else:
		augustus = "/vol/biotools/bin/augustus-3.2"
	
	if '--augustusseqex' in arguments:
		augustus_seqs_ex_script = arguments[ arguments.index('--augustusseqex')+1 ]
	else:
		augustus_seqs_ex_script = "/grp/pbb/tools/getAnnoFasta.pl"
	
	if '--busco' in arguments:
		busco = arguments[ arguments.index('--busco')+1 ]
	else:
		busco = "/grp/pbb/tools/BUSCO/scripts/run_BUSCO.py"
	
	if '--augustusbinpath' in arguments:
		augustus_path = arguments[ arguments.index('--augustusbinpath')+1 ]
	else:
		augustus_path = "/grp/pbb/tools/augustus-3.2.2/bin/"
	
	if '--augustusconfigpath' in arguments:
		augustus_config_path = arguments[ arguments.index('--augustusconfigpath')+1 ]
	else:
		augustus_config_path = "/grp/pbb/tools/augustus-3.2.2/config/"
	
	active = True
	
	filenames = glob.glob( input_dir + "*.pep.fasta" ) + glob.glob( input_dir + "*.pep.fa" )  + glob.glob( input_dir + "*.fa" )    + glob.glob( input_dir + "*.fasta" )  
	results = {}
	for filename in sorted( filenames ):
		ID = filename.split('/')[-1].split('.')[0]
		working_dir = output_dir + ID + '/'
		prefix = working_dir + "busco/"
		if not os.path.exists( working_dir ):
			os.makedirs( working_dir )
			os.makedirs( prefix )
			# --- run BUSCO v3 --- #
			if active:
				run_BUSCO( filename, prefix, busco, augustus_path, augustus_config_path )
		
		# ---- collect results --- #
		log_file = prefix + "log.txt"
		with open( log_file, "r" ) as f:
			content = f.read()
		try:
			result_string = re.findall( "C:\d+\.\d+%\[S:\d+\.\d+%,D:\d+\.\d+%\],F:\d+\.\d+%,M:[-]*\d+\.\d+%,n:\d+", content )[0]
			results.update( { ID: result_string } )
		except IndexError:
			print "ERROR: " + ID
	
	with open( output_file, "w" ) as out:
		for key in sorted( results.keys() ):
			out.write( key + '\t' + results[ key ] + '\n' )


if '--out' in sys.argv and '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
