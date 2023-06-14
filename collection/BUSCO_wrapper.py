### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python BUSCO_wrapper.py
					--out <FULL_PATH_TO_OUTPUT_DIR>
					--in <INPUT_ASSEMBLY_FILE>
					optional:
					--mode < prot | tran | geno >[geno]
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os, re, glob
from operator import itemgetter

# --- end of imports --- #

def run_BUSCO( input_file, prefix, busco_path, augustus_path, augustus_config_path, mode ):
	"""! @brief run BUSCO in genome mode on given assembly file """
	
	os.chdir( prefix )
	os.environ["PATH"] = augustus_path + ":" + os.environ["PATH"]
	print os.environ["PATH"]
	os.environ["AUGUSTUS_CONFIG_PATH"] = augustus_config_path
	print os.environ["AUGUSTUS_CONFIG_PATH"]
	tmp_path = prefix + "tmp"
	if not os.path.exists( tmp_path ):
		os.makedirs( tmp_path )
	cmd = "python " + busco_path + " --in " + input_file + " --tmp_path " + tmp_path + " --cpu 10 --mode " + mode + " --long --force --out busco_run > " + prefix +"log.txt"
	os.popen( cmd )


def main( parameters ):
	"""! @brief run all genome analysis methods """
	
	working_dir = parameters[ parameters.index('--out')+1 ]
	input_file = parameters[ parameters.index('--in')+1 ]
	if '--mode' in parameters:
		mode = parameters[ parameters.index('--mode')+1 ]
	else:
		mode = "genome"	#prot | tran | geno
		
	augustus = "/vol/biotools/bin/augustus-3.2"
	augustus_seqs_ex_script = "/grp/pbb/scripts/getAnnoFasta.pl"
	
	busco = "/grp/pbb/tools/BUSCO/scripts/run_BUSCO.py"
	augustus_path = "/vol/biotools/lib/augustus-3.2.2/bin/"
	augustus_config_path = "/grp/pbb/tools/augustus-3.2.2/config/"
	
	
	if working_dir[-1] != '/':
		working_dir += "/"
	if not os.path.exists( working_dir ):
		os.makedirs( working_dir )
	
	# --- transfer file to cluster --- #
	assembly_file = working_dir + input_file.split('/')[-1]
	os.popen( "cp " + input_file + " " + assembly_file )
	
	
	# --- run BUSCO v3 --- #
	prefix = working_dir + "busco/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	run_BUSCO( input_file, prefix, busco, augustus_path, augustus_config_path, mode )


if __name__ == '__main__':
	
	if '--out' in sys.argv and '--in' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
