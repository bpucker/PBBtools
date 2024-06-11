### Boas Pucker ###
### b.pucker@tu-bs.de ###
### v0.3 ###

__usage__ = """
				PE mode:
				python3 rename_reads_for_trinity.py
				--fwin <FW_READ_INPUT_FILENAME>
				--rvin <RV_READ_INPUT_FILENAME>
				--fwout <FW_READ_OUTPUT_FILENAME>
				--rvout <RV_READ_OUTPUT_FILENAME>
				
				SE mode:
				python3 rename_reads_for_trinity.py
				--fwin <FW_READ_INPUT_FILENAME>
				--fwout <FW_READ_OUTPUT_FILENAME>
				"""

import sys, gzip

# --- end of imports --- #

def main( arguments ):
	"""! @brief generate read names for fw and rv read based on given read number in file """
	
	input_fw_file = arguments[ arguments.index( '--fwin' )+1 ]
	if '--rvin' in arguments:
		input_rv_file = arguments[ arguments.index( '--rvin' )+1 ]
	else:
		input_rv_file = ""		

	output_fw_file = arguments[ arguments.index( '--fwout' )+1 ]
	if '--rvout' in arguments:
		output_rv_file = arguments[ arguments.index( '--rvout' )+1 ]
	else:
		output_rv_file = ""
	
	
	#example: @NS500530:22:HVLCMBGX5:1:11101:10266:1041 1:N:0:CTTGTA
	#x: 1-99999
	#y: 1-99999
	#tile numer:1000-9999
	#lane: 1-8
	#flowcell ID: diverse string => get list of some
	#run number
	
	x = 1
	y = 1
	tile = 1000
	lane = 1
	flowcell = "HVLCMBGX5"
	run = 1
	
	#PE mode
	if len( input_rv_file ) > 0:
		with gzip.open( output_fw_file, "wb" ) as out_fw:
			with gzip.open( output_rv_file, "wb" ) as out_rv:
				with gzip.open( input_fw_file, "rb" ) as f1:
					with gzip.open( input_rv_file, "rb" ) as f2:
						line = f1.readline()
						while line:
							# --- generate read ID and update all variables --- #
							ID = ":".join( list( map( str, [ "@NS500530", run, flowcell, lane, tile, x, y ] ) ) )
							x += 1
							if x == 99999:
								x = 1
								y += 1
								if y == 99999:
									y = 1
									tile += 1
									if tile == 9999:
										tile = 1000
										lane += 1
										if lane == 8:
											lane = 1
											run += 1
											
							
							# --- read block from read file1 --- #
							header1 =  ID + ' 1:N:0:ACACAG\n'
							seq1 = f1.readline()
							f1.readline()
							qual1 = f1.readline()
							
							# --- read block from read file2 --- #
							f2.readline()
							header2 =  ID + ' 2:N:0:ACACAG\n'
							seq2 = f2.readline()
							f2.readline()
							qual2 = f2.readline()
							
							if len( seq1 ) == len( qual1 ) and len( seq1 ) > 0:
								if len( seq2 ) == len( qual2 ) and len( seq2 ) > 0:
									out_fw.write( header1 + seq1 + "+\n" + qual1  )
									out_rv.write( header2 + seq2 + "+\n" + qual2  )
							
							line = f1.readline()
	#SE mode
	else:
		with gzip.open( output_fw_file, "wb" ) as out_fw:
			with gzip.open( input_fw_file, "rb" ) as f1:
				line = f1.readline()
				while line:
					# --- generate read ID and update all variables --- #
					ID = ":".join( list( map( str, [ "@NS500530", run, flowcell, lane, tile, x, y ] ) ) )
					x += 1
					if x == 99999:
						x = 1
						y += 1
						if y == 99999:
							y = 1
							tile += 1
							if tile == 9999:
								tile = 1000
								lane += 1
								if lane == 8:
									lane = 1
									run += 1
					# --- write read record into output file --- #
					header1 =  ID + ' 1:N:0:ACACAG\n'
					seq1 = f1.readline()
					f1.readline()
					qual1 = f1.readline()
					
					if len( seq1 ) == len( qual1 ) and len( seq1 ) > 0:
						out_fw.write( header1 + seq1 + "+\n" + qual1  )
					
					line = f1.readline()

if '--fwin' in sys.argv and '--fwout' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

