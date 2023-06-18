#Derived from a script collection provided by Ya Yang

__usage__ = """
					python mask_monophyly_bp.py
					--tre <TREE_FILE>
					--aln <ALIGNMENT_FILE>
					
					optional:
					--seq <SEQ_FILE>
					--black <BLACK_LIST_FILE>
					"""

import sys, os, argparse, newick3, phylo3, re
from tree_utils import get_name,remove_kink
from utils import parse_fasta

# --- end of imports --- #

def get_names_to_exclude(ignoref):
    ignore = [l.strip() for l in open(ignoref,"r").readlines()]
    return ignore


def mask_monophyletic_tips(curroot,unamb_chrDICT,ignore=[]):
	going = True
	while going and curroot != None and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): # walk through nodes
			if not node.istip: continue # only look at tips
			name = get_name(node.label)
			if name in ignore: continue # do not mask the genomes
			for sister in node.get_sisters():
				if sister.istip and name==get_name(sister.label): # mask
					if unamb_chrDICT[node.label] > unamb_chrDICT[sister.label]:
						node = sister.prune()			
					else: node = node.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot


def mask_paraphyletic_tips(curroot,unamb_chrDICT,ignore=[]):
	going = True
	while going and curroot != None and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): #walk through nodes
			if not node.istip: continue #only look at tips
			name = get_name(node.label)
			if name in ignore: continue # do not mask the genomes
			parent = node.parent
			if node == curroot or parent == curroot or parent == None:
				continue #no paraphyletic tips for the root
			for para in parent.get_sisters():
				if para.istip and name==get_name(para.label): # mask
					if unamb_chrDICT[node.label] > unamb_chrDICT[para.label]:
						node = para.prune()
					else: node = node.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot


def mask(curroot, clnfile, para=True,ignore=[]):
	chrDICT = {} #key is seqid, value is number of unambiguous chrs
	for key, value in dict([x for x in parse_fasta(clnfile)]).items():
		for ch in ['-','X',"x","?","*"]:
			value = value.replace(ch,"") #ignore gaps, xs and Xs
		chrDICT[key] = len(value)
	curroot = mask_monophyletic_tips(curroot,chrDICT,ignore)
	if para:
		curroot = mask_paraphyletic_tips(curroot,chrDICT,ignore)
	return curroot


def mask_monophyly(tre,clnaln,para=True,ignore=[]):
    with open(tre,"r") as inf:
        intree = newick3.parse(inf.readline())
    curroot = mask(intree,clnaln,para,ignore)
    with open(tre+".mm","w") as outf:
        outf.write(newick3.tostring(curroot)+";\n")


def load_sequences( fasta_file, pattern ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				try:
					ID = re.findall( pattern, header )[0]
					sequences.update( { ID: "".join( seq ) } )
				except IndexError:
					print ( header )
				header = line.strip()[1:]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		try:
			ID = re.findall( pattern, header )[0]
			sequences.update( { ID: "".join( seq ) } )
		except IndexError:
			print ( header )
	return sequences


def main( arguments ):
	"""! @brief run everything """
	tre_file = arguments[ arguments.index('--tre')+1 ]
	aln_file = arguments[ arguments.index('--aln')+1 ]
	if '--black' in arguments:
		exclude_file = arguments[ arguments.index('--black')+1 ]
	else:
		exclude_file = ""
	if '--seq' in arguments:
		seq_file = arguments[ arguments.index('--seq')+1 ]
	else:
		seq_file = ""
	
	pattern = "[a-zA-Z]+[0-9a-zA-Z\-_\.@]*"	#"\w[0-9a-zA-Z\-_]+@\d+"
	
	if len( exclude_file ) > 0:
		ignore = get_names_to_exclude( exclude_file )
	else:
		ignore = []
	mask_monophyly( tre_file, aln_file, True, ignore )
	
	if len( seq_file ) > 0:
		seqs = load_sequences( seq_file, pattern )
		with open( tre_file + ".mm", "r" ) as f:
			content = f.read()
		IDs = re.findall( pattern, content )
		with open( seq_file + ".mm", "w" ) as out:
			for ID in IDs:
				try:
					out.write( '>' + ID + '\n' + seqs[ ID ] + "\n" )
				except KeyError:
					pass


if '--tre' in sys.argv and '--aln' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
