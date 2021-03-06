#!/anaconda3/bin/python3

import sys
import argparse
import seq_class as sc

#check whether the user get the command line right
if len(sys.argv) != 5:
	sys.exit("ERROR: please run the program as follows: ./use_seq_class.py -i INFILENAME -G OLIGO")

#parse the command arguments and options
parser = argparse.ArgumentParser(description = "Thi program will count \
	the number of the provided oligo in the sequences from the provided file")
parser.add_argument(
	'-g',
	type = str,
	metavar = 'OLIGO',
	dest = 'oligo',
	required = True
	)
parser.add_argument(
	'-i',
	type = argparse.FileType('r'),
	metavar = 'INFILE',
	dest = 'infile',
	required = True
	)
args = parser.parse_args()
oligo = args.oligo
infile = args.infile

count = 0
for line in infile:
	line = line.rstrip()
	if not line:
		continue
	elif not line.startswith('>'):
		#create a sequence object
		seq_obj = sc.NuclSequence(line)
		#count the sequence object and add it up
		count += seq_obj.getOligoCount(oligo)


print('Oligo occurences: {}'.format(count))


target_oligo = sys.argv