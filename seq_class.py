#!/anaconda3/bin/python3
'''
Description:
	This program will get the reverse complements for a sequence.

Usage:
	./seq_class
	
'''
import re

class NuclSequence:

	def __init__(self, sequence):
		self.seq = sequence

	def getRevComp(self):
		rev_low_seq = self.seq[::-1].lower()
		#make a translation table:
		translation = rev_low_seq.maketrans('atgc', 'tacg')
		#return the reverse complements
		return rev_low_seq.translate(translation)

	def getOligoCount(self, oligo):
		#count the number of oligo in the forward seqeuence
		fw_count = len(re.findall(oligo, self.seq, re.I))
		rev_seq = self.getRevComp()
		#count the number of oligo in the reverse sequence
		rev_count = len(re.findall(oligo, rev_seq, re.I))
		return fw_count + rev_count

	def __str__(self):
		return self.seq

#The codes below only be executed when you run the script directly:
#python myscript.py
#If you import this script to another model, they will not be executed.
if '__name__' == '__main__':
	seq_object = NuclSequence('ACGTATAGCTAG')
	print(seq_object.getRevComp())
	print(seq_object)
	print(seq_object.getOligoCount('Acg'))