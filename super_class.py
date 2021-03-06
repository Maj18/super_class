#!/anaconda3/bin/python3

import re

class Sequence:
	#super class for sequences. Contains basic sequences information.
	def __init__(self, sequence):
		self.seq = sequence

	def getSubSequenceCount(self, sequence, sub_sequence):
		#Takes a sequence and a subsequence and then returns the number of sub sequence occurences in the sequence.
		return len(re.findall(sub_sequence, sequence, re.I))

	def checkSeq(self, valid):
		valid_count = 0
		for item in valid.lower():
			valid_count += self.seq.lower().count(item)
		if valid_count != len(self.seq):
			#print("Error! The sequence does NOT only contain standard nucleotides or amino acids") 
			return 'ERROR'
		else:
			return ("All valid sequence!")

	def as_dnaORrna (self, standard, valid):
		if self.checkSeq(valid) == 'ERROR':
			print("")
		else: 
			count = 0
			for i in standard.lower():
				count += self.seq.lower().count(i)
			if count == len(self.seq):
				return self.seq

	def translation(self, seq, aaCODE):
		trans_result = []
		for n in range(3):
			fv_aa = ''
			while n+2 < len(seq):
				aa_code = seq[n] + seq[n+1] + seq[n+2]
				fv_aa += aaCODE[aa_code.upper()]
				n += 3 
			trans_result.append(fv_aa)
		return trans_result

	def __str__(self):
		return self.seq

class NuclSequence(Sequence):
	#Sub class to sequence class, representing a nucleotide sequence
	def __init__(self, sequence):
		super().__init__(sequence)

	def getRevComp(self, seq):
		#Retrieve the reverse complements for the sequence
		translation = seq.maketrans('atgc', 'tacg')
		return seq[::-1].translate(translation)

	def getOligoCount(self, oligo):
		#
		fw_count = super().getSubSequenceCount(self.seq, oligo)
		rv_seq = self.getRevComp()
		rv_count = super().getSubSequenceCount(rv_seq, oligo)
		print(super().checkSeq('ATGCU'))
		print("The DNA sequence: {}".format(super().as_dnaORrna('ATGC', 'ATGCU')))
		print("The RNA sequence: {}".format(super().as_dnaORrna('AUGC','ATGCU')))
		return fw_count + rv_count

	def Translateseq(self):
		code_dna =[(x+y+z) for x in ["T","C","A","G"] for y in ["T","C","A","G"] for z in ["T","C","A","G"]]
		code_rna =[(x+y+z) for x in ["U","C","A","G"] for y in ["U","C","A","G"] for z in ["U","C","A","G"]]
		aa=['F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*'\
			,'W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R'\
			,'R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R'\
			,'R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G'\
			,'G']
		aaCODE_dna=dict(zip(code_dna,aa))
		aaCODE_rna=dict(zip(code_rna,aa))

		dnaSeq = super().as_dnaORrna('ATGC', 'ATGCU')
		rnaSeq = super().as_dnaORrna('AUGC', 'AUGCU')
		#if (not dnaSeq):
			#print("The DNA seq can not be translated!")
		if dnaSeq:
			trans_fw = super().translation(dnaSeq, aaCODE_dna)
			rv_seq = self.getRevComp(dnaSeq)
			trans_rv = super().translation(rv_seq, aaCODE_dna)

			trans_dna = trans_fw + trans_rv
			print(trans_dna)

			for trans_version in trans_dna:
				if (trans_version.count('*') == 0) or ((trans_version.count('*') == 1) \
						and trans_version[len(trans_version)-1] == '*'):
					return trans_version

		#elif not rnaSeq:
			#print("The RNA seq can not be translated!")
		elif rnaSeq:
			trans_rna = super().translation(rnaSeq, aaCODE_rna)
			print(trans_rna)

			for trans_version_r in trans_rna:
				if trans_version_r.count('*') == 0 or ((trans_version_r.count('*') == 1) \
						and trans_version_r[len(trans_version_r)-1] == '*'):
					return trans_version_r


class AminoSequence(Sequence):
	#Sub class to sequence class representing an amino acid sequence
	def __init__(self, sequence):
		super().__init__(sequence)

	def getPeptideCount(self, peptide):
		super().checkSeq('ARNDCEQGHILKMFPSTWYV')
		return super().getSubSequenceCount(self.seq, peptide)




if __name__ == '__main__':
	Amino_seq = AminoSequence("DDAGTKYWCDDEE_")
	print(Amino_seq.getPeptideCount('DD'))
	nu_seq = NuclSequence("ATGCAAGGCCAAAATTGCGCTACGTCCGGATGCCGAATTGCACTTAATGCAAGGCCAAAATTGCGCTACGTCCGGATGCCGAATTGCACTTAATGCAAGGCCAAAATTGCGCTACGTCCGGATGCCGAATTGCACTTAATGCAAGGCCAAAATTGCGCTACGTCCGGATGCCGAATTGCACTTA")
	print(nu_seq.Translateseq())



