#!/usr/bin/python
from __future__ import print_function
import re,sys,os

#This script is for checking whether the lenght of seqs in the aligned database are the same

Efa = open('database_len_summary.txt', 'w')
# Etax = open('Eukaryota.tax', 'w')
# Bfa = open('Bacteria.fa', 'w')
# Btax = open('Bacteria.tax', 'w')
# Afa = open('Archaea.fa', 'w')
# Atax = open('Archaea.tax', 'w')

with open(sys.argv[1]) as infile: #Open the file $fa_f at the path given by sys.argv[1], and assign the name infile to this file object
	for i in infile:
		if re.match('>', i): #match the begining of the fasta
			isplit = str(i.rstrip()) #Strip the firt line of the fasta
		else:
			str_length = str(len(i.rstrip())) #count the length of striped sequence
			#print str_length
			print (isplit+" "+str_length, file = Efa)
			

Efa.close()
# Etax.close()
# Bfa.close()
# Btax.close()
# Afa.close()
# Atax.close()

#print(Ecounter)
#print(Bcounter)
#print(Acounter)