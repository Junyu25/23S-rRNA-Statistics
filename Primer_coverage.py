# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 19:51:23 2019

@author: Junyu
"""
import os
import pandas as pd
from Bio import SeqIO   
Arc = "~/Lab/23S-rRNA/Rawdata/Archaea.fasta"
Bac = "~/Lab/23S-rRNA/Rawdata/Bacteria.fasta"
Euk = "~/Lab/23S-rRNA/Rawdata/Eukaryota.fasta"
Fa = "~/Lab/23S-rRNA/Rawdata/SILVA_132_LSURef.fasta"
Full_ID = "~/Lab/23S-rRNA/Rawdata/SILVA_132_LSURef_full_ID.txt"
Taxa = "~/Lab/23S-rRNA/Rawdata/SILVA_132_LSURef.tax"

Primer = "/home/junyuchen/Lab/23S-rRNA/Primer/PrimerBac.fasta"
#Primer = "/home/junyuchen/Lab/23S-rRNA/Rawdata/Primer.fasta"
'''
for Seq in SeqIO.parse(Primer,"fasta"):
    analyze_primers1 = "analyze_primers.py -f "+str(Arc)+" -p "+str(Seq.id)+" -s "+str(Seq.seq)
    analyze_primers2 = "analyze_primers.py -f "+str(Bac)+" -p "+str(Seq.id)+" -s "+str(Seq.seq)
    analyze_primers3 = "analyze_primers.py -f "+str(Euk)+" -p "+str(Seq.id)+" -s "+str(Seq.seq)
    analyze_primers = "analyze_primers.py -f "+str(Fa)+" -p "+str(Seq.id)+" -s "+str(Seq.seq)

    #print(Seq.id)
    print(analyze_primers)
    os.system(analyze_primers1)
    os.system(analyze_primers2)
    os.system(analyze_primers3)
    os.system(analyze_primers)
'''
   
f = pd.DataFrame()

for Seq in SeqIO.parse(Primer,"fasta"):
    hits1 = Seq.id+"_Archaea"+"_hits.txt"
    hits2 = Seq.id+"_Bacteria"+"_hits.txt"
    hits3 = Seq.id+"_Eukaryota"+"_hits.txt"
    hits = Seq.id+"_SILVA_132_LSURef"+"_hits.txt"

    get_amplicons1 = "get_amplicons_and_reads.py -f "+str(Arc)+" -i "+str(hits1)+" -R 250 -o ./Arc"
    get_amplicons2 = "get_amplicons_and_reads.py -f "+str(Bac)+" -i "+str(hits2)+" -R 250 -o ./Bac"
    get_amplicons3 = "get_amplicons_and_reads.py -f "+str(Euk)+" -i "+str(hits3)+" -R 250 -o ./Euk"
    get_amplicons = "get_amplicons_and_reads.py -f "+str(Fa)+" -i "+str(hits)+" -R 250 -o ./Fa"
    print(hits)
    os.system(get_amplicons1)
    os.system(get_amplicons2)
    os.system(get_amplicons3)
    os.system(get_amplicons)
    
    a1 = SeqIO.index("./Arc/"+str(Seq.id)+"_amplicons.fasta", "fasta")
    a2 = SeqIO.index("./Bac/"+str(Seq.id)+"_amplicons.fasta", "fasta")
    a3 = SeqIO.index("./Euk/"+str(Seq.id)+"_amplicons.fasta", "fasta")
    a = SeqIO.index("./Fa/"+str(Seq.id)+"_amplicons.fasta", "fasta")
    
    f = f.append({'Arc':len(a1), "Bac":len(a2), 'Euk':len(a3), 'Fa':len(a)}, ignore_index=True)
    #print(get_amplicons)
    
    #print(wc)
f.to_csv("PrimerBac.csv", encoding = "utf-8")