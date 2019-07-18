#Split the fasta file base on domain in to multiple file
#Replace "U" to "T", and "." to "-"
#Requrement Biopython
from Bio import SeqIO
f = open("test.fasta", "r")
Bacteria = open("Bacteria.fasta", "w")
Eukaryota = open("Eukaryota.fasta", "w")
Archaea = open("Archaea.fasta", "w")
Other = open("Other.fasta", "w")
for seq_record in SeqIO.parse(f,"fasta"):
    seq = str(seq_record.seq)
    reseq = seq.replace("U", "T").replace(".", "-")
    desc = seq_record.description.split(" ")
    domain = desc[1].split(";")
    if domain[0] == "Bacteria":
        Bacteria.write(">"+str(seq_record.description)+"\n"+str(reseq))
    elif domain[0] == "Eukaryota":
        Eukaryota.write(">"+str(seq_record.description)+"\n"+str(reseq))
    elif domain[0] == "Archaea":
        Eukaryota.write(">"+str(seq_record.description)+"\n"+str(reseq))
    else:
        Other.write(">"+str(seq_record.description)+"\n"+str(reseq))  
        
f.close()        
Bacteria.close()
Eukaryota.close()
Archaea.close()
Other.close()