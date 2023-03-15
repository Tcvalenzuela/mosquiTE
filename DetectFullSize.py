from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description= "Select Kimura under param and create bed")
parser.add_argument("--S", "-SeqConse",  help="Sequence consensus")
parser.add_argument("--T", "-tableInsertion",  help="Out of RepeatMasker table")
arg = parser.parse_args()
ListSeq={}
Insertions={}
Divergence={}
ListNames=[]
consensus = SeqIO.parse(open(arg.S),"fasta")
TEorganizer = {}
for fasta in consensus:
	name = fasta.id
	seq = fasta.seq
	leng = len(seq)
	ListSeq[name] = str(leng)
	ListNames.append(name)

Reps = open(arg.T)
for insertion in Reps:
	if insertion.startswith("Scaff"):
		continue
	stripped=insertion.split("\t")
	Scaff=stripped[0]
	start=stripped[1]
	end=stripped[2]
	size=stripped[3]
	TEtype=stripped[4].split("#")[0]
	Kval=stripped[5]
	if float(ListSeq[TEtype]) * 0.9 <=  float(size) <= float(ListSeq[TEtype])+(float(ListSeq[TEtype]) * 0.1):
		if float (Kval) < 20:
			#print (str(ListSeq[TEtype])+"\t"+insertion.rstrip("\n"))
			TEorganizer[insertion] = TEtype



for name in ListNames:
	value={i for i in TEorganizer if TEorganizer[i]==name}
	for ele in value:
		print (ele.rstrip("\n"))
