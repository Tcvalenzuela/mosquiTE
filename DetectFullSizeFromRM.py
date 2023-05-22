import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description= "Select the insertions that are 90% of the lenght of consensus and k divergence from it")
parser.add_argument("--awd", "-align_with_div",  help="out of calcDivergenceFromAlign.pl")
parser.add_argument("--k", "-kimura",  help="Kimura threshold")
parser.add_argument("--S", "-SeqConse",  help="Sequence consensus")
arg = parser.parse_args()
filehandle = open(arg.awd)
Tranver={}
InsertionSize={}
name={}
ListSeq={}
Count={}
Id="asd"
transitions_list=[]
consensus = SeqIO.parse(open(arg.S),"fasta")
for fasta in consensus:
    seq_name = fasta.id
    seq_leng = len(fasta.seq)
    ListSeq[seq_name] = str(seq_leng)

for line in filehandle:
        if line.startswith(" ") or line.startswith("\t") or line.startswith("Matrix") or line.startswith("Gap_init") or line.startswith("C"): 
                continue
        elif line[0].isdigit():
            stripped=line.split(" ")
            if  stripped[8].startswith("C"):
                if stripped[9].split("#")[1] != "Unspecified":
                    continue
                TE=stripped[9].split("#")[0]
            else:
                if stripped[8].split("#")[1] != "Unspecified":
                    continue
                TE=stripped[8].split("#")[0]

            Start=int(stripped[5])
            End=int(stripped[6])
            Size=End-Start
            ScaffInsertion=stripped[4]
            Id = (stripped[4]+"\t"+ stripped[5]+"\t"+ stripped[6]+"\t"+str(Size)+"\t"+TE)
        elif line.startswith("Kimura"):
                stripped=line.split("=")
                name[Id]=(Id+"\t"+stripped[1])

            
for ele in name:
    #name[Id]=float(kimuraval)
    Splitter2=name[ele].split("\t")
    Scaff=Splitter2[0]
    start2=Splitter2[1]
    end2=Splitter2[2]
    size2=Splitter2[3]
    SeqID=Splitter2[4]
    Kimura=Splitter2[5]
    if float(Kimura) <= float(arg.k) and float(size2) > 0.9*float(ListSeq[SeqID]):
#        print (name[ele]+"\t"+str(ListSeq[SeqID])+"\t"+str(size2))
        if SeqID not in Count:
            Count[SeqID] = 1
        else:
            Count[SeqID] += 1
 


for ele in Count:
    print (ele+"\t"+str(Count[ele]))
