import argparse

parser = argparse.ArgumentParser(description= "Select Kimura under param and create bed")
parser.add_argument("--awd", "-align_with_div",  help="out of calcDivergenceFromAlign.pl")
parser.add_argument("--k", "-kimura",  help="Kimura threshold")
parser.add_argument("-out", "--o", help="output file",type=str)
parser.add_argument("-bed", "--b", help="output Bed file",type=str)
arg = parser.parse_args()
filehandle = open(arg.awd)
Tranver={}
name={}
MutationRates={}
Id="asd"
transitions_list=[]
for line in filehandle:
        if line.startswith(" "):
                continue
        elif line.startswith("\t"):
                continue
        elif line.startswith("Matrix"):
                continue
        elif line.startswith("Gap_init"):
                continue
        elif line.startswith("C"):
                continue
        elif line[0].isdigit():
            stripped=line.split(" ")
            if  stripped[8].startswith("C"):
                TE=stripped[9]
            else:
                TE=stripped[8]
            Start=int(stripped[5])
            End=int(stripped[6])
            Size=End-Start
            Id = (stripped[4]+"\t"+ stripped[5]+"\t"+ stripped[6]+"\t"+str(Size)+"\t"+TE)
        elif line.startswith("Kimura"):
                stripped=line.split("=")
                if float(stripped[1]) <= float(arg.k):
                    kimuraval = stripped[1].rstrip("\n")
                    name[Id]=kimuraval
        elif line.startswith("Transitions"):
                strippedt=line.split(" ")
                trans=strippedt[4]
                Tranver[Id]=trans
                Parentesis=strippedt[5]
                Transitions=Parentesis.split("/")[0].lstrip("(")
                Ph=Parentesis.split("/")[1]
                Transversions=Ph.strip(")\n")
                Mutations=int(Transitions)+int(Transversions)
                if Mutations==0:
                    MutationRate=0
                else:
                    MutationRate=str((int(Mutations)/int(Size)))
                MutationRates[Id]=str(MutationRate)
oWriter=open(arg.o, "w")
oWriter.write("Scaffold"+"\t"+"beg"+"\t"+"end"+"\t"+"Size"+"\t"+"type"+"\t"+"kimuraval"+"Transitions/transversions"+"\t"+"MutationRate"+"\n")
bWriter=open(arg.b, "w")
for rec in name:
        oWriter.write(rec+"\t" +name[rec]+"\t" + Tranver[rec]+"\t" +MutationRates[rec]+"\n")
        bWriter.write(rec+"\t" +name[rec]+"\t" + Tranver[rec]+"\t" +MutationRates[rec]+"\n")

