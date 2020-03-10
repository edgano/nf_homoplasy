import operator
import os.path
import sys

fil=sys.argv[1]


class Fam:
    def __init__(self,fam,aligner,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint):
        self.fam=fam
        self.aligner=aligner
        self.Nratio=Nratio
        self.Nnumber=Nnumber
        self.usedNumber=usedNumber
        self.unusedRatio=unusedRatio
        self.unusedNumber=unusedNumber
        self.totalpoint=totalpoint

class Aln:
    def __init__(self,aligner,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,nfam,minacc,maxacc,deltaacc):
        self.aligner=aligner
        self.Nratio=Nratio
        self.Nnumber=Nnumber
        self.usedNumber=usedNumber
        self.unusedRatio=unusedRatio
        self.unusedNumber=unusedNumber
        self.totalpoint=totalpoint
        self.nfam=nfam
        self.minacc=minacc
        self.maxacc=maxacc
        self.deltaacc=deltaacc

class Global:
    def __init__(self,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,minacc,maxacc,deltaacc):
        self.Nratio=Nratio
        self.Nnumber=Nnumber
        self.usedNumber=usedNumber
        self.unusedRatio=unusedRatio
        self.unusedNumber=unusedNumber
        self.totalpoint=totalpoint
        self.minacc=minacc
        self.maxacc=maxacc
        self.deltaacc=deltaacc


fams=[]
alns=[]
with open(fil) as f:
    for line in f:
        line=line.strip("\n")

        if line[0:3]=="FAM":
            fields=line.split(" ")
            
            # Get family name
            tmp=fields[0].split(":")
            tmp=[x for x in tmp if x]   
            fam=tmp[1]
            
            # Get aligner name
            aligner=fields[1][1:-2]

            # Get N / total used ratio
            Nratio=float(fields[2])

            # Get number of N
            tmp=fields[3].split("/")
            Nnumber=int(tmp[0][1:])
            usedNumber=int(tmp[1][:-1])

            # Get unused ratio
            unusedRatio=float(fields[6])
            tmp=fields[7].split("/")
            unusedNumber=int(tmp[0][1:])
            totalpoint=int(tmp[1][:-1])

            tmp=Fam(fam,aligner,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint)
            fams.append(tmp)

        elif line[0:3]=="ALN":
            fields=line.split(" ")

            # Get aligner name
            tmp=fields[0].split(":")
            tmp=[x for x in tmp if x]   
            aligner=tmp[1]

            # Get N / total used ratio
            Nratio=float(fields[2])

            # Get number of N
            tmp=fields[3].split("/")
            Nnumber=int(tmp[0][1:])
            usedNumber=int(tmp[1][:-1])

            # Get unused ratio
            unusedRatio=float(fields[6])
            tmp=fields[7].split("/")
            unusedNumber=int(tmp[0][1:])
            totalpoint=int(tmp[1][:-1])

            # Get number of used families
            nfam=int(fields[8][1:])

            # TC or SP score on average if the worst/best MSA is picked every time
            minacc=float(fields[11])
            maxacc=float(fields[14])
            deltaacc=float(fields[-1][:-1])

            tmp=Aln(aligner,Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,nfam,minacc,maxacc,deltaacc)
            alns.append(tmp)

        elif line[:6]=="GLOBAL":
            fields=line.split(" ")
            Nratio=float(fields[1])
            tmp=fields[2].split("/")
            Nnumber=int(tmp[0][1:])
            usedNumber=int(tmp[1][:-1])
            unusedRatio=float(fields[5])
            tmp=fields[6].split("/")
            unusedNumber=int(tmp[0][1:])
            totalpoint=int(tmp[1][:-1])
            minacc=float(fields[9])
            maxacc=float(fields[12])
            deltaacc=float(fields[15][:-1])
            glob=Global(Nratio,Nnumber,usedNumber,unusedRatio,unusedNumber,totalpoint,minacc,maxacc,deltaacc)


# Write output
pre=os.path.splitext(fil)[0]
o=open(pre+"_fam.tsv","w")
o.write("family\taligner\tNratio\tNnumber\tusedNumber\tunusedRatio\tunusedNumber\ttotal\n")
for c in sorted(fams, key=lambda x: (x.aligner, x.fam)):
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(c.fam, c.aligner, c.Nratio, c.Nnumber, c.usedNumber, c.unusedRatio, c.unusedNumber, c.totalpoint))
o.close()
o=open(pre+"_aln.tsv","w")
o.write("aligner\tNratio\tunusedRatio\tnfam\tminacc\tmaxacc\tdeltaacc\n")
for c in alns:
    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(c.aligner, c.Nratio, c.unusedRatio, c.nfam, c.minacc, c.maxacc, c.deltaacc))
o.close()
o=open(pre+"_global.tsv","w")
o.write("Nratio\tunusedRatio\tminacc\tmaxacc\tdeltaacc\n")
o.write("{}\t{}\t{}\t{}\t{}\n".format(glob.Nratio, glob.unusedRatio, glob.minacc, glob.maxacc, glob.deltaacc))
o.close()