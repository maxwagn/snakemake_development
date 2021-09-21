import sys

INPUT = sys.argv[1]
POPCOMB = sys.argv[2]

POP1 = str(POPCOMB.split("_")[0])
POP2 = str(POPCOMB.split("_")[1])

def read_fasta(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))



with open(INPUT, "r") as infasta, open("Input_files/{}_{}.IMa2.txt".format(POP1, POP2), "w") as outfile:

    pop1str = str()
    pop2str = str()
    noPOP1 = 0
    noPOP2 = 0
    seqlength = int()

    for seqID, sequence in read_fasta(infasta):
        if seqID.split("_")[-1] == POP1:
            seqID = "_".join(seqID.split("_")[0:-1])
            seqID = seqID.replace(">","")
            idlen = len(seqID)
            diflen = 10 - idlen
            pop1str += str(seqID + ' ' * diflen + sequence + "\n")
            noPOP1 += 1
            seqlength = len(sequence)
        elif seqID.split("_")[-1] == POP2:
            seqID = "_".join(seqID.split("_")[0:-1])
            seqID = seqID.replace(">","")
            idlen = len(seqID)
            diflen = 10 - idlen
            pop2str += str(seqID + ' ' * diflen + sequence + "\n")
            noPOP2 += 1
        else:
            pass


    text = """{pop1} vs {pop2}
#IMA2 input file for {pop1} vs {pop2}
2
{pop1} {pop2}
(0,1):2
1
COI {lenPOP1} {lenPOP2} {seqlen} H 0.25
{alignpop1}{alignpop2} """.format(pop1=POP1, pop2=POP2, lenPOP1=noPOP1 ,lenPOP2=noPOP2,
               seqlen=seqlength, alignpop1=pop1str, alignpop2=pop2str)


    outfile.write(text)
    #for k,v in POP1dict.items()


