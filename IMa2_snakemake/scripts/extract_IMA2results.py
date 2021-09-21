import sys

run1_file = sys.argv[1]
run2_file = sys.argv[2]
run3_file = sys.argv[3]
outfh = sys.argv[4]

#run1_file = "../results/Umdloti_Umkomaas/1mio_Umdloti_Umkomaas_run1_IMa2.txt"
#run2_file =  "../results/Umdloti_Umkomaas/1mio_Umdloti_Umkomaas_run2_IMa2.txt"
#run3_file =  "../results/Umdloti_Umkomaas/1mio_Umdloti_Umkomaas_run3_IMa2.txt"
#outfh = "../results/Umdloti_Umkomaas/1mio_Umdloti_Umkomaas_IMa2.tsv"

def extract_pop_names(fh):
    with open(fh, "r") as infile:
        pop0 = str()
        pop1 = str()
        for line in infile:
            # line = line.rstrip()
            if line.startswith("Population 0 :"):
                pop0 += str(line.split(": ")[1])
            elif line.startswith("Population 1 :"):
                pop1 += str(line.split(": ")[1])
        return("## Population 0: {}## Population 1: {}".format(pop0, pop1))


def extract_migration_table(fh, run):
    with open(fh, "r") as infile:
        r = run
        pop0 = str()
        pop1 = str()
        copy = False
        outstring = str()
        #start="\tHiPt"
        #end = " SumP"
        start = "HISTOGRAM GROUP 2: MARGINAL DISTRIBUTION VALUES AND HISTOGRAMS OF POPULATION SIZE AND MIGRATION PARAMETERS"
        end = "HISTOGRAM GROUP 3: POPULATION MIGRATION (2NM) POSTERIOR PROBABILITY HISTOGRAMS"
        text_block = str()
        for line in infile:
            #line = line.rstrip()
            if line.startswith("Population 0 :"):
                pop0 += str(line.split(": ")[1])
            elif line.startswith("Population 1 :"):
                pop1 += str(line.split(": ")[1])

            if line.strip().startswith(start):
                copy = True
                continue
            elif line.strip().startswith(end):
                copy = False
                continue
            elif copy:
                outstring += line

        outlist = outstring.split("\n")

        outstring2 = str()

        prefix = "\t   "

        header = list()


        return_string = str()

        for i in outlist:
            if i.startswith("	Parameter"):
                header = str(i)[1:].split("\t")
                #header = header.split("\t")

            for j in range(0,999):
                if i.startswith("{}{}".format(prefix,j)):
                    tes = i.split("\t")[1:]
                    return_string += tes[0] + "\t" + header[-4] + "\t" + tes[-4] + "\t" + tes[-3] + "\t" + r + "\n" + \
                          tes[0] + "\t" + header[-2] + "\t" + tes[-2] + "\t" + tes[-1] + "\t" + r + "\n"
                elif i.startswith("{}{}".format(prefix[:-len(str(j))],j)):
                    tes = i.split("\t")[1:]
                    #test = i
                    return_string += tes[0] + "\t" + header[-4] + "\t" + tes[-4] + "\t" + tes[-3] + "\t" + r + "\n" + \
                                     tes[0] + "\t" + header[-2] + "\t" + tes[-2] + "\t" + tes[-1] + "\t" + r + "\n"
                    #print(i.split("\t")[1:])


        return(return_string)




            #pops0 = "\t".join(i.split("\t")[2:4]) + "\t" + "2N0m0>1_{}_into_{}".format(pop0, pop1) + "\t" + r
            #pops0 = pops0.replace(" ", "")
            #pops0 = pops0.replace("\n", "")
            #pops1 = "\t".join(i.split("\t")[4:6]) + "\t" + "2N1m1>0_{}_into_{}".format(pop1, pop0) + "\t" + r
            #pops1 = pops1.replace(" ", "")
            #pops1 = pops1.replace("\n", "")

            #outstring2 += pops0 + "\n" + pops1 + "\n"

        #return(i)

preheader="""## curve height is an estimate of the posterior probability
## each term is the product of a population parameter (e.g. q0) and a migration rate (e.g.m0>1) 
## migration rates are in the coalescent (backwards in times), so that a population migration rate of 
## q1m0>1 is the population rate (forward in time) at which population 1 receives migrants from population 0\n"""

popnames = extract_pop_names(run1_file)

header = "parameter\tdirection\tmigration\tposterior\trun\n"

run1 = extract_migration_table(run1_file, "run1")
run2 = extract_migration_table(run2_file, "run2")
run3 = extract_migration_table(run3_file, "run3")

outfile = preheader + popnames + header + run1 + run2 + run3

with open(outfh, "w") as ouf:
    ouf.write(outfile)


