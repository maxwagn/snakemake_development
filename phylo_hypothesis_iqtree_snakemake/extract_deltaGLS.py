### get BUSCO sets into dictionary

BUSCO_sets = dict()
with open("../results/BUSCO_Partitions_anabantoid_concat.tsv", "r") as iqtree:
    for line in iqtree:
        if line.startswith("Subset"):
            pass
        else:
            line = line.rstrip().split("\t")
            BUSCO_sets[str(line[10])]=[int(line[8]),int(line[9])]

#print(BUSCO_sets)


### calculate deltaGLS
header = "BUSCO\tlogLTree1\tlogLTree2\tlogLTree3\tlogLTree4\tlogLTree5\tdeltaLTree1\tdeltaLTree2\tdeltaLTree3\tdeltaLTree4\tdeltaLTree4\n"


with open("deltaGLS_anabantoid.tsv", "w") as outfile:
    outfile.write(header)
    for k, v in BUSCO_sets.items():
        GLSvals = []
        deltaGLSvals = []
        with open("../results/hypothesis_tests_1000bs/anabantoid_phylo_tests.sitelh", "r") as siteLHfh:
            for line in siteLHfh:
                if line.startswith("5"):
                    pass
                else:
                    line = line.split()
                    #print([i for i in value_list[v[0]-1:v[1]]])
                    GLS = sum([float(i) for i in line[1:]][v[0]-1:v[1]])
                    GLSvals.append(GLS)

        for i in GLSvals:
            deltaGLSvals.append(GLSvals[0]-i)

        outfile.write(k + "\t" + "\t".join(str(i) for i in GLSvals) + "\t" + "\t".join(str(i) for i in deltaGLSvals) + "\n")


