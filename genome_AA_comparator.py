__author__ = 'Wietrack 2015'

import sys
from subprocess import call
import os
import os.path
from operator import attrgetter

class Obj(object):
    """__init__() functions as the class constructor"""
    def __init__(self, text, score):
        self.text = text
        self.score = score
    def __repr__(self):
            return repr((self.text, self.score))

def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring

sim = 50.0
cov = 70.0
eval = 1e-5
outputTable = "comparison-summary.tab"
rewriteBlast = False
cores = 1

#process arguments
inputFiles = []
annotFiles = {}
index = 1
inputFilesReading = False
annotationFilesReading = False

print "--------------------------------"
print "--- Genome comparator v 0.90 ---"
print "---   Tomas Vetrovsky 2015   ---"
print "--------------------------------"
print " "

for arg in sys.argv:
    ch = arg[0]
    if ch == '-':
        inputFilesReading = False
        annotationFilesReading = False
        if arg == "-s":
            sim = float(sys.argv[index])
            print "similarity treshold set to: "+sys.argv[index]
        if arg == "-c":
            cov = float(sys.argv[index])
            print "coverage treshold set to: "+sys.argv[index]
        if arg == "-e":
            eval = float(sys.argv[index])
            print "evalue treshold set to: "+sys.argv[index]
        if arg == "-o":
            outputTable = sys.argv[index]
            print "output table: "+sys.argv[index]
        if arg == "-R":
            rewriteBlast = True
            print "Blast out will be redone."
        if arg == "-n":
            cores = int(sys.argv[index])
            print "Number of threads set for blast: "+str(cores)
        if arg == "-i":
            inputFilesReading = True
        #annotation files in GenBank formate corespondint to input files...(same)
        #same name with ".ant" ending!!!
        if arg == "-a":
            annotationFilesReading = True
    else:
        if inputFilesReading:
            inputFiles.append(arg)
            print "input file added: "+arg
        if annotationFilesReading:
            base = rchop(arg, '.gbk')
            annotFiles[base] = arg
            print "anottation file added: "+arg+" base:("+base+")"
    index = index + 1

#checking annotation files
print " "
if len(annotFiles)>0:
    includeAnnot = True
    for f in inputFiles:
        if annotFiles.has_key(f):
            print "anotation file found: "+annotFiles[f]
        else:
            includeAnnot = False
            print "Error: Annotation files must have the same names as input files with '.gbk' as ending!! Ignoring annotation files. ("+f+")"
            for annotF in annotFiles:
                print annotFiles[annotF]
            print " "
            break
else:
    includeAnnot = False

if includeAnnot:
    print "Annotation files ok :)"

if len(inputFiles)>1:
    #renaming imput files
    renamedFiles = []
    outputFiles = []
    allGenes = "all_genes.fa"
    fall = open(allGenes, 'w')
    lengths = {}
    title = ""
    seq_len = 0

    for f in inputFiles:
        renamedFile = f+".renamed"
        renamedFiles.append(renamedFile)
        fp = open(renamedFile, 'w')
        for line in open(f):
            ch = line[0]
            if ch == '>':
                title = f+"@"+line[1:].strip()
                seq_len = 0
                lengths[title] = seq_len

                fp.write(">"+title+"\n")
                fall.write(">"+title+"\n")
            else:
                fp.write(line)
                fall.write(line)
                seq_len = seq_len + len(line.strip())
                if lengths.has_key(title):
                    lengths[title] = seq_len
        fp.close()

    fall.close()

    bashLine = "grep '>' -o "+allGenes+" | wc -l"
    print "sequences - total: "
    os.system(bashLine)
    print " "

    #generatin databases
    print "generating blast databases... "
    for f in renamedFiles:
        bashLine = "formatdb -i "+f+" -o F -p T"
        print "running... "+bashLine
        os.system(bashLine)

    #blastp
    print "comparing files... "
    for f in renamedFiles:
        fileName = "blastp."+f+".out"
        outputFiles.append(fileName)
        outExist = False
        if os.path.isfile(fileName) and os.access(fileName, os.R_OK):
            outExist = True
        if rewriteBlast:
            outExist = False
        if outExist:
            print "blast result already exists - using: "+fileName
        else:
            bashLine = "./ggsearch36 "+allGenes+" "+f+" -m 8 > "+fileName
            print "running... "+bashLine
            os.system(bashLine)

    #setup results table
    titles = {}
    evals = {}
    for title in lengths:
        #print "title: "+title+" seq_len: "+str(titles[title])
        titles[title] = {}
        evals[title] = {}
        for f in outputFiles:
            titles[title][f] = "no-hit"
            evals[title][f] = 1000.0
        #print title + " length " + str(lengths[title])
            #print title+"-"+f+" "+str(titles[title][f])

    #filling results table
    print "checking blast results... "
    for f in outputFiles:
        for line in open(f):
            values = line.split("\t")
            #coverage = (float(values[3])/float(lengths[values[0]]))*100
            if int(values[7])>int(values[6]):
                partLen = int(values[7])-(int(values[6])-1)
            else:
                partLen = int(values[6])-(int(values[7])-1)
            coverage = (float(partLen)/float(lengths[values[0]]))*100
            if float(values[2]) >= sim and coverage >= cov:
                #print line.strip()+" coverage: "+str(format(coverage, '.0f'))+" align.len: "+values[3]+" query.len: "+str(lengths[values[0]])+" subj.len: "+str(lengths[values[1]])
                #print "similarity "+values[2]+" is above the treshold "+str(sim)+": "+values[1]+" e-val: "+values[10]
                if titles.has_key(values[0]) and titles[values[0]].has_key(f):
                        if float(values[10])<evals[values[0]][f]:
                            evals[values[0]][f] = float(values[10])
                            titles[values[0]][f] = values[1]+"/"+"s:"+str(format(float(values[2]), '.0f'))+"/c:"+str(format(coverage, '.0f'))
                else:
                    print "error title not found: "+values[0]+"-"+f

    #clean evals
    evals.clear()

    #sorting data
    data_obj = []
    for title in titles:
        if titles.has_key(title):
            line = title
            score = 0
            for f in outputFiles:
                if titles[title].has_key(f):
                    if titles[title][f] != "no-hit":
                        score = score +1
                    line = line + "\t" + str(titles[title][f])
                else:
                    print "error file name not found: "+title+"-"+f
            data_obj.append(Obj(line+"\t"+str(score)+"\n",score))
        else:
            print "error title not found: "+title
    #data_obj = sorted(data_obj, key=lambda obj: obj.score)   # sort by score
    data_obj = sorted(data_obj, key=attrgetter('score'), reverse=True)   # sort by score

    #saving output
    fo = open(outputTable, 'w')
    line = "genome@gene"
    for f in outputFiles:
        line = line + "\t" + f
    fo.write(line+"\tscore"+"\n")
    for obj in data_obj:
        fo.write(obj.text)
    fo.close()

    #cleaning lists...
    data_obj = []
    titles.clear()

    #annotation files extraction
    annotations = {}
    if includeAnnot:
        print "extracting annotation files..."
        for f in annotFiles:
            cdsRead = False
            done = False
            ref = ""
            prod = ""
            i = 0
            prodRead = False
            print "checking... "+annotFiles[f]
            for line in open(annotFiles[f]):
                if prodRead:
                    pos = line.find('"')
                    if pos > -1:
                        prod = prod + ' '.join(line.split())[:pos-len(line)+1]
                        prodRead = False
                        cdsRead = False
                    else:
                        prod = prod + ' '.join(line.strip().split())

                #CHECKING
                if 'SEED:' in line:
                    if i>0:
                        #print f+"@"+ref +" prod: "+prod
                        #print prod
                        annotations[f+"@"+ref] = prod
                        #print annotations[f+"@"+ref]
                    ref = ""
                    prod = ""
                    cdsRead = True
                    prodRead = False
                    i = i + 1
                    pos = line.find('SEED:')
                    if pos > -1:
                        #refRead = True
                        ref = line[pos+5:].strip()
                        #print "whole:"+ref
                        pos = ref.find('"')
                        if pos > -1:
                            #print "pos: "+str(pos)
                            ref = ref[:pos-len(ref)]
                        #print "part:"+ref
                        #refRead = False
                if cdsRead:
                    pos = line.find('/product="')
                    if pos > -1:
                        prodRead = False
                        prod = line[pos+10:].strip()
                        #print "whole:"+prod
                        pos = prod.find('"')
                        if pos > -1:
                            prod = prod[:pos-len(prod)]
                            cdsRead = False
                        else:
                            prodRead = True
            annotations[f+"@"+ref] = prod

    #annotate the out table
    print "annotate the out table... "
    fo = open(outputTable+".annot.txt", 'w')
    i = 0
    rem = 0
    for line in open(outputTable):
        if i == 0:
            fo.write("annotation\t"+line)
        else:
            values = line.strip().split("\t")
            gene_genome_line = values[0]
            genome_line = gene_genome_line.split("@")[0]
            index = 1
            annot = "-"
            if includeAnnot:
                if annotations.has_key(gene_genome_line):
                    annot = annotations[gene_genome_line]
                else:
                    print "WARNING: missing anotation for "+gene_genome_line
            new_line = annot+"\t"+gene_genome_line
            for f in outputFiles:
                gene_genome_val = values[index].split("/")[0]
                genome_val = gene_genome_val.split("@")[0]
                if genome_val != "no-hit":
                    new_line = new_line +"\t"+ values[index].split("@")[1]
                else:
                    new_line = new_line +"\t"+ "no-hit"
                index = index + 1
            fo.write(new_line+"\t"+values[index]+"\n")
        i = i + 1
    fo.close()

    #cleaning the out table
    print "cleaning the out table... "
    fo = open(outputTable+".cleaned.txt", 'w')
    to_remove = {}
    i = 0
    rem = 0
    passed = 0
    for line in open(outputTable):
        if i == 0:
            fo.write("annotation\t"+line)
        else:
            values = line.strip().split("\t")
            gene_genome_line = values[0]

            if to_remove.has_key(gene_genome_line):
                rem = rem + 1
            else:
                genome_line = gene_genome_line.split("@")[0]
                index = 1
                hits = []
                toRemove = False

                annot = "-"
                if includeAnnot:
                    if annotations.has_key(gene_genome_line):
                        #print "annotation: "+annotations[gene_genome_line]
                        annot = annotations[gene_genome_line]
                    else:
                        print "WARNING: missing anotation for "+gene_genome_line

                new_line = annot+"\t"+gene_genome_line
                for f in outputFiles:
                    gene_genome_val = values[index].split("/")[0]
                    genome_val = gene_genome_val.split("@")[0]
                    if genome_val != "no-hit":
                        if genome_line != genome_val:
                            to_remove[gene_genome_val] = 1
                        #print values[index]
                        new_line = new_line +"\t"+ values[index].split("@")[1]
                    else:
                        new_line = new_line +"\t"+ "no-hit"
                    index = index + 1

                fo.write(new_line+"\t"+values[index]+"\n")
                passed = passed + 1
        i = i + 1
    fo.close()

    print str(rem)+" lines were removed from cleaned file. "+str(passed)+" lines passed..."


    print "done :)"
else:
    print "not enough arguments... "