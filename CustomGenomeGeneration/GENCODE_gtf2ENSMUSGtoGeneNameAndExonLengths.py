HELP_STRING = """
GENCODE_gtf2ENStoGeneNameAndExonLengths.py

Author: Mary Couvillion
Date: 3/13/21


     -h     print this help message
     -i     file input (required)
     -o    ENSMUSTGoGeneName output file name (required)
     -O    GeneName and exon length output file name (required)
"""
import sys
from getopt import getopt


def main(argv=None):
    if argv is None:
        argv = sys.argv

    inputFile = ""
    outputFile = ""
    outputFile2 = ""

    try:
        optlist, args = getopt(argv[1:], "hi:o:O:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        #print opt
        #print opt_arg
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-i':
            inputFile = opt_arg
        elif opt == '-o':
            outputFile = opt_arg
        elif opt == '-O':
            outputFile2 = opt_arg
    if inputFile == "" or outputFile == "" or outputFile2 == "":
        print HELP_STRING
        sys.exit(1)
    

    
    InFile = file(inputFile,'r')
    outFile = open(outputFile, 'w')
    outFile2 = open(outputFile2, 'w')
    
    # iterate through and process all lines in input file
    outFile.write('GeneID\tGeneName\tgene_biotype\n')
    outFile2.write('GeneName\ttxptID\texonLength\n')

    for i,line in enumerate(open(inputFile)):
        if i%100000 == 0:
            print i
        if line[0] != '#':
            clist = line.replace('\n','').split('\t')
            # process transcript entries
            if line[0] == 'm' and clist[2] == 'gene':
                GeneID = clist[8].split(';')[0].split('"')[1]
                GeneName = clist[8].split(';')[2].split('"')[1]
                gene_biotype = clist[8].split(';')[4].split('"')[1]
                outFile.write('%s\t%s\t%s\n' % (GeneID, GeneName, gene_biotype))
            if line[0] == 'm' and clist[2] == 'exon':
                GeneName = clist[8].split(';')[5].split('"')[1]
                txptID = clist[8].split(';')[2].split('"')[1]
                exonLength = int(clist[4]) - int(clist[3]) + 1
                outFile2.write('%s\t%s\t%s\n' % (GeneName, txptID, str(exonLength)))
            elif line[0] == 'e' and clist[2] == 'exon':
                GeneID = clist[8].split(';')[0].split('"')[1]
                GeneName = clist[8].split(';')[0].split('"')[3]
                gene_biotype = 'ERCC'
                txptID = GeneName
                exonLength = int(clist[4]) - int(clist[3]) + 1
                outFile.write('%s\t%s\t%s\n' % (GeneID, GeneName, gene_biotype))
                outFile2.write('%s\t%s\t%s\n' % (GeneName, txptID, str(exonLength)))
            else:
                continue
        else:
            continue
            
    InFile.close()
    outFile.close()
    outFile2.close()


##############################################
if __name__ == "__main__":
    sys.exit(main())
