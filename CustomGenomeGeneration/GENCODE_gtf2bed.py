HELP_STRING = """
GENCODE_gtf2bed.py

Author: Mary Couvillion
Date: 5/7/18, updated 2/2021


     -h     print this help message
     -i     file input: gtf (required)
     -o     output file name: *.bed (required)
     -g     genome (human or mouse)

"""
import sys
from getopt import getopt


def main(argv=None):
    if argv is None:
        argv = sys.argv

    inputFile = ""
    outputFile = ""
    genomeType = ""

    try:
        optlist, args = getopt(argv[1:], "hi:o:g:")
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
        elif opt == '-g':
            genomeType = opt_arg
    if inputFile == "" or outputFile == "" or genomeType == "":
        print HELP_STRING
        sys.exit(1)
    

    
    InFile = file(inputFile,'r')
    outFile = open(outputFile, 'w')


    # Set fields for different gtf file formats
    if genomeType == 'human':
        gene_name_field = 4
        gene_biotype_field = 6
    elif genomeType == 'mouse':
        gene_name_field = 4
        gene_biotype_field = 6

    # iterate through and process all lines in input file

    for i,line in enumerate(open(inputFile)):
        if i%100000 == 0:
            print i
        if line[0] != '#':
            clist = line.replace('\n','').split('\t')
            # process transcript entries
            entry_type = clist[2]
            if line[0] != 'd' and line[0] != 'e' and entry_type == 'transcript':
                chr = clist[0]
                gene_biotype = clist[8].split(';')[gene_biotype_field].split('"')[1]
                gene_name = clist[8].split(';')[gene_name_field].split('"')[1]
                strand = clist[6]
                start = clist[3]
                end = clist[4]
                length = str(int(end) - int(start) + 1)
                if gene_biotype == 'protein_coding' or gene_biotype == 'snRNA' or gene_biotype == 'rRNA' or gene_biotype == 'tRNA':
                    outFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chr, start, end, gene_name, '0', strand))
                else: 
                    continue
#                 outFile.write('%s\t%s\t%s\t%s\t%s%s%s\n' % ('MT', gene_name, strand, length, clist[3], '-', clist[4]))

#                exon = int(clist[4])-int(clist[3])+1
#                exon_lengths.append(exon)
##                    exon_ranges.append(str(clist[3])+'-'+ str(clist[4]))
#            elif clist[8].split(';')[2].split('"')[1] == 'protein_coding' and clist[2] == 'stop_codon':
#                ORF_length = sum(exon_lengths)
#                exon_lengths=[]
#                exon_ranges=[]
#

            else:
                continue
        else:
            continue
            
    outFile.close()



##############################################
if __name__ == "__main__":
    sys.exit(main())
