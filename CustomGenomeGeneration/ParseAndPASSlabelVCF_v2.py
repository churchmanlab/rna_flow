# Written 1/28/21 to add 'PASS' to filter column of vcf over some quality threshold
# v2 added 8/23/23, adds a filtering step for removing very long indels and calls too close together

import sys

infilename = sys.argv[1]
QUALcutoff = int(sys.argv[2])
AFcutoff = float(sys.argv[3])

# J1-tot_snpcalls.vcf
invcf = open(infilename, 'r')
outvcf = open(infilename.replace('.vcf', '.QUALfilt' + str(QUALcutoff) + '.vcf'), 'w')

lastpos=0
lastchr=1
for line in invcf:
	if line[0] == '#':
		outvcf.write(line)
	else:
		col = line.replace('\n', '').split('\t')
		chrom = col[0]
		pos = int(col[1])
		id = col[2]
		ref = col[3]
		alt = col[4]
		qual = float(col[5])
		filter = col[6]
		info = col[7]
		format = col[8]
		s1 = col[9]
		s2 = col[10]
		fieldslist = info.split(';')
		# Get values from the DP4 field (first 2 are REF, second 2 are ALT)
		subs='DP4'
		DP4 = str([i for i in fieldslist if subs in i])
		DP4 = DP4.split("'")[1]
		DP4 = DP4.split('=')[1].split(',')
		DP4 = [float(i) for i in DP4] # Convert list elements to float
		# Calculate allele frequency
		tot = sum(DP4)
		AF = (DP4[2] + DP4[3])/tot

		# Add pass filter based on quality score
		if qual > QUALcutoff and tot > 4:
			if chrom == lastchr:
				if pos > (lastpos + 51) and len(alt) <= 50 and len(ref) <= 50 and AF > AFcutoff:
					outvcf.write(chrom +'\t'+ str(pos) +'\t'+ id +'\t'+ ref+'\t'+ alt+'\t'+ str(qual) +'\t'+ 'PASS' +'\t'+ info +'\t'+ format +'\t'+ s1 +'\t'+ s2 + '\n')
					lastpos = pos
					lastchr = chrom
				else:
					continue
			else:
				if len(alt) <= 50 and len(ref) <= 50 and AF > AFcutoff:
					outvcf.write(chrom +'\t'+ str(pos) +'\t'+ id +'\t'+ ref+'\t'+ alt+'\t'+ str(qual) +'\t'+ 'PASS' +'\t'+ info +'\t'+ format +'\t'+ s1 +'\t'+ s2 + '\n')
					lastpos = pos
					lastchr = chrom
				else:
					continue
		else:
			continue
			
		

outvcf.close()

