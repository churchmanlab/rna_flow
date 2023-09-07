# Written 1/28/21 separate vcf files into those that should be used for alternate allele and those that should be masked (N)
# Update 8/26/23 to add a quality filter (e.g. 20 (Phred score 20, means 99% confidence it is real)

import sys

infilename = sys.argv[1]
AFcutoff = float(sys.argv[2])
quality = float(sys.argv[3])
# J1-tot_snpcalls.vcf
invcf = open(infilename, 'r')
outvcf1 = open(infilename.replace('.vcf', '.AFatleast' + str(AFcutoff) + '.vcf'), 'w')
outvcf2 = open(infilename.replace('.vcf', '.AFbelow' + str(AFcutoff) + '.vcf'), 'w')

for line in invcf:
	if line[0] == '#':
		outvcf1.write(line)
		outvcf2.write(line)
	else:
		col = line.replace('\n', '').split('\t')
		fieldslist = col[7].split(';')
		# Get values from the DP4 field (first 2 are REF, second 2 are ALT)
		subs='DP4'
		DP4 = str([i for i in fieldslist if subs in i])
		DP4 = DP4.split("'")[1]
		DP4 = DP4.split('=')[1].split(',')
		DP4 = [float(i) for i in DP4] # Convert list elements to float
		# Calculate allele frequency
		tot = sum(DP4)
		AF = (DP4[2] + DP4[3])/tot
		qual = float(col[5])
	
		# Apply conditions and write
		if tot > 4 and AF >= AFcutoff and qual >= quality:
			outvcf1.write(line)
		elif qual >= quality:
			outvcf2.write(line)
		else:
			continue




outvcf1.close()
outvcf2.close()

