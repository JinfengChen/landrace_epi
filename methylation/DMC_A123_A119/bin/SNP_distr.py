#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python SNP_distr.py --bed A123_A119.DMC_Q0.05.core.bed --gff ../input/MSU_r7.gene_TE.bed

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#select distribution from collection of overlap
#1. if transposon just return transposon
#2. if mRNA which is from 5utr to 3utr, we continue to exon, intron and utrs
#3. others, intergenic but not transposon just return what it is
def position(pos):
    if len(pos) > 1:
        for i in pos:
            if i == 'Transposon':
                return i 
            elif i == 'mRNA':
                continue
            else:
                return i
    else:
        return pos[0]

#Chr1    RelocaTE        mPing   1041521 1041523 .       .       .       ID=Chr1_1041521;Strains=A119_2;GT=homozygous
#Chr1    MSU_osa1r7      mRNA    1046604 1053166 .       +       .       ID=LOC_Os01g02890.1;Name=phosphatidylserine synthase,putative,expressed;        5081
'''
Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous  Chr1    MSU_osa1r7      mRNA    10048264        10050309        .       -       .       ID=LOC_Os
01g17470.1
'''
def readmrna(infile):
    data = defaultdict(list)
    dist5 = defaultdict(lambda : list)
    dist3 = defaultdict(lambda : list)
    total = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0] + '_' + unit[3] + '_' + unit[4]
                pt = unit[11]
                strand = unit[15]
                dist = int(unit[18])
                total += 1
                if dist == 0:
                    comp = 1 if abs(int(unit[3])-int(unit[12])) <= abs(int(unit[4])-int(unit[13])) else 0 # mping in gene and close to start
                    if strand == '+':
                        if comp == 1:
                            dist5[mping] = abs(int(unit[3])-int(unit[12]))
                            #print 'Five\t%s' %(line)
                        else:
                            dist3[mping] = abs(int(unit[4])-int(unit[13]))
                            #print 'Three\t%s' %(line)
                    else:
                        if comp == 1:
                            dist3[mping] = abs(int(unit[3])-int(unit[12]))
                            #print 'Three\t%s' %(line)
                        else:
                            dist5[mping] = abs(int(unit[4])-int(unit[13]))
                            #print 'Five\t%s' %(line)
                if pt == 'mRNA' and dist > 0:
                    comp = 1 if int(unit[3]) < int(unit[12]) else 0
                    if strand == '+':
                        if comp == 1:
                            dist5[mping] = -dist
                        else:
                            dist3[mping] = -dist
                    else:
                        if comp == 1:
                            dist3[mping] = -dist
                        else:
                            dist5[mping] = -dist
    return dist5, dist3


'''
Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous  Chr1    MSU_osa1r7      intergenic
Chr1    5406    5406    23      24      2       9       0.0143876138892852      Chr1    MSU_osa1r7      mRNA    2903    10817
'''
def readtable(infile, prefix):
    data = defaultdict(list)
    mping_inf = defaultdict(lambda : str())
    #dist5, dist3 = readmrna(mrnafile)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0] + '_' + unit[1] + '_' + unit[2]
                pt = unit[10]
                data[mping].append(pt)

    sumx = defaultdict(int)
    tt = 0
    for m in data.keys():
        #print m, data[m]
        p = position(data[m])
        ###for intergenic we classify these into upstream and downstream
        #if p == 'intergenic':
        #    if dist5.has_key(m):
        #        if dist5[m] >= -5000 and dist5[m] < 0:
        #            p = '5k_upstream'
        #    if dist3.has_key(m):
        #        if dist3[m] >= -5000 and dist3[m] < 0:
        #            p = '5k_downstream'
        ###for mRNA which is acturaly in intron of UTR we classify these in intron
        if p == 'mRNA':
            p = 'intron'
        mping_inf[m] = p
        sumx[p] += 1
        tt += 1
    types = ['intergenic', 'Transposon', 'three_prime_UTR', 'CDS', 'intron', 'five_prime_UTR']
    #prefix = os.path.splitext(infile)[0]
    ofile = open (prefix + '.position.distr', 'w')
    count = 0
    for p in types:
        count +=1
        print >> ofile, '%s\t%s\t%s\t%s' %(str(count),p,sumx[p],str(float(sumx[p])/int(tt)))

    ofile.close()
    #
    #ofile = open (prefix + '.position.distr.list', 'w')
    #print >> ofile, 'Total mPing: %s' %(tt)
    #for m in sorted(mping_inf.keys()):
    #    print >> ofile, '%s\t%s' %(m, mping_inf[m])
    #ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed')
    parser.add_argument('--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bed) > 0
    except:
        usage()
        sys.exit(2)
   
    if not args.output:
        args.output = 'A123_A119'

    cmd = 'bedtools intersect -a %s -b %s -wao > %s.overlap_gff.bed' %(args.bed, args.gff, args.output)
    if not os.path.exists('%s.overlap_gff.bed' %(args.output)):
        os.system(cmd) 
    readtable('%s.overlap_gff.bed' %(args.output), args.output)

if __name__ == '__main__':
    main()

