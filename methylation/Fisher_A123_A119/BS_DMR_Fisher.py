#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import gzip
import multiprocessing as mp
from scipy.stats import binom_test
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
Define DMR from bisufite sequencing data from two sample without replicate.

python BS_DMR_Fisher.py --control A119R1.BSseeker.CGmap.gz --treat A123R1.BSseeker.CGmap.gz

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#get gziped file size: line numbers
def cal_filesize(filename):
    f = gzip.open(filename, 'rb')                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

#chunk file by line into "numbers" of files
def chunk_files(infile, numbers):
    filesize  = cal_filesize(infile)
    splitline = filesize//numbers
    unzip_infile = os.path.splitext(infile)[0]
    if not os.path.exists(unzip_infile):
        os.system('gunzip -c %s > %s' %(infile, unzip_infile))
    if not os.path.exists('%s_part00' %(unzip_infile)): 
        os.system('split -l %s %s %s_part -d' %(splitline, unzip_infile, unzip_infile)) 


#read non_conversion_rate from existing file
def read_rate(infile):
    rate = 0.00
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                rate = float(line)
    return rate

#Chr1    G       1118    CHG     CT      1.0     2       2
def non_conversion_rate_estimator(infile, outfile):
    total_read = 0
    non_conversion_read = 0
    ofile = open (outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                total_read += int(unit[7])
                non_conversion_read += int(unit[6])
    non_conversion_rate = '%.3f' %(float(non_conversion_read)/float(total_read))
    print >> ofile, non_conversion_rate
    ofile.close()
    return non_conversion_rate

##Chr1    G       1118    CHG     CT      1.0     2       2
def binomial_test_file(infile, non_conversion_rate):
    ofile = open('%s.methylation_call' %(infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                n    = int(unit[7])
                x    = int(unit[6])
                p    = binom_test(x, n, p=non_conversion_rate)
                #p    = 'NA'
                unit.append(str(p))
                print >> ofile, '\t'.join(unit)
    ofile.close()
    return 1

def binomial_test_helper(args):
    return binomial_test_file(*args)

##run function with parameters using multiprocess of #cpu
def multiprocess_pool(parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(binomial_test_helper, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list

#infile                                  rate    mC      totalC
#Chr1    G       1118    CHG     CT      1.0     2       2
#Chr1    G       1121    CHG     CA      1.0     2       2
#outfile							Status	P-value
#Chr1    G       1118    CHG     CT      1.0     2       2	C	1	
#Chr1    G       1121    CHG     CA      1.0     2       5	mC	0.001
#no q value is needed for this. Use p<=0.01 or p<1e-5 as suggest in 
#Takuno and Gaut, Gene-body methylation is conserved between plant orthologs
#or
#Genome-Wide Analysis of DNA Methylation in Soybean
def methylated_C(infile, cpu):
    infile = os.path.abspath(infile)    

    #non conversion rate
    chloroplast_file = '%s.chloroplast.CGmap' %(infile)
    non_conversion_file = '%s.chloroplast.no_conversion_rate.txt' %(infile)
    p = 0.00
    if not os.path.exists(non_conversion_file):
        os.system('zcat %s | grep "chrC" > %s' %(infile, chloroplast_file))
        p = non_conversion_rate_estimator(chloroplast_file, non_conversion_file)
    else:
        p = read_rate(non_conversion_file)
    print '%s, non conversion rate: %s' %(infile, p)

    #split file into chunks
    chunk_files(infile, 100)

    #call methylation
    chunks = glob.glob('%s_part*' %(os.path.splitext(infile)[0]))
    parameters = []
    for chunk in chunks:
        parameters.append([chunk, p])
    if not os.path.exists('%s_part00.methylation_call' %(os.path.splitext(infile)[0])):
        collect = multiprocess_pool(parameters, cpu)
    if not os.path.exists('%s.methylation_call.gz' %(os.path.splitext(infile)[0])):
        os.system('cat %s_part*.methylation_call | sort -k1,1 -k3,3n | awk '{print $1"\t"$3"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > %s.methylation_call' %(os.path.splitext(infile)[0], os.path.splitext(infile)[0]))
        os.system('gzip %s.methylation_call' %(os.path.splitext(infile)[0]))
   
    #DMR
    #bedtools intersect -a test.bed -b 1.1.bed -wao | less -S 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--control')
    parser.add_argument('-t', '--treat')
    parser.add_argument('--cpu')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.control) > 0 and len(args.treat) > 0
    except:
        usage()
        sys.exit(2)

    if not args.cpu:
        args.cpu = 8

    #Step1. Call methylated Cytosine by binomial test
    if not os.path.exists('%s.meth_C.gz' %(args.control)):
        methylated_C(args.control, args.cpu)
    if not os.path.exists('%s.meth_C.gz' %(args.treat)):
        methylated_C(args.treat, args.cpu)

if __name__ == '__main__':
    main()

