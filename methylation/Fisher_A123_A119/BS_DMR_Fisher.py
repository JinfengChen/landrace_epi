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


#split gzip file into chromosomes base on chromosome name
#Chr1    0       200
def split_chr_files_gzip(infile):
    bed_files = []
    data = defaultdict(lambda : str())
    with gzip.open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'):
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    bed = '%s.%s.bed' %(os.path.splitext(infile)[0], unit[0])
                    bed_files.append(bed)
                    data[unit[0]] = open(bed, 'w')
                    print >> data[unit[0]], line
                else:
                    print >> data[unit[0]], line
    for chrs in data.keys():
        data[chrs].close()
    return bed_files

#split file into chromosomes base on chromosome name
#Chr1    0       200
def split_chr_files(infile):
    bed_files = []
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'):
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    bed = '%s.%s.bed' %(os.path.splitext(infile)[0], unit[0])
                    bed_files.append(bed)
                    data[unit[0]] = open(bed, 'w')
                    print >> data[unit[0]], line
                else:
                    print >> data[unit[0]], line
    for chrs in data.keys():
        data[chrs].close()
    return bed_files

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
def multiprocess_pool(function_helper, parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(function_helper, tuple(parameters))
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
#no q value is needed for this. Use p<=0.01 or p<1e-5 as suggest in? 
#Takuno and Gaut, Gene-body methylation is conserved between plant orthologs
#or
#Genome-Wide Analysis of DNA Methylation in Soybean
def methylated_C(infile, cpu):
    print 'Step1. Call methylated cytosine by binomial test'
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
        collect = multiprocess_pool(binomial_test_helper, parameters, cpu)
    if not os.path.exists('%s.methylation_call.gz' %(os.path.splitext(infile)[0])):
        cmd = '''cat %s_part*.methylation_call | sort -k1,1 -k3,3n | awk '{print $1"\t"$3"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > %s.methylation_call''' %(os.path.splitext(infile)[0], os.path.splitext(infile)[0])
        print cmd
        os.system(cmd)
        os.system('gzip %s.methylation_call' %(os.path.splitext(infile)[0]))

#summary mC vs. C for each window in each context
#Chr1    850     1050    Chr1    1001    1001    C       CHH     CT      1.0     1       1       0
#Chr1    850     1050    Chr1    1006    1006    C       CHH     CC      0.0     0       1       0
def sum_window_cytosine(infile):
    #p value for methylated C for each context
    p = defaultdict(lambda : float())
    p['CG'] = 0.001
    p['CHG'] = 0.001
    p['CHH'] = 0.001

    outfile = '%s.window_sum' %(infile)
    ofile = open(outfile, 'w')
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    current_win = ''
    last_win = ''
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                current_win = '%s_%s_%s' %(unit[0], unit[1], unit[2])
                #no overlapping cytosine
                if unit[3] == '.':
                    sum_line = '%s\tNA\tNA\tNA\tNA\tNA\tNA' %('\t'.join(re.split(r'_', current_win)))
                    print >> ofile, sum_line
                    continue
                #have overlapping cytosine
                if not current_win == last_win:
                    #print 'change: %s' %(line)
                    if data.has_key(last_win):
                        for c in ('CG', 'CHG', 'CHH'):
                            if data[last_win][c][2] > 0:
                                if data[last_win][c][1] > 0:
                                    data[last_win][c][0] = data[last_win][c][2] - data[last_win][c][1]
                                else:
                                    data[last_win][c][0] = data[last_win][c][2]
                                    data[last_win][c][1] = 0
                            else:
                                data[last_win][c][0] = 'NA'
                                data[last_win][c][1] = 'NA'
                        sum_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s' %('\t'.join(re.split(r'_', last_win)), str(data[last_win]['CG'][0]), str(data[last_win]['CG'][1]), str(data[last_win]['CHG'][0]), str(data[last_win]['CHG'][1]), str(data[last_win]['CHH'][0]), str(data[last_win]['CHH'][1]))
                        print >> ofile, sum_line
                        del data[last_win]
                    last_win    = current_win
                    #mC vs. C
                    data[current_win][unit[7]][2] += 1
                    if float(unit[12]) < p[unit[7]]:
                        data[current_win][unit[7]][1] += 1
                else:
                    #mC vs. C
                    #total cytosine
                    data[current_win][unit[7]][2] += 1
                    #mC
                    if float(unit[12]) < p[unit[7]]:
                        data[current_win][unit[7]][1] += 1
    #last window
    for c in ('CG', 'CHG', 'CHH'):
        if data[last_win][c][2] > 0:
            if data[last_win][c][1] > 0:
                data[last_win][c][0] = data[last_win][c][2] - data[last_win][c][1]
            else:  
                data[last_win][c][0] = data[last_win][c][2]
                data[last_win][c][1] = 0
        else:
            data[last_win][c][0] = 'NA'
            data[last_win][c][1] = 'NA'
    sum_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s' %('\t'.join(re.split(r'_', last_win)), str(data[last_win]['CG'][0]), str(data[last_win]['CG'][1]), str(data[last_win]['CHG'][0]), str(data[last_win]['CHG'][1]), str(data[last_win]['CHH'][0]), str(data[last_win]['CHH'][1]))
    print >> ofile, sum_line 
    return outfile

 

def chr_summary(bed, min_depth, control_methC, treat_methC):
    bed_control_olp = '%s.control.overalp' %(os.path.splitext(bed)[0])
    bed_treat_olp = '%s.treat.overalp' %(os.path.splitext(bed)[0])
    cmd1 = 'bedtools intersect -a %s -b %s -wao > %s 2> %s.log' %(bed, control_methC, bed_control_olp, bed_control_olp)
    cmd2 = 'bedtools intersect -a %s -b %s -wao > %s 2> %s.log' %(bed, treat_methC, bed_treat_olp, bed_control_olp)
    print cmd1
    print cmd2
    if not os.path.exists(bed_control_olp):
        os.system(cmd1)
    if not os.path.exists(bed_treat_olp):
        os.system(cmd2)
    control_sum = sum_window_cytosine(bed_control_olp) 
    treat_sum   = sum_window_cytosine(bed_treat_olp)
    cmd3 = 'paste %s %s > %s.control_treat.window_sum' %(control_sum, treat_sum, os.path.splitext(bed)[0]) 
    os.system(cmd3)

def chr_summary_helper(args):
    return chr_summary(*args)
 
def call_DMR(window, step, genome, mini_depth, control_methC, treat_methC):
    print 'Step2. Call DMR by Fisher exact test'
    #bed files for chromosomes
    genome_window = 'MSU7_w%s_s%s.bed' %(window, step)
    bed_chr_files = []
    if not os.path.exists(genome_window):
        os.system('bedtools makewindows -g %s -w %s -s %s | grep "^Chr" > %s' %(genome, window, step, genome_window))
        bed_chr_files = split_chr_files(genome_window)
    else:
        bed_chr_files = glob.glob('MSU7_w*_s*.Chr*.bed')
    #methC file for chromosomes
    tester = '%s.Chr1.bed' %(os.path.splitext(control_methC)[0])
    control_methC_files = []
    treat_methC_files = []
    if not os.path.exists(tester):
        control_methC_files = split_chr_files_gzip(control_methC)
        treat_methC_files = split_chr_files_gzip(treat_methC)
    else:
        control_methC_files = glob.glob('%s.Chr*.bed' %(os.path.splitext(control_methC)[0]))
        treat_methC_files  = glob.glob('%s.Chr*.bed' %(os.path.splitext(treat_methC)[0]))
    #bedtools intersect -a test.bed -b 1.1.bed -wao | less -S 
    
    #summary cytosine for window on each chromosome
    #12 chromosome for rice genome
    cpu = 2
    parameters = []
    bed_chr_files = sorted(bed_chr_files)
    control_methC_files = sorted(control_methC_files)
    treat_methC_files   = sorted(treat_methC_files)
    for i in range(0, len(bed_chr_files)):
        parameters.append([bed_chr_files[i], mini_depth, control_methC_files[i], treat_methC_files[i]])
    if 1:
        collect = multiprocess_pool(chr_summary_helper, parameters, cpu) 
    #sum_window_cytosine('MSU7_w200_s50.Chr10.control.overalp')

    #cal DMR using Fisher exact test

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--control')
    parser.add_argument('-t', '--treat')
    parser.add_argument('-g', '--genome')
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

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/seqlib/GFF/MSU7/MSU_r7.ALL.chrlen'

    #Step1. Call methylated Cytosine by binomial test
    control_methC = '%s.methylation_call.gz' %(os.path.splitext(args.control)[0])
    treat_methC = '%s.methylation_call.gz' %(os.path.splitext(args.treat)[0])
    if not os.path.exists(control_methC):
        methylated_C(args.control, args.cpu)
    if not os.path.exists(treat_methC):
        methylated_C(args.treat, args.cpu)

    #Step2. Call DMPR by fisher's exact test
    window = 200
    step   = 50
    mini_depth = 4
    call_DMR(window, step, args.genome, mini_depth, control_methC, treat_methC)

if __name__ == '__main__':
    main()

