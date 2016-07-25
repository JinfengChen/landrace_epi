#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import gzip
import time
import multiprocessing as mp
from scipy.stats import binom_test, fisher_exact
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
def cal_filesize_gzip(filename):
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
def chunk_files_gzip(infile, numbers):
    #filesize  = cal_filesize_gzip(infile)
    #splitline = filesize//numbers
    unzip_infile = os.path.splitext(infile)[0]
    if not os.path.exists(unzip_infile):
        os.system('gunzip -c %s > %s' %(infile, unzip_infile))
    if not os.path.exists('%s_part00' %(unzip_infile)):
        filesize  = cal_filesize_gzip(infile)
        splitline = filesize//numbers 
        os.system('split -l %s %s %s_part -d' %(splitline, unzip_infile, unzip_infile)) 

#get gziped file size: line numbers
def cal_filesize(filename):
    f = open(filename, 'r')                  
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
    #filesize  = cal_filesize(infile)
    #splitline = filesize//numbers
    unzip_infile = os.path.splitext(infile)[0]
    #if not os.path.exists(unzip_infile):
    #    os.system('gunzip -c %s > %s' %(infile, unzip_infile))
    if not os.path.exists('%s_part00' %(infile)):
        filesize  = cal_filesize(infile)
        splitline = filesize//numbers 
        os.system('split -l %s %s %s_part -d' %(splitline, infile, infile)) 


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
                p    = binom_test(x, n, p=float(non_conversion_rate))
                #p    = 'NA'
                unit.append(str(p))
                print >> ofile, '\t'.join(unit)
    ofile.close()
    return 1

def binomial_test_helper(args):
    return binomial_test_file(*args)

##run function with parameters using multiprocess of #cpu
#def multiprocess_pool(function_helper, parameters, cpu):
#    pool = mp.Pool(int(cpu))
#    imap_it = pool.map(function_helper, tuple(parameters))
#    collect_list = []
#    for x in imap_it:
        #print 'status: %s' %(x)
#        collect_list.append(x)
#    return collect_list

##run function with parameters using multiprocess of #cpu
def multiprocess_pool_fisher(parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(fisher_test_helper, tuple(parameters))
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
    chunk_files_gzip(infile, 100)

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
#output:
#Chr1    850     1050	notmCG	mCG	#CG	notmCHG     mCHG    #CHG  notmCHH    mCHH    #CHH 
def sum_window_cytosine(infile, mini_depth):
    #p value for methylated C for each context
    #p = defaultdict(lambda : float())
    #p['CG'] = 0.001
    #p['CHG'] = 0.001
    #p['CHH'] = 0.001
    outfile = '%s.window_sum' %(infile)
    if os.path.exists(outfile):
        return outfile
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
                #if unit[3] == '.':
                #    sum_line = '%s\tNA\tNA\tNA\tNA\tNA\tNA' %('\t'.join(re.split(r'_', current_win)))
                #    print >> ofile, sum_line
                #    continue
                #have overlapping cytosine
                if not current_win == last_win:
                    #current window is not last win
                    #write last window into file if there is last win record
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
                        sum_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('\t'.join(re.split(r'_', last_win)), str(data[last_win]['CG'][0]), str(data[last_win]['CG'][1]), str(data[last_win]['CG'][3]), str(data[last_win]['CHG'][0]), str(data[last_win]['CHG'][1]), str(data[last_win]['CHG'][3]), str(data[last_win]['CHH'][0]), str(data[last_win]['CHH'][1]), str(data[last_win]['CHH'][3]))
                        print >> ofile, sum_line
                        del data[last_win]

                    #update record to current windows
                    last_win    = current_win
                    #mC vs. C
                    #no overlapping cytosine
                    if unit[3] == '.':
                        data[current_win]['CG'][2] = 0
                        data[current_win]['CHG'][2] = 0
                        data[current_win]['CHH'][2] = 0
                        data[current_win]['CG'][1] = 0
                        data[current_win]['CHG'][1] = 0
                        data[current_win]['CHH'][1] = 0
                        continue
                    #have overlapping cytosine
                    #base has mini_depth coverage
                    if int(unit[11]) >= mini_depth:
                        data[current_win][unit[7]][3] += 1
                    #total
                    data[current_win][unit[7]][2] += int(unit[11])
                    #mC
                    data[current_win][unit[7]][1] += int(unit[10])
                else:
                    #current window is last win, add value to existing key
                    #mC vs. C
                    #only consider base has mini_depth coverage
                    if int(unit[11]) >= mini_depth:
                        data[current_win][unit[7]][3] += 1
                    data[current_win][unit[7]][2] += int(unit[11])
                    data[current_win][unit[7]][1] += int(unit[10])
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
    sum_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('\t'.join(re.split(r'_', last_win)), str(data[last_win]['CG'][0]), str(data[last_win]['CG'][1]), str(data[last_win]['CG'][3]), str(data[last_win]['CHG'][0]), str(data[last_win]['CHG'][1]), str(data[last_win]['CHG'][3]), str(data[last_win]['CHH'][0]), str(data[last_win]['CHH'][1]), str(data[last_win]['CHH'][3]))
    print >> ofile, sum_line 
    return outfile

 

def chr_summary(bed, mini_depth, control_methC, treat_methC):
    bed_control_olp = '%s.control.overalp' %(os.path.splitext(bed)[0])
    bed_treat_olp = '%s.treat.overalp' %(os.path.splitext(bed)[0])
    cmd1 = 'bedtools intersect -a %s -b %s -wao -sorted > %s 2> %s.log' %(bed, control_methC, bed_control_olp, bed_control_olp)
    cmd2 = 'bedtools intersect -a %s -b %s -wao -sorted > %s 2> %s.log' %(bed, treat_methC, bed_treat_olp, bed_treat_olp)
    print cmd1
    print cmd2
    if not os.path.exists(bed_control_olp):
        os.system(cmd1)
    if not os.path.exists(bed_treat_olp):
        os.system(cmd2)
    control_sum = sum_window_cytosine(bed_control_olp, mini_depth) 
    treat_sum   = sum_window_cytosine(bed_treat_olp, mini_depth)
    if not os.path.exists('%s.control_treat.window_sum' %(os.path.splitext(bed)[0])):
        cmd3 = 'paste %s %s > %s.control_treat.window_sum' %(control_sum, treat_sum, os.path.splitext(bed)[0]) 
        os.system(cmd3)
    return '%s.control_treat.window_sum' %(os.path.splitext(bed)[0])

def chr_summary_helper(args):
    return chr_summary(*args)
#			3        4      5      6       7       8         9       10     11      12      13     14     15       16      17     18      19    20     21      22    23
#			 CG	 mCG	#CG    CHG     mCHG    #CHG    	 CHH	 mCHH   #CHH
#Chr12   23100   23300   1       21  	N      11      6       N         50      0       N     Chr12   23100   23300   6       16      N      5       12     N     52      0     N
def fisher_test_file(infile, mini_site):
    #print 'in fisher function'
    #print infile
    ofile = open('%s.fisher_test.txt' %(infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #print unit
                #CG
                if not 'NA' in [unit[3], unit[4], unit[15], unit[16]]:
                    if int(unit[5]) <= mini_site or int(unit[17]) <= mini_site:
                        print >> ofile, '%s\t%s\t%s\tCG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[3], unit[4], unit[15], unit[16], 'NA')  
                    else:
                        c1   = int(unit[3])
                        mc1  = int(unit[4])
                        c2   = int(unit[15])
                        mc2  = int(unit[16])
                        oddsratio, pvalue  = fisher_exact([[c1, mc1], [c2, mc2]])
                        print >> ofile, '%s\t%s\t%s\tCG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], str(c1), str(mc1), str(c2), str(mc2), str(pvalue))
                else:
                    print >> ofile, '%s\t%s\t%s\tCG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[3], unit[4], unit[15], unit[16], 'NA')
                #CHG
                if not 'NA' in [unit[6], unit[7], unit[18], unit[19]]:
                    if int(unit[8]) <= mini_site or int(unit[20]) <= mini_site:
                        print >> ofile, '%s\t%s\t%s\tCHG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[6], unit[7], unit[18], unit[19], 'NA')
                    else:
                        chg1   = int(unit[6])
                        mchg1  = int(unit[7])
                        chg2   = int(unit[18])
                        mchg2  = int(unit[19])
                        oddsratio, pvalue  = fisher_exact([[chg1, mchg1], [chg2, mchg2]])
                        print >> ofile, '%s\t%s\t%s\tCHG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], str(chg1), str(mchg1), str(chg2), str(mchg2), str(pvalue))
                else:
                    print >> ofile, '%s\t%s\t%s\tCHG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[6], unit[7], unit[18], unit[19], 'NA')
                #CHH
                if not 'NA' in [unit[9], unit[10], unit[21], unit[22]]:
                    if int(unit[11]) <= mini_site or int(unit[23]) <= mini_site:
                        print >> ofile, '%s\t%s\t%s\tCHH\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[9], unit[10], unit[21], unit[22], 'NA')
                    else:
                        chh1   = int(unit[9])
                        mchh1  = int(unit[10])
                        chh2   = int(unit[21])
                        mchh2  = int(unit[22])        
                        oddsratio, pvalue  = fisher_exact([[chh1, mchh1], [chh2, mchh2]])
                        print >> ofile, '%s\t%s\t%s\tCHH\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], str(chh1), str(mchh1), str(chh2), str(mchh2), str(pvalue))                
                else:
                    print >> ofile, '%s\t%s\t%s\tCHH\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[9], unit[10], unit[21], unit[22], 'NA')
    ofile.close()
    return 1

def fisher_test_helper(args):
    return fisher_test_file(*args)

#Chr1    1006    1006    C       CHH     CC      1.0     1       1       0.006   Chr1    1006    1006    C       CHH     CC      0.0     0       1       1.0     0
def DMC_fisher_test_file(infile, mini_depth):
    ofile = open('%s.fisher_test.txt' %(infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not 'NA' in [unit[7], unit[8], unit[17], unit[18]]:
                    c1   = int(unit[8]) - int(unit[7])
                    mc1  = int(unit[7])
                    c2   = int(unit[18]) - int(unit[17])
                    mc2  = int(unit[17])
                    #pvalue = 0.01
                    oddsratio, pvalue  = fisher_exact([[c1, mc1], [c2, mc2]])
                    unit[20] = str(pvalue)
                    print >> ofile, '\t'.join(unit)
                elif int(unit[8]) < int(mini_depth) or int(unit[18]) < int(mini_depth):
                    unit[20] = 'NA'
                    print >> ofile, '\t'.join(unit)
                else:
                    unit[20] = 'NA'
                    print >> ofile, '\t'.join(unit)
    ofile.close()
    return 1 

def DMC_fisher_test_helper(args):
    return DMC_fisher_test_file(*args)

##run function with parameters using multiprocess of #cpu
def multiprocess_pool(function_helper, parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(function_helper, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list

def P_adjust_BH(infile, p_column):
    outfile = '%s.P_adjusted_BH.txt' %(os.path.splitext(infile)[0])
    Rscript = '%s.P_adjusted_BH.R' %(os.path.splitext(infile)[0])
    Rcmd='''
#library(qvalue)
x <- read.table("%s")
p <- x[,%s]
q <- p.adjust(p, "BH")
#qStorey <- qvalue(p = p)
x_1 <- cbind(x, q)
write.table(x_1, file="%s", sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

''' %(infile, p_column, outfile)
    ofile = open(Rscript, 'w')
    print >> ofile, Rcmd
    ofile.close()
    os.system('cat %s | R --slave' %(Rscript))




#DMR info: 10 col
#Chr1    950     1150    CG      2       4       0       6       0.454545454545  1
#DMC info: 22 col
#Chr1    1128    1128    C       CG      CG      1       7       7       1e-21   Chr1    1128    1128    C       CG      CG      0.88    14      16      9.29867447152e-30       1       1	
def sum_window_diff_methylated_cytosine(infile, mini_diff):
    outfile = '%s.window_sum_diffmC' %(infile)
    if os.path.exists(outfile):
        return outfile
    ofile = open(outfile, 'w')
    data = defaultdict(lambda : defaultdict(lambda : int()))
    current_win = ''
    last_win = ''
    last_unit = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                current_win = '%s_%s_%s' %(unit[0], unit[1], unit[2])
                #no overlapping cytosine
                if unit[10] == '.':
                    sum_line = '%s\t%s' %('\t'.join(unit[:10]), 'NA')
                    continue
                #have overlapping cytosine
                if not current_win == last_win:
                    #print 'change: %s' %(line)
                    if data.has_key(last_win):
                        ##whether this is DMR: 1. number of differential methylated C in the windows
                        #if data[last_win] >= mini_diff:
                        sum_line = '%s\t%s\t%s' %('\t'.join(unit[:10]), str(data[last_win]['mC']), str(data[last_win]['C']))
                        print >> ofile, sum_line
                        del data[last_win]
                    last_win    = current_win
                    #print unit[10], unit[20], unit[30]
                    #only consider methylated C with Pvalue <= 0.05, which mean differential methylated in two sample
                    try:
                        data[current_win]['C'] += 1
                        if float(unit[30]) <= 0.05:
                            data[current_win]['mC'] += 1
                    except:
                        print 'error: %s' %(line)
                else:
                    #print unit[10], unit[20], unit[30]
                    try:
                        data[current_win]['C'] += 1
                        if float(unit[30]) <= 0.05:
                            data[current_win]['mC'] += 1
                    except:
                        print 'error: %s' %(line)
                last_unit = unit
    #last window
    sum_line = '%s\t%s\t%s' %('\t'.join(last_unit[:10]), str(data[last_win]['mC']), str(data[last_win]['C'])) 
    print >> ofile, sum_line
    return outfile



##run function with parameters using multiprocess of #cpu
#def multiprocess_pool_fisher(parameters, cpu):
#    pool = mp.Pool(int(cpu))
#    imap_it = pool.map(fisher_test_helper, tuple(parameters))
#    collect_list = []
#    for x in imap_it:
        #print 'status: %s' %(x)
#        collect_list.append(x)
#    return collect_list
 
def call_DMR(window, step, genome, mini_depth, control_methC, treat_methC, cpu, prefix, diff_mC_per_win, mini_site):
    print 'Step2. Call DMR by Fisher exact test'
    #bed files for chromosomes
    genome_window = '%s_w%s_s%s.bed' %(prefix, window, step)
    bed_chr_files = []
    if not os.path.exists(genome_window):
        os.system('bedtools makewindows -g %s -w %s -s %s | grep "^Chr" > %s' %(genome, window, step, genome_window))
        bed_chr_files = split_chr_files(genome_window)
    else:
        bed_chr_files = glob.glob('%s_w*_s*.Chr*.bed' %(prefix))
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
    #12 chromosome for rice genome, use 2 cpu because of large mem used by bedtools
    parameters = []
    bed_chr_files = sorted(bed_chr_files)
    control_methC_files = sorted(control_methC_files)
    treat_methC_files   = sorted(treat_methC_files)
    collect = []
    for i in range(0, len(bed_chr_files)):
        parameters.append([bed_chr_files[i], mini_depth, control_methC_files[i], treat_methC_files[i]])
    if 1:
        collect = multiprocess_pool(chr_summary_helper, parameters, 2) 
    #sum_window_cytosine('MSU7_w200_s50.Chr10.control.overalp')

    #cal DMR using Fisher exact test
    for table in sorted(collect):
        #fisher_test_file(table)
        chunk_files(table, 100)
        chunks = glob.glob('%s_part*' %(table)) 
        parameters = []
        for chunk in chunks:
            parameters.append([chunk, mini_site])
        tester = '%s_part00.fisher_test.txt' %(table)
        if not os.path.exists(tester):
            multiprocess_pool(fisher_test_helper, parameters, cpu)
        #merge_chr = 'cat %s_part*.fisher_test.txt > %s.fisher_test.txt' %(table, table)
        #os.system(merge_chr)
    merge_all = 'cat %s.Chr*.fisher_test.txt > %s.fisher_test.txt' %(os.path.splitext(genome_window)[0], os.path.splitext(genome_window)[0])
    merge_CG  = 'grep "CG" %s.fisher_test.txt > %s.fisher_test.CG.txt' %(os.path.splitext(genome_window)[0], os.path.splitext(genome_window)[0])
    merge_CHG = 'grep "CHG" %s.fisher_test.txt > %s.fisher_test.CHG.txt' %(os.path.splitext(genome_window)[0], os.path.splitext(genome_window)[0])
    merge_CHH = 'grep "CHH" %s.fisher_test.txt > %s.fisher_test.CHH.txt' %(os.path.splitext(genome_window)[0], os.path.splitext(genome_window)[0])
    if not os.path.exists('%s.fisher_test.txt' %(os.path.splitext(genome_window)[0])):
        os.system(merge_all)
    if not os.path.exists('%s.fisher_test.CG.txt' %(os.path.splitext(genome_window)[0])):
        time.sleep(60)
        os.system(merge_CG)
        os.system(merge_CHG)
        os.system(merge_CHH)
        P_adjust_BH('%s.fisher_test.CG.txt' %(os.path.splitext(genome_window)[0]), 9)
        P_adjust_BH('%s.fisher_test.CHG.txt' %(os.path.splitext(genome_window)[0]), 9)
        P_adjust_BH('%s.fisher_test.CHH.txt' %(os.path.splitext(genome_window)[0]), 9)

def Filer_DMR(window, step, genome, mini_depth, control_methC, treat_methC, cpu, prefix, diff_mC_per_win):
    genome_window = '%s_w%s_s%s.bed' %(prefix, window, step)
    #filter DMR using DMC
    if not os.path.exists('%s.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.txt' %(os.path.splitext(genome_window)[0])):
        for context in ["CG", "CHG", "CHH"]:
            DMR_beds = []
            DMC_beds = []
            if not os.path.exists('%s.fisher_test.%s.P_adjusted_BH.Chr1.bed' %(os.path.splitext(genome_window)[0], context)):
                DMR_beds = sorted(split_chr_files('%s.fisher_test.%s.P_adjusted_BH.txt' %(os.path.splitext(genome_window)[0], context)))
            else:
                DMR_beds = sorted(glob.glob('%s.fisher_test.%s.P_adjusted_BH.Chr*.bed' %(os.path.splitext(genome_window)[0], context)))
            if not os.path.exists('%s.cytosine_table.fisher_test.%s.P_adjusted_BH.Chr1.bed' %(prefix, context)):
                DMC_beds = sorted(split_chr_files('%s.cytosine_table.fisher_test.%s.P_adjusted_BH.txt' %(prefix, context)))
            else:
                DMC_beds = sorted(glob.glob('%s.cytosine_table.fisher_test.%s.P_adjusted_BH.Chr*.bed' %(prefix, context))) 
            #print 'DMR beds: %s' %(DMR_beds)
            #print 'DMC beds: %s' %(DMC_beds)
            for i in range(len(DMR_beds)):
                #print i
                overlap_chr_file = '%s.DMC_overlap.bed' %(os.path.splitext(DMR_beds[i])[0])
                overlap = 'bedtools intersect -a %s -b %s -wao > %s' %(DMR_beds[i], DMC_beds[i], overlap_chr_file)
                if not os.path.exists(overlap_chr_file):
                    os.system(overlap)
                sum_file = sum_window_diff_methylated_cytosine(overlap_chr_file, diff_mC_per_win)
            merge_chr = 'cat %s.fisher_test.%s.P_adjusted_BH.Chr*.DMC_overlap.bed.window_sum_diffmC > %s.fisher_test.%s.P_adjusted_BH.Filter_by_DMC.txt' %(os.path.splitext(genome_window)[0], context, os.path.splitext(genome_window)[0], context)
            pvalue    = '''awk '$9 <= 0.05 && $11 >= %s' %s.fisher_test.%s.P_adjusted_BH.Filter_by_DMC.txt > %s.fisher_test.%s.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt''' %(diff_mC_per_win, os.path.splitext(genome_window)[0], context, os.path.splitext(genome_window)[0], context)
            qvalue    = '''awk '$10 <= 0.05 && $11 >= %s' %s.fisher_test.%s.P_adjusted_BH.Filter_by_DMC.txt > %s.fisher_test.%s.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt''' %(diff_mC_per_win, os.path.splitext(genome_window)[0], context, os.path.splitext(genome_window)[0], context)
            os.system(merge_chr)
            os.system(pvalue)
            os.system(qvalue)

def call_DMC(genome, mini_depth, control_methC, treat_methC, cpu, prefix):

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
    control_methC_files = sorted(control_methC_files)
    treat_methC_files   = sorted(treat_methC_files)
 
    #DMC for each chromosome
    for i in range(len(control_methC_files)):
        chrs = re.split(r'\.', control_methC_files[i])[-2]
        #if not chrs == 'Chr1':
        #    continue
        table = '%s_%s.cytosine_table' %(prefix, chrs)
        overlap = '''bedtools intersect -a %s -b %s -wao -f 1 -sorted | awk '$12!~/\-1/' > %s''' %(control_methC_files[i], treat_methC_files[i], table)
        if not os.path.exists(table):
            os.system(overlap)
        chunk_files(table, 100)
        chunks = glob.glob('%s_part*' %(table))
        parameters = []
        for chunk in chunks:
            parameters.append([chunk, mini_depth])
        tester = '%s_part00.fisher_test.txt' %(table)
        if not os.path.exists(tester):
            multiprocess_pool(DMC_fisher_test_helper, parameters, cpu)  
    merge_all = 'cat %s_Chr*.*.fisher_test.txt > %s.cytosine_table.fisher_test.txt' % (prefix, prefix)
    merge_CG   = 'grep "CG" %s.cytosine_table.fisher_test.txt > %s.cytosine_table.fisher_test.CG.txt' %(prefix, prefix)
    merge_CHG  = 'grep "CHG" %s.cytosine_table.fisher_test.txt > %s.cytosine_table.fisher_test.CHG.txt' %(prefix, prefix)
    merge_CHH  = 'grep "CHH" %s.cytosine_table.fisher_test.txt > %s.cytosine_table.fisher_test.CHH.txt' %(prefix, prefix)
    if not os.path.exists('%s.cytosine_table.fisher_test.txt' %(prefix)): 
        os.system(merge_all)
    if not os.path.exists('%s.cytosine_table.fisher_test.CG.txt' %(prefix)):
        time.sleep(60)
        print merge_CG
        print merge_CHG
        print merge_CHH
        os.system(merge_CG)
        os.system(merge_CHG)
        os.system(merge_CHH)
        P_adjust_BH('%s.cytosine_table.fisher_test.CG.txt' %(prefix), 21)
        P_adjust_BH('%s.cytosine_table.fisher_test.CHG.txt' %(prefix), 21)
        P_adjust_BH('%s.cytosine_table.fisher_test.CHH.txt' %(prefix), 21)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--control')
    parser.add_argument('-t', '--treat')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-p', '--project')
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

    if not args.project:
        args.project = 'MSU7'

    ##cutoff setting
    #minimum read coverage required to consider cytosine to analysis: used in DMC and DMR, not in methylation call
    mini_depth      = 4
    #minimum cytosine site in window
    mini_site       = 4
    #minimum number of differetial methylated cytosine required to call DMR
    diff_mC_per_win = 4
    #fold changes of methylation level(average level among all cytosine in the window) between comparsion control vs. treat in a window
    fold_change     = 1.5

    #Step1. Call methylated Cytosine by binomial test
    control_methC = '%s.methylation_call.gz' %(os.path.splitext(args.control)[0])
    treat_methC = '%s.methylation_call.gz' %(os.path.splitext(args.treat)[0])
    if not os.path.exists(control_methC):
        methylated_C(args.control, args.cpu)
    if not os.path.exists(treat_methC):
        methylated_C(args.treat, args.cpu)

    #Step2. Call DMC by fisher's exact test
    if not os.path.exists('%s.cytosine_table.fisher_test.CG.P_adjusted_BH.txt' %(args.project)):
        call_DMC(args.genome, mini_depth, control_methC, treat_methC, args.cpu, args.project)

    #Step3. Call DMR by fisher's exact test
    window = 200
    step   = 50
    if not os.path.exists('%s_w%s_s%s.fisher_test.CG.P_adjusted_BH.txt' %(args.project, window, step)):
        call_DMR(window, step, args.genome, mini_depth, control_methC, treat_methC, args.cpu, args.project, diff_mC_per_win, mini_site)

    #Step4. Filter DMR
    if not os.path.exists('%s_w%s_s%s.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.txt' %(args.project, window, step)):
        Filer_DMR(window, step, args.genome, mini_depth, control_methC, treat_methC, args.cpu, args.project, diff_mC_per_win)

if __name__ == '__main__':
    main()

