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
    filesize  = cal_filesize_gzip(infile)
    splitline = filesize//numbers
    unzip_infile = os.path.splitext(infile)[0]
    if not os.path.exists(unzip_infile):
        os.system('gunzip -c %s > %s' %(infile, unzip_infile))
    if not os.path.exists('%s_part00' %(unzip_infile)): 
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
    filesize  = cal_filesize(infile)
    splitline = filesize//numbers
    unzip_infile = os.path.splitext(infile)[0]
    #if not os.path.exists(unzip_infile):
    #    os.system('gunzip -c %s > %s' %(infile, unzip_infile))
    if not os.path.exists('%s_part00' %(infile)): 
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
                p    = binom_test(x, n, p=non_conversion_rate)
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
def sum_window_cytosine(infile):
    #p value for methylated C for each context
    p = defaultdict(lambda : float())
    p['CG'] = 0.001
    p['CHG'] = 0.001
    p['CHH'] = 0.001
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
    cmd2 = 'bedtools intersect -a %s -b %s -wao > %s 2> %s.log' %(bed, treat_methC, bed_treat_olp, bed_treat_olp)
    print cmd1
    print cmd2
    if not os.path.exists(bed_control_olp):
        os.system(cmd1)
    if not os.path.exists(bed_treat_olp):
        os.system(cmd2)
    control_sum = sum_window_cytosine(bed_control_olp) 
    treat_sum   = sum_window_cytosine(bed_treat_olp)
    if not os.path.exists('%s.control_treat.window_sum' %(os.path.splitext(bed)[0])):
        cmd3 = 'paste %s %s > %s.control_treat.window_sum' %(control_sum, treat_sum, os.path.splitext(bed)[0]) 
        os.system(cmd3)
    return '%s.control_treat.window_sum' %(os.path.splitext(bed)[0])

def chr_summary_helper(args):
    return chr_summary(*args)

#			 CG	 mCG	 CHG	 mCHG	 CHH	 mCHH
#Chr12   23100   23300   1       21      11      6       50      0       Chr12   23100   23300   6       16      5       12      52      0
def fisher_test_file(infile):
    print 'in fisher function'
    ofile = open('%s.fisher_test.txt' %(infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #CG
                if not 'NA' in [unit[3], unit[4], unit[12], unit[13]]:
                    c1   = int(unit[3])
                    mc1  = int(unit[4])
                    c2   = int(unit[12])
                    mc2  = int(unit[13])
                    oddsratio, pvalue  = fisher_exact([[c1, mc1], [c2, mc2]])
                    print >> ofile, '%s\t%s\t%s\tCG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], str(c1), str(mc1), str(c2), str(mc2), str(pvalue))
                else:
                    print >> ofile, '%s\t%s\t%s\tCG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[3], unit[4], unit[12], unit[13], 'NA')
                #CHG
                if not 'NA' in [unit[5], unit[6], unit[14], unit[15]]:
                    chg1   = int(unit[5])
                    mchg1  = int(unit[6])
                    chg2   = int(unit[14])
                    mchg2  = int(unit[15])
                    oddsratio, pvalue  = fisher_exact([[chg1, mchg1], [chg2, mchg2]])
                    print >> ofile, '%s\t%s\t%s\tCHG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], str(chg1), str(mchg1), str(chg2), str(mchg2), str(pvalue))
                else:
                    print >> ofile, '%s\t%s\t%s\tCHG\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[5], unit[6], unit[14], unit[15], 'NA')
                #CHH
                if not 'NA' in [unit[7], unit[8], unit[16], unit[17]]:
                    chh1   = int(unit[7])
                    mchh1  = int(unit[8])
                    chh2   = int(unit[16])
                    mchh2  = int(unit[17])        
                    oddsratio, pvalue  = fisher_exact([[chh1, mchh1], [chh2, mchh2]])
                    print >> ofile, '%s\t%s\t%s\tCHH\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], str(chh1), str(mchh1), str(chh2), str(mchh2), str(pvalue))                
                else:
                    print >> ofile, '%s\t%s\t%s\tCHH\t%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[7], unit[8], unit[16], unit[17], 'NA')
    ofile.close()
    return 1

def fisher_test_helper(args):
    print 'in helper'
    return fisher_test_file(*args)

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
x <- read.table("%s")
p <- x[,%s]
q <- p.adjust(p, "BH")
x_1 <- cbind(x, q)
write.table(x_1, file="%s", sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

''' %(infile, p_column, outfile)
    ofile = open(Rscript, 'w')
    print >> ofile, Rcmd
    ofile.close()
    os.system('cat %s | R --slave' %(Rscript))

##run function with parameters using multiprocess of #cpu
#def multiprocess_pool_fisher(parameters, cpu):
#    pool = mp.Pool(int(cpu))
#    imap_it = pool.map(fisher_test_helper, tuple(parameters))
#    collect_list = []
#    for x in imap_it:
        #print 'status: %s' %(x)
#        collect_list.append(x)
#    return collect_list
 
def call_DMR(window, step, genome, mini_depth, control_methC, treat_methC, cpu):
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
        print table
        #fisher_test_file(table)
        chunk_files(table, 100)
        chunks = glob.glob('%s_part*' %(table)) 
        parameters = []
        for chunk in chunks:
            print chunk
            parameters.append([chunk])
        tester = '%s_part00.fisher_test.txt' %(table)
        print tester
        if not os.path.exists(tester):
            print 'running fisher exact test'
            print parameters, cpu
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
        os.system(merge_CG)
        os.system(merge_CHG)
        os.system(merge_CHH)
        P_adjust_BH('%s.fisher_test.CG.txt' %(os.path.splitext(genome_window)[0]), 9)
        P_adjust_BH('%s.fisher_test.CHG.txt' %(os.path.splitext(genome_window)[0]), 9)
        P_adjust_BH('%s.fisher_test.CHH.txt' %(os.path.splitext(genome_window)[0]), 9)

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
    call_DMR(window, step, args.genome, mini_depth, control_methC, treat_methC, args.cpu)

if __name__ == '__main__':
    main()

