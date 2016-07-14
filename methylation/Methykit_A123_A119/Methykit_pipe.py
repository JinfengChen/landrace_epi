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

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

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


#A119R1.BSseeker.CGmap.gz	1
#A119R2.BSseeker.CGmap.gz	1
#A123R1.BSseeker.CGmap.gz	0
#A123R2.BSseeker.CGmap.gz	0
def read_meta(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[1]
    return data

#Chr1    C       1001    CHH     CT      1.0     1       1
def convert_bsseek2_file(infile):
    ofile = open('%s.methykit.table' %(infile), 'w') 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                chrbase = '%s.%s' %(unit[0], unit[2])
                chrs    = unit[0]
                base    = unit[2]
                strand  = 'F' if unit[1] == 'C' else 'R'
                coverage= unit[7]
                freqC   = 100 * float(unit[5])
                freqT   = 100 - freqC
                print >> ofile, '\t'.join(map(str,[chrbase, chrs, base, strand, coverage, freqC, freqT]))
    ofile.close()


def convert_bsseek2_helper(args):
    convert_bsseek2_file(*args)


##run function with parameters using multiprocess of #cpu
def multiprocess_pool(function_helper, parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(function_helper, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--meta')
    parser.add_argument('-c', '--cpu')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)


    if not args.cpu:
        args.cpu = 8


    #read meta
    meta = read_meta(args.meta)

    #Convert table
    for filename in sorted(meta.keys()):
        if os.path.exists('%s.methykit.table' %(os.path.splitext(filename)[0])):
            #already have methykit table, skip all step
            continue
        chunk_files_gzip(filename, 100)
        chunks = glob.glob('%s_part*' %(os.path.splitext(filename)[0]))
        parameters = []
        for chunk in chunks:
            parameters.append([chunk])
        if not os.path.exists('%s_part00.methykit.table' %(os.path.splitext(filename)[0])):
            collect = multiprocess_pool(convert_bsseek2_helper, parameters, args.cpu)
        if not os.path.exists('%s.methykit.table' %(os.path.splitext(filename)[0])):
            merge_all = 'cat %s_part*.methykit.table > %s.methykit.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
            os.system(merge_all)

    #methykit
    
if __name__ == '__main__':
    main()

