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
    #filesize  = cal_filesize_gzip(infile)
    #splitline = filesize//numbers
    unzip_infile = os.path.splitext(infile)[0]
    if not os.path.exists(unzip_infile):
        os.system('gunzip -c %s > %s' %(infile, unzip_infile))
    if not os.path.exists('%s_part00' %(unzip_infile)):
        filesize  = cal_filesize_gzip(infile)
        splitline = filesize//numbers 
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
                name = re.split(r'\.', unit[0])[0]
                data[unit[0]] = [name, unit[1]]
    return data

#Chr1    C       1001    CHH     CT      1.0     1       1
def convert_bsseek2_file(infile):
    ofile1 = open('%s.methykit.CG.table' %(infile), 'w')
    ofile2 = open('%s.methykit.CHG.table' %(infile), 'w')
    ofile3 = open('%s.methykit.CHH.table' %(infile), 'w') 
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
                if unit[3] == 'CG':
                    print >> ofile1, '\t'.join(map(str,[chrbase, chrs, base, strand, coverage, freqC, freqT]))
                elif unit[3] == 'CHG':
                    print >> ofile2, '\t'.join(map(str,[chrbase, chrs, base, strand, coverage, freqC, freqT]))
                elif unit[3] == 'CHH':
                    print >> ofile3, '\t'.join(map(str,[chrbase, chrs, base, strand, coverage, freqC, freqT]))
    ofile1.close()
    ofile2.close()
    ofile3.close()


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

def run_methykit(meta, cpu, context, prefix):
    files  = []
    sample = []
    design = []
    for f in sorted(meta.keys()):
        files.append('%s.methykit.%s.table' %(os.path.splitext(f)[0], context))
        sample.append(meta[f][0])
        design.append(meta[f][1])

    Rcmd='''
library(methylKit)
pdf("%s.%s.pdf")
file.list=list("%s",
               "%s")

# read the files to a methylRawList object: myobj
myobj=read(file.list, sample.id=list("%s","%s"), assembly="MSU7",treatment=c(%s))

#Get descriptive stats on methylatio
# plot methylation statistics on samples in myobj which is a class of methylRawList
getMethylationStats(myobj[[1]],plot=F,both.strands=F)
getMethylationStats(myobj[[2]],plot=F,both.strands=F)
#getMethylationStats(myobj[[3]],plot=F,both.strands=F)
#getMethylationStats(myobj[[4]],plot=F,both.strands=F)
getMethylationStats(myobj[[1]],plot=T,both.strands=F)
getMethylationStats(myobj[[2]],plot=T,both.strands=F)
#getMethylationStats(myobj[[3]],plot=T,both.strands=F)
#getMethylationStats(myobj[[4]],plot=T,both.strands=F)

getCoverageStats(myobj[[1]], plot = T, both.strands = F)
getCoverageStats(myobj[[2]], plot = T, both.strands = F)
#getCoverageStats(myobj[[3]], plot = T, both.strands = F)
#getCoverageStats(myobj[[4]], plot = T, both.strands = F)

#filter and normalize
#filter base with fewer than 4X coverage and more than 99.9 percetile of coverage in each sample?
#filtered.myobj = filterByCoverage(myobj, lo.count = 4, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
#normalized.myobj = normalizeCoverage(filtered.myobj)
#myobj = normalized.myobj

#cluster of sample
meth=unite(myobj)
#a correlation matrix
#getCorrelation(meth, plot = T)
# cluster all samples using correlation distance and plot hiarachical clustering
#clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
# principal component anlaysis of all samples.
#PCASamples(meth)

#DMC
# calculate differential methylation p-values and q-values
myDiff=calculateDiffMeth(meth)
# get differentially methylated regions with 25 percent difference and qvalue<0.01
myDiff25p=get.methylDiff(myDiff,difference=25,qvalue=0.01)
head(myDiff25p)
write.table(myDiff, file="%s.%s.DMC.txt", sep="\\t", quote=FALSE)

#DMR
tiles = tileMethylCounts(myobj, win.size = 200, step.size = 50, cov.bases = 5)
meth_region = unite(tiles)
myDiff_region =calculateDiffMeth(meth_region)
myDiff_region25p=get.methylDiff(myDiff_region,difference=25, qvalue=0.01)
head(myDiff_region25p)
write.table(myDiff_region, file="%s.%s.DMR.txt", sep="\\t", quote=FALSE)

dev.off()
''' %(prefix, context, files[0], files[1], sample[0], sample[1], ','.join(design), prefix, context, prefix, context)
    Rscript = '%s.%s.methykit.R' %(prefix, context)
    ofile = open(Rscript, 'w')
    print >> ofile, Rcmd
    ofile.close()
    os.system('cat %s | R --slave' %(Rscript))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--meta')
    parser.add_argument('-c', '--cpu')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)


    if not args.cpu:
        args.cpu = 8
    if not args.project:
        args.project = 'A123_A119.Methykit'


    #read meta
    meta = read_meta(args.meta)

    #Convert table
    for filename in sorted(meta.keys()):
        if os.path.exists('%s.methykit.CG.table' %(os.path.splitext(filename)[0])):
            #already have methykit table, skip all step
            continue
        chunk_files_gzip(filename, 100)
        chunks = glob.glob('%s_part*' %(os.path.splitext(filename)[0]))
        parameters = []
        for chunk in chunks:
            parameters.append([chunk])
        if not os.path.exists('%s_part00.methykit.CG.table' %(os.path.splitext(filename)[0])):
            collect = multiprocess_pool(convert_bsseek2_helper, parameters, args.cpu)
        if not os.path.exists('%s.methykit.CG.table' %(os.path.splitext(filename)[0])):
            merge_CG = 'cat %s_part*.methykit.CG.table > %s.methykit.CG.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
            merge_CHG = 'cat %s_part*.methykit.CHG.table > %s.methykit.CHG.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
            merge_CHH = 'cat %s_part*.methykit.CHH.table > %s.methykit.CHH.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
            os.system(merge_CG)
            os.system(merge_CHG)
            os.system(merge_CHH)
        #if not os.path.exists('%s.methykit.CG.table' %(os.path.splitext(filename)[0])):
        #    merge_CG  = 'grep "CG" %s.methykit.table > %s.methykit.CG.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
        #    merge_CHG  = 'grep "CHG" %s.methykit.table > %s.methykit.CHG.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
        #    merge_CHH  = 'grep "CHH" %s.methykit.table > %s.methykit.CHH.table' %(os.path.splitext(filename)[0], os.path.splitext(filename)[0])
        #    os.system(merge_CG)
        #    os.system(merge_CHG)
        #    os.system(merge_CHH)

    #methykit
    run_methykit(meta, args.cpu, 'CG', args.project)   
 
if __name__ == '__main__':
    main()

