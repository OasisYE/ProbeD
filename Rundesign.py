#!/usr/bin/env python
# -*-coding:utf-8 -*-
# ** author:Oasis
# *****************************
import pysam
import sys,os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt

def file_exists(path):
    if not os.path.exists(path):
        raise Exception('Error: File ' + path + ' does not exist!!!')

def repeatstat(s):
    '''
    统计repeat比例
    '''
    rcount = 0;
    if len(s) == 0:
        print ("序列为空，请输入正确的序列信息")
        return
    for i in range(97, 123):
        if s.count(chr(i)) == 0:
            continue;
        else:
            rcount += s.count(chr(i))
    return rcount / len(s)

def nstat(s):
    '''
    统计N比例
    '''
    if len(s) == 0:
        print ("序列为空，请输入正确的序列信息")
        return
    return s.upper().count('N') / len(s)

def faout(s,chr,start,end):
    '''
    输出探针序列信息
    '''
    print('>{}:{}-{} Repeat:{:.3f} GC:{:.3f} Nrate:{:.3f}'.format(chr,
    start, end, repeatstat(s), GC(s),nstat(s)))
    print(s)

def probeTm(seq1, sal, form):
    """Calculates the melting temperature of a given sequence under the
    specified salt and formamide conditions."""
    tmval = ('%0.2f' % mt.Tm_NN(seq1, Na=sal))
    fcorrected = ('%0.2f' % mt.chem_correction(float(tmval), fmd=form))
    return fcorrected


def main(argv):
    fout = open("./test.fa", 'w+')
    fafile = '/Users/yeweijian/Downloads/data/hg19.fa'
    bedfile = '/Users/yeweijian/Downloads/data/test.bed'

    parser = argparse.ArgumentParser(description='python Rundesign.py ')
    parser.add_argument('--FA', type=str, default=fafile, help='the reference fasta file')
    parser.add_argument('--BED', type=str, default=bedfile, help='the target region file')
    args = parser.parse_args()

    fafile = args.FA
    bedfile = args.BED

    file_exists(fafile)
    file_exists(bedfile)

    #读取fasta文件
    fh = pysam.Fastafile(fafile)

    #sal = 390 #The mM Na+ concentration to be used for Tm
    #form = 50  #The percent formamide to be used for Tm

    #读取区间文件，提取序列信息
    for line in open(bedfile):
        chrom, start, end = line.rstrip().split('\t')
        start = int(start)
        end = int(end)
        regionsize = int(end) - int(start)

        #区间太小，直接取区间序列
        if regionsize <= 120:
            seq = Seq(fh.fetch(reference=chrom, start=start, end=end), IUPAC.unambiguous_dna)
            faout(seq, chrom, start, end)
            #Tm = probeTm(seq, sal, form)

            #print("%0.2f" % mt.Tm_NN(seq))

            print('>{}:{}-{} Repeat:{:.3f} GC:{:.3f} Nrate:{:.3f} Tm:{}'.format(chrom, start, end, repeatstat(seq), GC(seq),
                                                                          nstat(seq),mt.Tm_NN(seq)),file=fout)
            print (seq,file=fout)

        else:
            for p1 in range(start,end):
                p2 = start + 119
                if p2>=end:
                    break
                else:
                    seq = Seq(fh.fetch(reference=chrom, start=p1, end=p2), IUPAC.unambiguous_dna)
                    faout(seq,chrom,p1,p2)

if __name__ == "__main__":
    main(sys.argv[1:])

