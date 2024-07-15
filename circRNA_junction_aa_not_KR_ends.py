#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2016-6-13
E-mail: fanxiaojuan@picb.ac.cn
Description: Convert DNA sequence around back-spliced junction to protein sequence (only annotated circRNA)
"""

import argparse

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--input file', dest='fnIn', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_database.fa', type=str,help='input file')
    parser.add_argument('-l', '--length', dest='l', default=90, type=int,help='input single end length')
    parser.add_argument('-junc_l', '--junc-length', dest='junc_l', default=1, type=int,help='input single end length')
    parser.add_argument('-o', '--output file', dest='fnOut', default='F:/experiment/Yangyun/IRES/CIRI2/open-pFind/circBase_ZFQ_junc_protein.fa', type=str,help='output file')
    op=parser.parse_args()
    return op

def DNA_to_pro(DNA,amino_acid_dic,op):
    """
    Convert DNA sequence to protein sequence
    """
    pro_list = []
    for i in xrange(0,len(DNA),2):
        seq_header = DNA[i].strip()
        DNA_seq = DNA[i+1].strip()
        if len(DNA_seq) >= 2*op.l:
            junc_seq = DNA_seq[-op.l:].upper() + DNA_seq[0:op.l].upper()
            #print junc_seq
            pro_list.append(seq_header + '@frame_1')
            pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(0,len(junc_seq),3)]))
            pro_list.append(seq_header + '@frame_2')
            pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(1,len(junc_seq)-2,3)]))
            pro_list.append(seq_header + '@frame_3')
            pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(2,len(junc_seq)-1,3)]))
        else:
            junc_seq = DNA_seq[len(DNA_seq)/2:].upper() + DNA_seq[0:len(DNA_seq)/2].upper()
            #print junc_seq
            if len(junc_seq) % 3 == 0:
                pro_list.append(seq_header + '@frame_1')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(0,len(junc_seq),3)]))
                pro_list.append(seq_header + '@frame_2')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(1,len(junc_seq)-2,3)]))
                pro_list.append(seq_header + '@frame_3')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(2,len(junc_seq)-1,3)]))
            if len(junc_seq) % 3 == 1:
                pro_list.append(seq_header + '@frame_1')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(0,len(junc_seq)-1,3)]))
                pro_list.append(seq_header + '@frame_2')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(1,len(junc_seq),3)]))
                pro_list.append(seq_header + '@frame_3')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(2,len(junc_seq)-2,3)]))
            if len(junc_seq) % 3 == 2:
                pro_list.append(seq_header + '@frame_1')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(0,len(junc_seq)-2,3)]))
                pro_list.append(seq_header + '@frame_2')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(1,len(junc_seq)-1,3)]))
                pro_list.append(seq_header + '@frame_3')
                pro_list.append(''.join([amino_acid_dic[junc_seq[n:n+3]] for n in xrange(2,len(junc_seq),3)]))
            #print pro_dic
    #print pro_list
    return pro_list
    
def filter_MS_dataset(junc_list,op):
    """
    Keep sequence that longer than 7aa, without '-'
    """
    filter_list = []
    for i in xrange(0,len(junc_list),2):
        if '>' in junc_list[i]:
            if '-' not in junc_list[i+1]:
                if len(junc_list[i+1]) >= 7:
                    filter_list.append(junc_list[i])
                    filter_list.append(junc_list[i+1])
                    #print junc_list[i+1],junc_list[i+1][0:j+1]
                    continue
            elif junc_list[i+1].count('-') == 1:
                if junc_list[i+1].index('-') < len(junc_list[i+1])/2 - op.junc_l + 1:
                    #print junc_list[i+1][junc_list[i+1].index('-')+1:j+1]
                    if len(junc_list[i+1]) >= 7:
                        filter_list.append(junc_list[i])
                        filter_list.append(junc_list[i+1][junc_list[i+1].index('-')+1:])
                        continue
                elif junc_list[i+1].index('-') > len(junc_list[i+1])/2 + op.junc_l:
                    #print junc_list[i+1]
                    #print junc_list[i+1][len(junc_list[i+1])/2-15:j+1]
                    if len(junc_list[i+1]) >= 7:
                        filter_list.append(junc_list[i])
                        filter_list.append(junc_list[i+1][0:junc_list[i+1].index('-')])
                        #print junc_list[i+1],junc_list[i+1][0:j+1]
                        continue
            else:
                temp = []
                temp.append(-1) #temp[n] is the index of '-', if append 0 to temp, if temp[0] is not '-', that's wrong
                for j in xrange(0,len(junc_list[i+1])):
                    if junc_list[i+1][j] == '-':
                        temp.append(j)
                temp.append(len(junc_list[i+1]))
                #print temp
                for n in xrange(0,len(temp)-1):
                    if len(junc_list[i+1])/2 - op.junc_l + 1 > temp[n] and len(junc_list[i+1])/2 + op.junc_l < temp[n+1]:
                        #print temp[n+1]-1,len(junc_list[i+1])/2 + op.junc_l -1
                        if temp[n+1] - temp[n] - 1 >= 7:
                            filter_list.append(junc_list[i])
                            filter_list.append(junc_list[i+1][temp[n]+1:temp[n+1]])
                            break
    return filter_list

if __name__ == '__main__':
    op = createHelp()
    
    amino_acid_dic = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"-", "TAG":"-",
    "TGT":"C", "TGC":"C", "TGA":"-", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    DNA_seq = open(op.fnIn).readlines()
    junc_list = DNA_to_pro(DNA_seq,amino_acid_dic,op)
    
    junc_filter = filter_MS_dataset(junc_list,op)
    
    output = open(op.fnOut,'w')
    output.write('\n'.join(junc_filter) + '\n')
    
    print 'Done'