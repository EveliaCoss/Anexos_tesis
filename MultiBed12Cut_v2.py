#!/bin/python

''' 
Program MultiBed12Cut.py
This program take BED files and do cutting transcripts in genomes
the input information for this program is the output of Liftover.sh 
in adition this program can procese multi BEDFiles and create a big data base
Powered by Python

$ python  MultiBed12Cut.py Genomes/ Div_coexpresion/
J. Antonio Corona
08/01/2019
H:\\PilotoCMF
'''

import os
import re
import sys
import numpy as np
import pandas as pd
import csv
import Bio
from Bio import SeqIO

#Read BED files and list of species

try:
    GENOMES = sys.argv[1]
    
except:
    sys.exit("Can't read directory of genomes")

try:
    BEDS = sys.argv[2]
except:
    sys.exit("Can't read directory of liftover bed files")
'''
BEDS='Coolair/'
GENOMES='Genomes/'
'''
#list of sp

def ReadGenome(genome):
    genomedb = {}
    for seq_record in SeqIO.parse(genome, "fasta"):
        CHR=seq_record.id
        CHR=str(CHR)
        genomedb[CHR] = seq_record.seq
    return(genomedb)

def CutSeq(row,genome):
    x=row['start'][0]
    y=row['end'][0]
    z=row['chrom'][0]
    z=str(z)
    seq_record = genome[z]
    a=row['bstarts'][0]
    a=a.split(",")
    w=row['bsizes'][0]
    w=w.split(",")
    s=row['stran'][0]
    c=row['bcount'][0]
    #print(z,x,y,a,w,c)
    seq=''
    for i in range(0,c):
        k=int(a[i])+x
        j=k+int(w[i])
        if s == '+':
            n=seq_record[k:j]
        if s == '-':
            n=seq_record[k:j]
            n = n.reverse_complement()
        if s == '.':
            n=seq_record[k:j]
        n=str(n)
        seq=seq+n    
    return(seq)

def SumExons(row):
    x=row['bsizes']
    x=x.split(",")
    y=0
    for i in x:
        i = int(i)
        y = i+y
    suma = y
    return (suma)

Head = ('chrom','start','end','name','score','stran','tstar','tend','rgb','bcount','bsizes','bstarts')
BED12=pd.DataFrame(columns =  Head)
ESPECIES=["Arabidopsis_thaliana","Aethionema_arabicum","Arabidopsis_halleri","Arabidopsis_lyrata","Arabis_alpina","Boechera_stricta","Brassica_oleracea","Brassica_rapa","Camelina_sativa","Capsella_rubella","Eutrema_parvulum","Eutrema_salsugineum"]
outhead=('sp','gen','seq')
OUT=pd.DataFrame(columns=outhead)
c=0
for i in ESPECIES:
    #print (i)
    #read bed12 file for specie
    bedpath=BEDS+i+'.bed'
    print(bedpath)
    BED=pd.read_table(bedpath,index_col=None)
    BED.columns =  Head
    #read genome
    genome=GENOMES+ i+'.fa'
    genome=ReadGenome(genome)
    #read genes
    genes=set(BED['name'])
    genes=list(genes)
    for l1 in genes:
        gen = BED[BED['name']==l1]
        #filter form with more nucletides aligned and cut sequence from geno
        gen['ntsizes']=gen.apply(SumExons, axis=1)
        genmax=gen[gen['ntsizes']==gen['ntsizes'].max()]
        genmax=genmax.reset_index(drop=True)
        seq=CutSeq(genmax,genome)
        #write result
        #print(i,l1,seq)
        OUT.loc[c] = [i,l1,seq]
        c=c+1

#Use OUT to create Fasta file for each gen
listg=OUT['gen'].tolist()
listg=set(listg)
for g in listg:
    print (g)
    dfg=(OUT[OUT['gen']==g])
    #print (dfg)
    name = (g+'.fa')
    f = open (name,'w')
    for index, row in dfg.iterrows():
        #print(row['sp'], row['seq'])
        sp = row['sp']
        seq =row['seq']
        fasta='>'+sp+"\n"+seq+"\n"
        f.write(fasta)
    f.close()
