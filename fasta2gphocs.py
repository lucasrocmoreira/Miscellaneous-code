#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 19:47:44 2018

@author: lexlu

Python 3
"""

import sys # used for arguments
import os # used for working with file paths
from Bio import SeqIO # used to input fasta files
#from Bio import pairwise2
from Bio.Align.Applications import MafftCommandline
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from Bio import AlignIO

def align(fasta):
    in_file = os.path.relpath(fasta)
    mafft_cline = MafftCommandline(input=in_file)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    sequence1=str(align[0].seq)
    sequence2=str(align[1].seq)
    return [sequence1,sequence2]
	
def main(argv):
    # are arguments given correctly by user?
    if len(argv)<3:
        print("Arguments: Fasta1, Fasta2, Output")
        return
    
    # if yes, continue on
    output = argv[-1]
    
    name1 = argv[0].split(".")[0]
    name2 = argv[1].split(".")[0]
    longname = max(map(len,[name1,name2]))+4
    
    fasta1 = os.path.relpath(argv[0])
    fasta2 = os.path.relpath(argv[1])
    
    # use SeqRecord object to hold file data
    seq_1 = SeqIO.to_dict(SeqIO.parse(fasta1, "fasta"), lambda record: record.id.split('_')[0])
    seq_2 = SeqIO.to_dict(SeqIO.parse(fasta2, "fasta"), lambda record: record.id.split('_')[0])
    
    #get intersection
    intersect = seq_1.keys() & seq_2.keys()
    
    number_of_loci = len(intersect)
    
    #Write G-PhoCS file
    outfile = open(output, 'w')
    
    # Print the header, the number of loci in this file
    outfile.write(str(number_of_loci) + "\n\n")
    
    #Loop through loci
    for i in range(number_of_loci):
        locus = list(intersect)[i]
        
        #Aligned direct in biopython (not ideal)
#        seq1 = str(seq_1[locus].seq)
#        seq2 = str(seq_2[locus].seq)        
#        alignments = pairwise2.align.globalms(seq1, seq2, 1, -2, -15, -6.6666)
#        sequence1 = alignments[0][0]
#        sequence2 = alignments[0][1]

        seq1 = seq_1[locus]
        seq2 = seq_2[locus]
        #write locus fasta for alignment
        SeqIO.write([seq1,seq2], locus+".fasta", "fasta")
        
        #aligning with MAFFT
        alignments = align(locus+".fasta")
        sequence1 = alignments[0]
        sequence2 = alignments[1]

        #Delete alignment file
        os.remove(locus+".fasta")
        
        # Print out the header for this locus
        # <locus_name> <n_samples> <locus_length>
        outfile.write('{} {} {}\n'.format(list(intersect)[i], '2', len(sequence1)))
        
        #write sequences
        outfile.write(name1 + " "*(longname-len(name1))+ sequence1 + "\n")
        outfile.write(name2 + " "*(longname-len(name2))+ sequence2 + "\n")
        
        ## Separate loci with so it's prettier
        outfile.write("\n")

        ###print status
        print("Locus "+locus+" done!"+"\n"
              +str(i+1)+" of "+str(number_of_loci)+" loci completed")
        
    return

if __name__ == "__main__":
    main(sys.argv[1:])
