###################################################################################################
### Description: 
# This program counts k-mers for multiple specified values of k and saves them
# to a binary file that countains a dictionary of dictionaries

### Details:
# this function takes in fasta file such as query seq or background seqs 
# it counts the kmer frequency for each kmer size in kvec, kvec can also be a single kmer size
# the output of the function is a dictionary of dictionaries
# with the outer dictionary keys as kmer size and the inner dictionary keys as kmer sequences and values as kmer counts


### Input:
# fadir: path to the input fasta file
# kvec: Comma delimited string of possible k-mer values. For example, '3,4,5' or just '4'
# alphabet: String, Alphabet to generate k-mers, default=ATCG
# outputname: name of output count file, default is 'out'
# outputdir: path of output directory to save output count file, default is current directory


### Output:
# a dictionary of dictionaries with the outer dictionary keys as kmer size and the inner dictionary keys as kmer sequences and values as kmer counts
# the output is saved as a binary file with .dict extension

### Example:
# from hmseekr import kmers

# testdict = kmers(fadir='./fastaFiles/mXist_rA.fa',kvec='2,3,4',
#                  alphabet='ATCG',outputname='repeatA',
#                  outputdir='./counts/')

########################################################################################################


import kmersc
import pickle
from itertools import starmap
from itertools import product
from hmseekr import corefunctions


def kmers(fadir,kvec,alphabet='ATCG',outputname='out', outputdir='./'):

    # Read in specified values of k, and the alphabet
    kVals = [int(i) for i in kvec.split(',')]
    alphabet = alphabet.upper()

    #Read in fasta file
    seqs = corefunctions.getCookedFasta(fadir)[1::2]


    #Join sequences together using $ delimiter character
    fString = '$'.join(seqs)
    #lenFString = sum([len(i) for i in fS])

    # Need to figure out how to deal with very long fasta files (~ 2-3X the size of the transcriptome in mice)
    # if lenFString >= 2147483647:
    #     fString='$'.join(fS[::10]).upper()

    #Split jobs onto processors and call kmers.pyx cython file
    # with pool.Pool(args.n) as multiN:
    #     jobs = multiN.starmap(kmers.main,product(*[[fString],kVals,[a]]))
    #     dataDict = dict(jobs)
    #     print(dataDict.keys())

    # call kmers.pyx cython file to conduct kmer calculation and get kmer count dictionary as return
    dataDict = dict(starmap(kmersc.main,product(*[[fString],kVals,[alphabet]])))


    #Save data
    kDir = outputdir
    if not kDir.endswith('/'):
        kDir+='/'
    pickle.dump(dataDict,open(f'{kDir}{outputname}.dict','wb'))

    return dataDict





