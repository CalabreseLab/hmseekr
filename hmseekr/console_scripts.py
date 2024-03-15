import sys
import argparse



import kmers
import train
import findhits

from seekr.__version__ import __version__


KMERS_DOC = """
Description: 
This program counts k-mers for multiple specified values of k and saves them
to a binary file that countains a dictionary of dictionaries

Details:
this function takes in fasta file such as query seq or background seqs 
it counts the kmer frequency for each kmer size in kvec, kvec can also be a single kmer size
the output of the function is a dictionary of dictionaries
with the outer dictionary keys as kmer size and the inner dictionary keys as kmer sequences and values as kmer counts

Example:
generate count files for kmer size 2, 3, 4 using mXist repeat A sequence as input fasta file
    $ hmseekr_kmers -fd './fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -dir ./counts/ -name repeatA

minimal code with all settings to default
    $ hmseekr_kmers -fd './fastaFiles/mXist_rA.fa' -k 4 

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""

TRAIN_DOC = """
Description: 
This program calcuate the emission matrix based on kmer counts 
Prepare transition matrix, states and saves them to a binary file 
that countains a dictionary of dictionaries

Details:
this function takes in kmer count file for sequences of interest (e.g. query seq, functional regions of a ncRNA)
and kmer count file for background sequences (e.g. null seq, transcriptome, genome)
these kmer count files are generated using the kmers function
with the k specified (which should be calculted in the kmers function)
and query to query transition rate (qT) and null to null transition rate (nT) specified
it calculates the Hidden state transition matrix, Hidden state emission matrix, states and Starting probability of each hidden state 
and saves them to a binary file

Example:
train a model using previously generated kmer count files for repeatA and all lncRNA (hmseekr_kmers function) with kmer size 4
and transition rates of 0.9999 for both query to query and null to null
    $ hmseekr_train -qd './counts/repeatA.dict' -nd './counts/all_lncRNA.dict' -k 4 -a ATCG -qT 0.9999 -nT 0.9999 -qPre repeatA -nPre lncRNA -dir ./models/

minimal code with all settings to default
    $ hmseekr_train -qd './counts/repeatA.dict' -nd './counts/all_lncRNA.dict' -k 4

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


FINDHITS_DOC = """
Description: 
This program use precalculated model (emission matrix, prepared transition matrix, pi and states) from train function
to find out HMM state path through sequences of interest
therefore return sequences that has high similarity to the query sequence -- hits sequneces

Details:
this function takes in a fasta file which defines the region to search for potential hits (highly similar regions) to the query seq
also takes in the precalculated model (hmseekr_train function)
along the searchpool fasta sequences, similarity scores to query seq will be calculated based on the model
hits segments (highly similar regions) would be reported along with the sequence header from the input fasta file, start and end location of the hit segment
kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment


Example:
use the previously trained model (hmm.dict by hmseekr_train function) to search for highly similar regions to query seq (repeatA)
within the pool.fa files (area of interest region to find sequences similar to query, could be all lncRNAs or just chromosome 6) 
with kmer size 4 and save the hit sequences while showing progress bar
    $ hmseekr_findhits -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -name 'hits' -a 'ATCG' -fa -pb

minimal code with all settings to default
    $ hmseekr_findhits -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -fa

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


def _parse_args_or_exit(parser):
    if len(sys.argv) == 1:
        # this means no arguments given
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()



def console_hmseekr_kmers():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=KMERS_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-fd","--fadir", type=str,help='Path to input fasta file', required=True)
    parser.add_argument("-k","--kvec",type=str,help="Comma delimited string of possible k-mer values. For example, 3,4,5 or just 4", required=True)
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
    parser.add_argument("-dir","--outputdir'",type=str,help='Directory to save output count file',default='./')
    parser.add_argument("-name","--outputname",type=str,help='Desired output name for count file',default='out')
    
    args = _parse_args_or_exit(parser)

    kmers.kmers(
        args.fadir,
        args.kvec,
        args.alphabet,
        args.outputdir,
        args.outputname
        )
    

def console_hmseekr_train():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=TRAIN_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-qd","--querydir",type=str,help='Path to kmer count file for query sequence or sequence of interest (e.g. functional regions of a ncRNA)', required=True)
    parser.add_argument("-nd","--nulldir", type=str,help='Path to kmer count file that compose null model or bakground sequences (e.g. transcriptome, genome, etc.)', required=True)
    parser.add_argument("-k","--kvec",type=str,help='Comma delimited string of possible k-mer values, must be found in the k-mer count file', required=True)
    parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-qT","--queryT",type=float,help='Probability of query to query transition', default=0.9999)  
    parser.add_argument("-nT","--nullT",type=float,help='Probability of null to null transition', default=0.9999) 
    parser.add_argument("-qPre","--queryPrefix",type=str,help='prefix file name for query', default='query')
    parser.add_argument("-nPre","--nullPrefix",type=str,help='prefix file name for null', default='null')
    parser.add_argument("-dir","--outputdir",type=str,help='path of output directory to save output trained model file',default='./')
    
    args = _parse_args_or_exit(parser)

    train.train(
        args.querydir,
        args.nulldir,
        args.kvec,
        args.alphabet,
        args.queryT,
        args.nullT,
        args.queryPrefix,
        args.nullPrefix,
        args.outputdir
        )


def console_hmseekr_findhits():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=FINDHITS_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-pool","--searchpool",type=str,help='Path to fasta file from which the similarity scores to query seq will be calculated and hits seqs be found', required=True)
    parser.add_argument("-m","--modeldir",type=str,help='Path to precalculated model .dict file output from train.py', required=True)
    parser.add_argument("-k","--knum",type=int,help='Value of k to use as an integer. Must be the same as the k value used in training (train function)', required=True)
    parser.add_argument("-name","--outputname",type=str,help='File name for output, useful to include information about the experiment', default='hits')
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-fa","--fasta",action='store_true',help='FLAG: save sequence of hit, ignored if --wt is passed')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    findhits.findhits(
        args.searchpool,
        args.modeldir,
        args.knum,
        args.outputname,
        args.alphabet,
        args.fasta,
        args.progressbar
        )
    

