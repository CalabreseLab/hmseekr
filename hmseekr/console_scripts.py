import sys
import argparse



from hmseekr.kmers import kmers
from hmseekr.train import train
from hmseekr import findhits
from hmseekr import findhits_condE
from hmseekr import findhits_nol
from hmseekr.gridsearch import gridsearch
from hmseekr.fastarev import fastarev
from hmseekr.genbed import genbed
from hmseekr.genbedrev import genbedrev
from hmseekr.hitseekr import hitseekr
from hmseekr.seqstosummary import seqstosummary

from hmseekr.__version__ import __version__


KMERS_DOC = """
Description: 
This program counts k-mers for multiple specified values of k and saves them
to a binary file that countains a dictionary of dictionaries

Details:
this function takes in fasta file such as query seq or background seqs 
if input fasta contains more than one sequence, all the sequences will be merged to one sequence before processing
it counts the kmer frequency for each kmer size in kvec, kvec can also be a single kmer size
the output of the function is a dictionary of dictionaries
with the outer dictionary keys as kmer size and the inner dictionary keys as kmer sequences and values as kmer counts

Example:
generate count files for kmer size 2, 3, 4 using mXist repeat A sequence as input fasta file
    $ hmseekr_kmers -fd './fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -name repeatA -dir './counts/' 

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
this function takes in kmer count file for sequences of interest (e.g. query seq, functional regions of a lncRNA)
and kmer count file for background sequences (e.g. null seq, transcriptome, genome)
these kmer count files are generated using the kmers function
with the k specified (which should be calculted in the kmers function)
and query to query transition rate (qT) and null to null transition rate (nT) specified
it calculates the Hidden state transition matrix, Hidden state emission matrix, states and Starting probability of each hidden state 
and saves them to a binary file
qT and nT should be between 0 and 1 but not equal to 0 or 1

Example:
train a model using previously generated kmer count files for repeatA and all lncRNA (hmseekr_kmers function) with kmer size 4
and transition rates of 0.9999 for both query to query and null to null
    $ hmseekr_train -qd './counts/repeatA.dict' -nd './counts/all_lncRNA.dict' -k 4 -a ATCG -qT 0.9999 -nT 0.9999 -qPre repeatA -nPre lncRNA -dir './models/'

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
kmerLLR is defined as the sum of the log likelihood of each k-mer in the hit sequence being in the Q state minus the log likelihood of them being in the N state

Example:
use the previously trained model (hmm.dict by hmseekr_train function) to search for highly similar regions to query seq (repeatA)
within the pool.fa files (area of interest region to find sequences similar to query, could be all lncRNAs or just chromosome 6) 
with kmer size 4 and save the hit sequences while showing progress bar
    $ hmseekr_findhits -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -name 'hits' -dir './models/'  -a 'ATCG' -fa -pb

minimal code with all settings to default
    $ hmseekr_findhits -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -fa

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


FINDHITS_CONDE_DOC = """
Description: 
this function is different from the basic findhits function in that it calculates the emission probability of the next word given the current word

Details:
this is a variant of findhits/hmseekr_findhits function
in findhits, as the shift is by 1 nt, the next word has k-1 overlap with the current word
for example, if the current word is 'TAGC', the next possible words are 'AGCA', 'AGCT', 'AGCC', 'AGCG'
then the emission probability of the next word ('AGCA') given the current word as 'TAGC' is calculated as 
np.log2(emission probability of 'AGCA' in the original E) - np.log2(sum of emission probability of all four possible worlds in the original E)

Example:
same example as hmseekr_findhits but uses the conditioned emission probability of the next word given the current word
    $ hmseekr_findhits_condE -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -name 'hits' -dir './models/'  -a 'ATCG' -fa -pb

minimal code with all settings to default
    $ hmseekr_findhits_condE -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -fa

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


FINDHITS_NOL_DOC = """
Description: 
this function is different from the basic findhits in that it uses non-overlapping kmers and shifts the sequence by 1nt each time until the kmer size is reached
this reduces the kmer dependency of overlapping kmers

Details:
this is a variant of findhits/hmseekr_findhits function
if sequence is 'ATGCTTTTGCGC' the kmers would be 'ATGC','TTTT','GCGC'
then it chops off 1nt from the start of each sequence in the searchpool and re-scan with non-overlapping kmers
so the sequence would be 'TGCTTTTGCGC' and the kmers would be 'TGCT','TTTG'
the last bit of sequences that is less than the kmer size would be discarded
this is done until the kmer size is reached
in this way, it generates k result txt files each starts at different position of the sequence with non-overlapping kmers
this function reduces the kmer dependency of overlapping kmers in findhits, but also considered the different combinations of kmers at different start positions
the outcome txt files can be further process to find the best hit regions: for example, only use the regions that are hits in all k result txt files

Example:
same example as hmseekr_findhits but uses the non-overlapping kmers 
    $ hmseekr_findhits_nol -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -name 'hits' -dir './models/'  -a 'ATCG' -fa -pb

minimal code with all settings to default
    $ hmseekr_findhits_nol -pool './fastaFiles/pool.fa' -m './markovModels/hmm.dict' -k 4 -fa

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


GRIDSEARCH_DOC = """
Description: 
This function performs a grid search to find the best trasnition probabilities for query to query and null to null states
which is used in the train function

Details:
this function takes in query sequence, null squence and background sequence (see Inputs for details) fasta files
ranges and steps or values for query to query transition rate (qT) and null to null transition rate (nT)
a specific kmer number and performs the train function and findhits function for each combination of qT and nT
within the hits sequences (findhits function results), only keep the sequence with length greater than lengthfilter 
then calculates pearson correlation r score (seekr.pearson) between the filtered hit sequences (findhits results) and the query sequence
it returns a dataframe (.csv file) containing the qT, nT, kmer number, the total number of hits sequences and the mean, median, standard deviation of the hits sequences' pearson correlation r score to the query sequence
and the mean, median, standard deviation of the length of the hits sequences
if query fasta contains more than one sequence all the sequences in query fasta file will be merged to one sequence 
for calculating kmer count files for hmseekr (hmseekr.kmers) and for calculating seekr.pearson 
this function requires the seekr package to be installed
as there are iterations of train and findhits functions, which could take long time, it is recommended to run this function on a high performance computing cluster
variants of findhits functions can be specified to run

Example:
perform a grid search to find the best transition probabilities for qT and nT each within the range of 0.9 to 0.99 with step of 0.01 with lengthfilter set to 25
which only keep the hit sequences with length greater than 25 for stats calculation. the regular findhits function is used
    $ hmseekr_gridsearch -qf './fastaFiles/repeatA.fa' -nf './fastaFiles/all_lncRNA.fa' -pool './fastaFiles/pool.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -ql 0.9,0.99,0.01 -nl 0.9,0.99,0.01 -step -fc 'findhits' -lf 25 -name 'gridsearch_results' -dir './gridsearch/' -a 'ATCG' -pb

perform a grid search to find the best transition probabilities for qT and nT each exactly as 0.9,0.99,0.999 with lengthfilter set to 25
which only keep the hit sequences with length greater than 25 for stats calculation. the regular findhits function is used
    $ hmseekr_gridsearch -qf './fastaFiles/repeatA.fa' -nf './fastaFiles/all_lncRNA.fa' -pool './fastaFiles/pool.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -ql 0.9,0.99,0.999 -nl 0.9,0.99,0.999 -fc 'findhits' -lf 25 -name 'gridsearch_results' -dir './gridsearch/' -a 'ATCG' -pb


For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


FASTAREV_DOC = """
Description: 
This function reverses a fasta file and save it to a new file
'ATGC' will be reversed to 'CGTA'

Details:
this function takes in fasta file and reverse the sequence.
if input fasta contains more than one sequence, each sequence will be reversed
the headers will be kept the same, only the sequences will be reversed
keeping the headers the same will allow the user to match the reversed sequences to the original sequences easily
the output is a fasta file with the reversed sequences

Example:
reverse the sequence of mXist repeat A fasta file and save it to a new file while keeping the headers the same
    $ hmseekr_fastarev -i '../fastaFiles/mXist_rA.fa' -o '../fastaFiles/mXist_rA_rev.fa'

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""

GENBED_DOC = """
Description: 
This function firstly filter the hits output from hmseekr_findhits
and then generate a bed file for the filtered hits
this function only applies to the output from hmseekr_findhits with the regular fasta file as input
for reversed fasta file, please use hmseekr_genbedrev

Details:
this function takes in the output from hmseekr_findhits and filter the hits based on the hit length and normalized kmerLLR score
kmerLLR is the log likelihood ratio of of the probability 
that the set of k-mers y within a hit derived from the QUERY versus the NULL state
it is the sum of the log2(Q/N) ratio for each kmer within a hit
the normalized kmerLLR (normLLR) is the kmerLLR divided by the hit length
so the normLLR is the log2 of the Q/N ratio, i.e. if set normLLR > 0.5, the Q/N ratio is ~ 1.41
then all the hits after the filtering will be converted to a bedfile


Example:
generate bedfile for filtered hits from the hmseekr_findhits output file mm10expmap_queryA_4_viterbi.txt
    $ hmseekr_genbed -hd '../mm10expmap_queryA_4_viterbi.txt' -o '../mm10expmap_queryA_4_viterbi' -len 25 -llr 0.5 -pb

minimal code with all settings to default
    $ hmseekr_genbed -hd '../mm10expmap_queryA_4_viterbi.txt' -o '../mm10expmap_queryA_4_viterbi' 

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


GENBEDREV_DOC = """
Description: 
This function firstly filter the hits output from hmseekr_findhits with the reversed fasta file as input
and then generate a bed file for the filtered hits
this function only applies to the output from hmseekr_findhits with the reversed fasta file as input
for regular or forward fasta file, please use hmseekr_genbed

Details:
this function takes in the output from hmseekr_findhits with the reversed fasta file as input
and filter the hits based on the hit length and normalized kmerLLR score
kmerLLR is the log likelihood ratio of of the probability 
that the set of k-mers y within a hit derived from the QUERY versus the NULL state
it is the sum of the log2(Q/N) ratio for each kmer within a hit
the normalized kmerLLR (normLLR) is the kmerLLR divided by the hit length
so the normLLR is the log2 of the Q/N ratio, i.e. if set normLLR > 0.5, the Q/N ratio is ~ 1.41
then all the hits after the filtering will be converted to a bedfile


Example:
generate bedfile for filtered hits from the hmseekr_findhits output file FLIPmm10expmap_queryA_4_viterbi.txt 
which is generated using the reversed fasta file as input
    $ hmseekr_genbedrev -hd '../FLIPmm10expmap_queryA_4_viterbi.txt' -o '../FLIPmm10expmap_queryA_4_viterbi' -len 25 -llr 0.5 -pb

minimal code with all settings to default
    $ hmseekr_genbedrev -hd '../FLIPmm10expmap_queryA_4_viterbi.txt' -o '../FLIPmm10expmap_queryA_4_viterbi' 

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


HITSEEKR_DOC = """
Description: 
This function takes in the results of findhits function and calculates the seekr pearson correlation p value between the hit sequences and the query sequence

Details:
this function takes in the output of findhits function and the query sequence fasta file with a background fasta file
it fit the background sequences to the common10 distributions and takes the best ranked distribution
it calculates the seekr pearson correlation p values between the hit sequences and the query sequence based on the best ranked distribution
it adds the seekr p value to the hits dataframe
on top of the existing kmer log likelihood score (kmerLLR)
the seekr p value could provide additional information about the similarity between the hit sequences and the query sequence
outputname is automatically generated as the input findhits filename with '_seekr' appended to it

Example:
add seekr p value to the hits dataframe generated by findhits function
    $ hmseekr_hitseekr -hd './mm10_queryA_4_viterbi.txt' -qf './fastaFiles/repeatA.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -lf 25 -dir './' -pb


For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues

"""


SEQSTOSUMMARY_DOC = """
Description: 
This function takes in a fasta file of multiple query sequences, a transtion probabilty dataframe, a search pool fasta file, a null fasta file and a background fasta file.
and generates a summary dataframe where each row contains a sequence in the search pool fasta file
and five columns: seqName, feature, counts, sum_normLLR, sum_len

Details:
this function is designed to get the overall likeliness of each search pool sequence to the query sequences
the function takes in a fasta file of multiple query sequences, a transtion probabilty dataframe, a search pool fasta file, a null fasta file (for hmseekr) and a background fasta file (for seekr)
here the transition probability dataframe must have the same rows as the query fasta file
the transition prbability for each query sequence can be different and can be optimized by the gridsearch function
users can also choose to set the transition probability to be the same for all query sequences
the function will run the kmers, train, findhits and hitseekr functions for each query sequence
then the results can be filtered by the length of the hit regions, the normalized kmer log likelihood ratio (kmerLLR) and the seekr pearson correlation p value
finally for each search pool sequence, the function will calculate the counts of filtered hit regions with a specific query sequence
and also the sum of length normalized kmerLLR and length for all the counts of a search pool sequence with each the query sequences
for long format: each row of the output dataframe contains a sequence in the search pool fasta
and has five columns: seqName, feature, counts, sum_normLLR, sum_len
seqName corresponds to the header in the search pool fasta file
feature corresponds to the header in the query fasta file
counts is the counts of filtered hit regions of the search pool sequences with the query sequences
sum_normLLR is the sum of length normalized kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences
sum_len is the sum of the length of all counts of a search pool sequence with the query sequences
for wide format: each row of the output dataframe corresponds to a sequence in the search pool fasta
and columns are eachfeature_counts, eachfeature_sum_normLLR, eachfeature_sum_len
the output dataframe can then be used to generalize an overall likeliness of each search pool sequence to all the query sequences


Example:
search all genes on chr16 for the potential hit counts and similarities to the query sequences include mXist repeat A, B, C and E
    $ hmseekr_seqstosummary -qf './fastaFiles/mXist_repeats.fa' -td './transdf.csv' -nf './fastaFiles/all_lncRNA.fa' -pool './fastaFiles/chr16.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -fc 'findhits_condE' -lenf 25 -llrf 0 -pf 1 -name 'seqstosummary_results' -dir './seqstosummary/' -format long -a 'ATCG' -pb


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
    parser.add_argument("-name","--outputname",type=str,help='Desired output name for count file',default='out')
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output count file',default='./')
    
    args = _parse_args_or_exit(parser)

    kmers(
        args.fadir,
        args.kvec,
        args.alphabet,
        args.outputname,
        args.outputdir
        )
    

def console_hmseekr_train():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=TRAIN_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-qd","--querydir",type=str,help='Path to kmer count file for query sequence or sequence of interest (e.g. functional regions of a ncRNA)', required=True)
    parser.add_argument("-nd","--nulldir", type=str,help='Path to kmer count file that compose null model or bakground sequences (e.g. transcriptome, genome, etc.)', required=True)
    parser.add_argument("-k","--kvec",type=str,help='Comma delimited string of possible k-mer values, must be found in the k-mer count file', required=True)
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-qT","--queryT",type=float,help='Probability of query to query transition', default=0.9999)  
    parser.add_argument("-nT","--nullT",type=float,help='Probability of null to null transition', default=0.9999) 
    parser.add_argument("-qPre","--queryPrefix",type=str,help='prefix file name for query', default='query')
    parser.add_argument("-nPre","--nullPrefix",type=str,help='prefix file name for null', default='null')
    parser.add_argument("-dir","--outputdir",type=str,help='path of output directory to save output trained model file',default='./')
    
    args = _parse_args_or_exit(parser)

    train(
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
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output dataframe',default='./')
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-fa","--fasta",action='store_true',help='FLAG: save sequence of hit, ignored if --wt is passed')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    findhits.findhits(
        args.searchpool,
        args.modeldir,
        args.knum,
        args.outputname,
        args.outputdir,
        args.alphabet,
        args.fasta,
        args.progressbar
        )
    

def console_hmseekr_findhits_condE():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=FINDHITS_CONDE_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-pool","--searchpool",type=str,help='Path to fasta file from which the similarity scores to query seq will be calculated and hits seqs be found', required=True)
    parser.add_argument("-m","--modeldir",type=str,help='Path to precalculated model .dict file output from train.py', required=True)
    parser.add_argument("-k","--knum",type=int,help='Value of k to use as an integer. Must be the same as the k value used in training (train function)', required=True)
    parser.add_argument("-name","--outputname",type=str,help='File name for output, useful to include information about the experiment', default='hits')
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output dataframe',default='./')
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-fa","--fasta",action='store_true',help='FLAG: save sequence of hit, ignored if --wt is passed')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    findhits_condE.findhits_condE(
        args.searchpool,
        args.modeldir,
        args.knum,
        args.outputname,
        args.outputdir,
        args.alphabet,
        args.fasta,
        args.progressbar
        )


def console_hmseekr_findhits_nol():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=FINDHITS_NOL_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-pool","--searchpool",type=str,help='Path to fasta file from which the similarity scores to query seq will be calculated and hits seqs be found', required=True)
    parser.add_argument("-m","--modeldir",type=str,help='Path to precalculated model .dict file output from train.py', required=True)
    parser.add_argument("-k","--knum",type=int,help='Value of k to use as an integer. Must be the same as the k value used in training (train function)', required=True)
    parser.add_argument("-name","--outputname",type=str,help='File name for output, useful to include information about the experiment', default='hits')
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output dataframe',default='./')
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-fa","--fasta",action='store_true',help='FLAG: save sequence of hit, ignored if --wt is passed')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    findhits_nol.findhits_nol(
        args.searchpool,
        args.modeldir,
        args.knum,
        args.outputname,
        args.outputdir,
        args.alphabet,
        args.fasta,
        args.progressbar
        )
   


def console_hmseekr_gridsearch():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=GRIDSEARCH_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-qf","--queryfadir",type=str,help='Path to the fasta file of query sequence or sequence of interest (e.g. functional regions of a ncRNA)', required=True)
    parser.add_argument("-nf","--nullfadir", type=str,help='Path to the fasta file of null model (e.g. transcriptome, genome, etc.)', required=True)
    parser.add_argument("-pool","--searchpool",type=str,help='Path to fasta file from which the similarity scores to query seq will be calculated and hits seqs be found', required=True)
    parser.add_argument("-bkgf","--bkgfadir", type=str,help='Path to the fasta file of bakground sequences for seekr normalization vectors(see manual for details)', required=True)
    parser.add_argument("-k","--knum",type=int,help='Value of k to use as an integer. Must be one single integer', required=True)
    parser.add_argument("-ql","--qTlist",type=str,help="Comma delimited string of possible qT (query to query transition) values. For example, 0.1,0.9,0.05 or 0.9,0.99,0.999", required=True)  
    parser.add_argument("-nl","--nTlist",type=str,help="Comma delimited string of possible nT (null to null transition) values. For example, 0.1,0.9,0.05 or 0.9,0.99,0.999", required=True)
    parser.add_argument("-step","--stepaction",action='store_true',help='when called, stepping mode will be applied in generating qT and nT values from qTlist and nTlist (min, max, step); if omitted, qT and nT values will be directly used from qTlist and nTlist')
    parser.add_argument("-fc","--func",type=str,help='which findhits function to use, options are findhits and findhits_condE', default='findhits_condE')
    parser.add_argument("-lf","--lengthfilter",type=int,help='only keep hits sequences that have length > lengthfilter for calculating stats. must be one single integer', default=25)
    parser.add_argument("-name","--outputname",type=str,help='File name for output dataframe', default='gridsearch_results')
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output dataframe and intermediate files',default='./gridsearch/')
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    gridsearch(
        args.queryfadir,
        args.nullfadir,
        args.searchpool,
        args.bkgfadir,
        args.knum,
        args.qTlist,
        args.nTlist,
        args.stepaction,
        args.func,
        args.lengthfilter,
        args.outputname,
        args.outputdir,
        args.alphabet,
        args.progressbar
        )
    

def console_hmseekr_fastarev():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=FASTAREV_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-i","--inputdir", type=str,help='Path to input fasta file', required=True)
    parser.add_argument("-o","--outputdir",type=str,help="Path and name to the output fasta file", required=True)
    
    args = _parse_args_or_exit(parser)

    fastarev(
        args.inputdir,
        args.outputdir
        )
    

def console_hmseekr_genbed():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=GENBED_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-hd","--hitsdir",type=str,help='path to the input hits file, should be the output from findhits with the regular fasta file as input', required=True)
    parser.add_argument("-o","--outputdir",type=str,help='path and name to the output bedfile, do not need to add .bed at the end', required=True)
    parser.add_argument("-len","--lenfilter",type=int,help='the minimum length of the hit, only keep hits that has a length > lenfilter',  default=25)
    parser.add_argument("-llr","--llrfilter",type=float,help='the minimum normalized kmerLLR score, only keep hits that has a normLLR > llrfilter', default=0.5)
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    genbed(
        args.hitsdir,
        args.outputdir,
        args.lenfilter,
        args.llrfilter,
        args.progressbar
        )


def console_hmseekr_genbedrev():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=GENBEDREV_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-hd","--hitsdir",type=str,help='path to the input hits file, should be the output from findhits with the reversed fasta file as input', required=True)
    parser.add_argument("-o","--outputdir",type=str,help='path and name to the output bedfile, do not need to add .bed at the end', required=True)
    parser.add_argument("-len","--lenfilter",type=int,help='the minimum length of the hit, only keep hits that has a length > lenfilter',  default=25)
    parser.add_argument("-llr","--llrfilter",type=float,help='the minimum normalized kmerLLR score, only keep hits that has a normLLR > llrfilter', default=0.5)
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    genbedrev(
        args.hitsdir,
        args.outputdir,
        args.lenfilter,
        args.llrfilter,
        args.progressbar
        )


def console_hmseekr_hitseekr():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=HITSEEKR_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-hd","--hitsdir",type=str,help='Path to the hits .txt file that is the output of findhits function', required=True)
    parser.add_argument("-qf","--queryfadir",type=str,help='Path to the fasta file of query sequence or sequence of interest (e.g. functional regions of a ncRNA)', required=True)
    parser.add_argument("-bkgf","--bkgfadir", type=str,help='Path to the fasta file of bakground sequences for seekr normalization vectors(see manual for details)', required=True)
    parser.add_argument("-k","--knum",type=int,help='Value of k to use as an integer. Must be one single integer', required=True)
    parser.add_argument("-lf","--lengthfilter",type=int,help='only keep hits sequences that have length > lengthfilter for calculating stats. must be one single integer', default=25)
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output dataframe and intermediate files',default='./')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    hitseekr(
        args.hitsdir,
        args.queryfadir,
        args.bkgfadir,
        args.knum,
        args.lengthfilter,
        args.outputdir,
        args.progressbar
        )
    

def console_hmseekr_seqstosummary():
    assert sys.version_info >= (3, 9), "Python version must be 3.9 or higher"
    parser = argparse.ArgumentParser(usage=SEQSTOSUMMARY_DOC, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-qf","--queryfadir",type=str,help='Path to the fasta file of query sequences or sequences of interest (e.g. all repeats of Xist)', required=True)
    parser.add_argument("-td","--transdf",type=str,help='Path to the transition probability dataframe in csv format, columns should be [qT,nT], rows = num of seqs in query fasta', required=True)
    parser.add_argument("-nf","--nullfadir", type=str,help='Path to the fasta file of null model (e.g. transcriptome, genome, etc.)', required=True)
    parser.add_argument("-pool","--searchpool",type=str,help='Path to fasta file from which the similarity scores to query seq will be calculated and hits seqs be found', required=True)
    parser.add_argument("-bkgf","--bkgfadir", type=str,help='Path to the fasta file of bakground sequences for seekr normalization vectors(see manual for details)', required=True)
    parser.add_argument("-k","--knum",type=int,help='Value of k to use as an integer. Must be one single integer', required=True)
    parser.add_argument("-fc","--func",type=str,help='which findhits function to use, options are findhits and findhits_condE', default='findhits_condE')
    parser.add_argument("-lenf","--lenfilter",type=int,help='only keep hits sequences that have length > lenfilter for calculating stats. must be one single integer, default=25', default=25)
    parser.add_argument("-llrf","--llrfilter",type=float,help='only keep hits sequences that have length normalized kmerLLR > llrfilter for calculating stats in the output, default=0', default=0.0)
    parser.add_argument("-pf","--pfilter",type=float,help='only keep hits sequences that have seekr pearson correlation p value < pfilter for calculating stats in the output, default=1', default=1.0)
    parser.add_argument("-name","--outputname",type=str,help='File name for output dataframe', default='seqstosummary_results')
    parser.add_argument("-dir","--outputdir",type=str,help='Directory to save output dataframe and intermediate files',default='./seqstosummary/')
    parser.add_argument("-format","--outdfformat",type=str,help="the format of the output dataframe, default='long', other option is 'wide'",default='long')
    parser.add_argument("-a","--alphabet",type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
    parser.add_argument("-pb","--progressbar",action='store_true',help='when called, progress bar will be shown; if omitted, no progress bar will be shown')

    args = _parse_args_or_exit(parser)

    seqstosummary(
        args.queryfadir,
        args.transdf,
        args.nullfadir,
        args.searchpool,
        args.bkgfadir,
        args.knum,
        args.func,
        args.lenfilter,
        args.llrfilter,
        args.pfilter,
        args.outputname,
        args.outputdir,
        args.outdfformat,
        args.alphabet,
        args.progressbar
        )
    

def _run_console_hmseekr_help(version):
    if version:
        print(__version__)
        sys.exit()

    intro = (
        f"Welcome to hmseekr! ({__version__})\n"
        "Below is a description of all hmseekr commands.\n"
        "For additional help see the README at: \n"
        "https://github.com/CalabreseLab/hmseekr.\n\n"
    )
    print(intro)
    cmds2doc = {
        "hmseekr_kmers": KMERS_DOC,
        "hmseekr_train": TRAIN_DOC,
        "hmseekr_findhits": FINDHITS_DOC,
        "hmseekr_findhits_condE": FINDHITS_CONDE_DOC,
        "hmseekr_findhits_nol": FINDHITS_NOL_DOC,
        "hmseekr_gridsearch": GRIDSEARCH_DOC,
        "hmseekr_fastarev": FASTAREV_DOC,
        "hmseekr_genbed": GENBED_DOC,
        "hmseekr_genbedrev": GENBEDREV_DOC,
        "hmseekr_hitseekr": HITSEEKR_DOC,
        "hmseekr_seqstosummary": SEQSTOSUMMARY_DOC,
    }
    for c, d in cmds2doc.items():
        print(f"{'='*25}\n{c}\n{'='*25}\n{d}")
    conclusion = (
        "To see a full description of flags and defaults, "
        "run any of the commands listed above, without any parameters "
        '(e.g. "$ hmseekr_kmers").'
    )
    print(conclusion)


def console_hmseekr_help():
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--version", action="store_true", help="Print current version and exit.")
    args = parser.parse_args()
    _run_console_hmseekr_help(args.version)
