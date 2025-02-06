###################################################################################################
### Description: 
# This program calculate the emission matrix based on kmer counts 
# Prepare transition matrix, states and saves them to a binary file 
# that countains a dictionary of dictionaries

### Details:
# this function takes in kmer count file for sequences of interest (e.g. query seq, functional regions of a ncRNA)
# and kmer count file for background sequences (e.g. null seq, transcriptome, genome)
# these kmer count files are generated using the kmers function
# with the k specified (which should be calculted in the kmers function)
# and query to query transition rate (qT) and null to null transition rate (nT) specified
# it calculates the emission matrix, transition matrix, states and saves them to a binary file

### Input:
# querydir: Path to kmer count file for sequence of interest or query seq (e.g. functional regions of a ncRNA)
# nulldir: Path to kmer count file that compose null model/background sequences (e.g. transcriptome, genome, etc.)
# kvec: Comma delimited string of possible k-mer values. For example, '3,4,5' or just '4'
# numbers in kvec must be found in the k-mer count file (precalculated by kmers function)
# alphabet: String, Alphabet to generate k-mers, default=ATCG
# queryT: Probability of query to query transition, default=0.99, should be between 0 and 1 but not equal to 1 or 0
# nullT: Probability of null to null transition, default=0.93, should be between 0 and 1 but not equal to 1 or 0
# queryPrefix: prefix file name for query, defualt='query'
# nullPrefix: prefix file name for null, defualt='null'
# outputdir: path of output directory to save output trained model file in .dict format, default is current directory


### Output:
# a dictionary containing models for each value of k specified
# it saves Hidden state transition matrix, Hidden state emission matrix, states and Starting probability of each hidden state
# the output is saved as a binary file with .dict extension

### Example:
# from hmseekr.train import train

# testmodel = train(querydir='./counts/repeatA.dict', 
#                   nulldir='./counts/all_lncRNA.dict',
#                   kvec='4', alphabet='ATCG', 
#                   queryT=0.9999, nullT=0.9999,
#                   queryPrefix='repeatA', 
#                   nullPrefix='lncNRA', outputdir='./')


########################################################################################################


'''

------------------------------------------------------------------------
query: str
    query sequence. Represented by '+' in certain code
null: str
    null sequence. Represented by '-' in certain code
qT: float
    transition rate from query to query so transition rate from query to null would be 1 - qT
nT: float
    transition rate from null to null so transition rate from null to query would be 1 - nT
alphabet: list
    a list of base pairs. for example ['A', 'T', 'C', 'G']
qCount: dict
    query sequence reads' kmer count. 
    data structure example: 
    {2: {'AA': 30, 'AT': 24, 'AC': 18,.....}, 3: {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....}, 4: {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}}
nCount: dict
    null sequence reads' kmer count.
    data structure example:
    {2: {'AA': 30, 'AT': 24, 'AC': 18,.....}, 3: {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....}, 4: {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}}
kVals: list
    a list of k values - for example [2, 3, 4]
qKCount: dict
    query sequence's specific k's kmer counts dictionary
    data structure example:
    {'AA': 30, 'AT': 24, 'AC': 18,.....} or {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....} or {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}
nKCount: dict
    null sequence's specific k's kmer counts dictionary
    data structure example:
    {'AA': 30, 'AT': 24, 'AC': 18,.....} or {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....} or {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}
A: dict
    Hidden state transition matrix
    initial value is: 
    {'+':{'+':np.log2(args.qT),'-':np.log2(1-args.qT)},'-':{'+':np.log2(1-args.nT),'-':np.log2(args.nT)}}
E: dict
    Hidden state emission matrix
    Format example:
    {'+': {'AAAA': -5.082989364671032, 'AAAT': -8.330916878114618, 'AAAC': -7.330916878114617,.....'GGGC': -6.523561956057013, 'GGGG': -6.523561956057013}, '-': {'AAAA': -6.735347642028761, 'AAAT': -7.242499465916189, 'AAAC': ...}}
states: dict
    states = ('+','-')
pi: dict
    Starting probability of each hidden state (+ or -)
    initial value is: 
    {'+':np.log2(.5),'-':np.log2(.5)} 
kmers: list
    a list of specific k length kmers
    Format example:
    k=3 ['AAA', 'AAT', 'AAC', 'AAG', 'ATA'....]
    k=4 ['AAAA', 'AAAT', 'AAAC', 'AAAG', 'AATA',....]
'''


from hmseekr import corefunctions
import os
import pickle

def train(querydir, nulldir, kvec, alphabet='ATCG', queryT=0.99, nullT=0.93, queryPrefix='query', nullPrefix='null', outputdir='./'):


    kVals = [int(i) for i in kvec.split(',')]

    # Check if specified directory exists
    # create new directory if not existing
    if not outputdir.endswith('/'):
        outputdir+='/'
    newDir = f'{outputdir}{queryPrefix}_{nullPrefix}/'

    if not os.path.exists(newDir):
        os.mkdir(newDir)

    #alphabet = list(alphabet)  # like ['A', 'T', 'C', 'G']

    # Load k-mer counts
    qCount = pickle.load(open(querydir,'rb'))
    nCount = pickle.load(open(nulldir,'rb'))


    # Loop through specified values of k
    # Check if they exist in the counts file,
    # and call corefunctions.HMM to generate the HMM matrices
    for k in kVals:
        if (k in qCount.keys()) and (k in nCount.keys()):
            qKCount = qCount[k]
            nKCount = nCount[k]
            kDir = newDir+f'{k}/'
            if not os.path.exists(kDir):
                os.mkdir(kDir)
            A,E,states,pi = corefunctions.HMM(qKCount,nKCount,k,alphabet,queryT,nullT)
            # kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
            # queryMkv = corefunctions.transitionMatrix(qKCount,k,alphabet)
            # nullMkv = corefunctions.transitionMatrix(nKCount,k,alphabet)
            # lgTbl = corefunctions.logLTbl(queryMkv,nullMkv)
        else:
            print(f'Missing {k}-mer counts in count file... skipping')

        # np.savetxt(f'{kDir}logtbl.mkv',lgTbl)
        pickle.dump({'A':A,'E':E,'pi':pi,'states':states},open(f'{kDir}hmm.dict','wb'))
