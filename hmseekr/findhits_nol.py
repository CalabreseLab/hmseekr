###################################################################################################
### Description: 
# This program use precalculated model (emission matrix, prepared transition matrix, pi and states) from train function
# to find out HMM state path through sequences of interest
# therefore return sequences that has high similarity to the query sequence -- hits sequneces
# this function is different from the basic findhits in that it uses non-overlapping kmers and shifts the sequence by 1nt each time until the kmer size is reached
# this reduces the kmer dependency of overlapping kmers

### Details:
# this function takes in a fasta file which defines the region to search for potential hits (highly similar regions) to the query seq
# also takes in the precalculated model (train function)
# along the searchpool fasta sequences, similarity scores to query seq will be calculated based on the model
# hits segments (highly similar regions) would be reported along with the sequence header from the input fasta file, start and end location of the hit segment
# kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment
# different from the basic findhits, this function scan along the sequences with non-overlapping kmers
# if sequence is 'ATGCTTTTGCGC' the kmers would be 'ATGC','TTTT','GCGC'
# then it chops off 1nt from the start of each sequence in the searchpool and re-scan with non-overlapping kmers
# so the sequence would be 'TGCTTTTGCGC' and the kmers would be 'TGCT','TTTG'
# the last bit of sequences that is less than the kmer size would be discarded
# this is done until the kmer size is reached
# in this way, it generates k result txt files each starts at different position of the sequence with non-overlapping kmers
# this function reduces the kmer dependency of overlapping kmers in findhits, but also considered the different combinations of kmers at different start positions
# the outcome txt files can be further process to find the best hit regions
# for example, only use the regions that are hits in all k result txt files


### Input:
# searchpool: Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
# modeldir: Path to precalculated model .dict file output from train.py'
# knum: Value of k to use as an integer. Must be the same as the k value used in training (train function) that produced the model
# outputname: File name for output, useful to include information about the experiment, default='hits'
# outputdir: path of output directory to save output dataframe, default is current directory
# alphabet: String, Alphabet to generate k-mers default='ATCG'
# fasta: whether to save sequence of hit in the output dataframe, default=True: save the actual sequences
# progressbar: whether to show progress bar, default=True: show progress bar

### Output:
# k dataframes containing information about the hits regions: highly similar regions to query seq based on the precalculated model within the input fasta file
# information about the hits regions includes: the sequence header from the input fasta file, start and end location of the hit segment
# kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment if fasta=True

### Example:
# from hmseekr.findhits_nol import findhits_nol

# testhits = findhits_nol(searchpool='../fastaFiles/pool.fa',
#                         modeldir='../markovModels/hmm.dict',
#                         knum=4,outputname='hits',outputdir='./',
#                         alphabet='ATCG',fasta=True,
#                         progressbar=True)


########################################################################################################


'''
-------------------------------------------------------------------------------------------------------------------------------------------------
alphabet: list
    a list of base pairs. for example ['A', 'T', 'C', 'G']
model: str
    Path to .dict file output from train.py 
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
hmm: dict
    Format example:
    {
    'A':{'+': {'+': -0.00014427671804501932, '-': -13.287712379549609}, '-': {'+': -13.287712379549609, '-': -0.00014427671804501932}},
    'E':{'+': {'AAAA': -5.082989364671032, 'AAAT': -8.330916878114618, 'AAAC': -7.330916878114617, ....}, '-': {'AAAA': -6.735347642028761, 'AAAT': -7.242499465916189, 'AAAC': ...}}
    'pi':{'+': -1.0, '-': -1.0},
    'states':('+', '-')
    }
O: list
    observed sequence of k-mers. (doesn't have 'N' in kmers)
kmers: list
    a list of specific k length kmers
    Format example:
    k=3 ['AAA', 'AAT', 'AAC', 'AAG', 'ATA'....]
    k=4 ['AAAA', 'AAAT', 'AAAC', 'AAAG', 'AATA',....]
target: Reader object
    an object of Reader class from package seekr
targetSeqs: list
    a list of sequence reads without header
    Format example:
    ['AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCC', 'CTCCGCGTGGTCTATGATGGTGCATTTTGGTCCAGTCAGGCCCGGTGTGG', 'TCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGGCGTGGTTACG',....]
targetHeaders: list
    a list of headers
    Format example
    ['>ENSMUST00000193812.1|ENSMUSG00000102693.1|OTTMUSG00000049935.1|OTTMUST00000127109.1|RP23-271O17.1-001|RP23-271O17.1|1070|', '>ENSMUST00000195335.1|ENSMUSG00000103377.1|OTTMUSG00000049960.1|OTTMUST00000127145.1|RP23-317L18.4-001|RP23-317L18.4|2819|'......]
dataDict: dict
    key is a sequence header
    value is a dataframe
    Format example:
    {'>mm10kcnq1ot1':     Start    End    kmerLLR        seqName                                           Sequence
    0     835    999  84.173049  >mm10kcnq1ot1  ATTCGTGCCGCGCTTTCGCGGCTGGGCTCCATCTTCGTTTTGCCGC...
    1   79184  79257  72.086333  >mm10kcnq1ot1  ATTATTTTGTGTCTTTTTTTGTTTGTTTGTTTTTTGTTTTTTGTTT...
    2   69556  69605  63.331314  >mm10kcnq1ot1  TCCCAACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA.....}
O: list
    a list of k length kmers extracted from sequence and exclude all kmers with 'N' inside.
    Format example:
    when k=4, O is something like ['GGAC', 'GACA', 'ACAG', 'CAGC', 'AGCA',...]
oIdx: list
    a list of index. These index is the real location of O's kmers in the original sequence
    Format example:
    [0, 1, 2, 3, 6, 7, 8, 9...]
nBP: list
    a list of tuples.
    Format example:
    [(36, 'N'), (111, 'N'), (300, 'N').....]
bTrack: list
    a list of query/null states
    Format example:
    [-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',....]
hmmCalc
Run several functions including viterbi algorithm, log-likelihood, and generate output dataframes
Input: fasta file information
Output: dataframe object

'''



from hmseekr import corefunctions
#from itertools import product
#from itertools import starmap

import pickle
from math import log
import pandas as pd
from operator import itemgetter
from tqdm import tqdm
import numpy as np



''' kmersWithAmbigIndex
Return list of kmers, indices of k-mers without ambiguity, and indices of those
with ambiguity
Input: string
Output: List of string, list of indices, list of indices
'''
def kmersWithAmbigIndex_nol(tSeq,k):
    O = [tSeq[i:i+k].upper() for i in range(0, len(tSeq) - k + 1, k)]
    O = [o for o in O if 'N' not in o]  # example: when k=4, O is something like ['GGAC', 'GACA', 'ACAG', 'CAGC', 'AGCA',...]
    # Match k-mers without ambig char to index in original string
    oIdx = [i for i in range(0,len(tSeq)-k+1,k) if 'N' not in tSeq[i:i+k]] # example: when k=4, oIdx is something like [0, 1, 2, 3, 4, 5, 6, 7, 8, 9...]
    # Match k-mers with ambig char to index in original string
    nBP = [i for i in range(0,len(tSeq)-k+1,k) if 'N' in tSeq[i:i+k]]
    # zip the indices with marker character N
    nBP = list(zip(nBP,['N']*len(nBP))) # example. nBP is something like [(36, 'N'), (111, 'N'), (300, 'N').....]
    return O, oIdx, nBP


'''
Process raw results to find out the break point of the original sequence then break it to corresponding query reads or null reads fragments based on groupedHits. 
Also the start position and end position of each fragments in the original sequence are recorded.

From:
groupedHits: ['-----','++++++++++','-','++++','------------']
tSeq: 'AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCCCTCCGCGTGG......'

to:
seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
starts: [   88   498   835  5614 14919 15366 19473 22753 69402 69556.....]
ends: [  177   627   999  5636 14971 15388 19502 22776 69436 69605.......]
'''
def formatHits_nol(groupedHits,k,tSeq):
    idx = 0
    indexGroupHits = []
    # Loop below formats the hmm output as such:
    # [([0,1,2]),'---'),([3,4],'++'),([5],'-'),...]
    # Grouping HMM states with their correct index in the list of k-mers
    for i,group in enumerate(groupedHits):   # groupedHits example ['-----','++++++++++','-','++++','------------']
        indexGroupHits.append([])
        for kmer in group:
            indexGroupHits[i].append(idx)
            idx+=1
    hits = list(zip(indexGroupHits,groupedHits)) # hits example [([0,1,2]),'---'),([3,4],'++'),([5],'-'),...]
    
    seqHits = []
    seqHitCoords = []
    for group in hits:
        if '+' in group[1]:
            start,end = group[0][0],group[0][-1]+1 #convert k-mer coord to bp coord
            seqHitCoords.append(f'{start}:{end}')
            seqHits.append(tSeq[start:end])
    # seqHits - squence hits base pair version - example: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
    # tSeq example 'AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCCCTCCGCGTGG......'
    # seqHitCoords example ['88:177', '498:627', '835:999', '5614:5636', '14919:14971', '15366:15388', '19473:19502', '22753:22776', '69402:69436', '69556:69605'......]
    starts = np.array([int(c.split(':')[0]) for c in seqHitCoords]) # starts example [   88   498   835  5614 14919 15366 19473 22753 69402 69556.....]
    ends = np.array([int(c.split(':')[1]) for c in seqHitCoords]) # ends example [  177   627   999  5636 14971 15388 19502 22776 69436 69605.......]
    return seqHits,starts,ends


''' LLR
calculate LLR of non-overlapping k-mers
Return log-likelihood ratio between two models in HMM for + k-mers
Input: sequnce of hits, value of k, k-mer frequencies in HMM emmission matrix
Output: Array of LLRs for each hit

hits = seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
'''
def LLR_nol(hits,k,E):
    arr = np.zeros(len(hits))
    for i,hit in enumerate(hits):
        LLRPos,LLRNeg=0,0
        for j in range(0,len(hit)-k+1,k):
            kmer=hit[j:j+k]
            LLRPos += E['+'][kmer]
            LLRNeg += E['-'][kmer]
        llr = LLRPos-LLRNeg
        arr[i] = llr
    return arr

'''
Combine all the data to create a dataframe
E: emission matrix
seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
'''

def hitOutput_nol(seqHits,starts,ends,k,E,tHead,tSeq):
    info = list(zip(seqHits,starts,ends)) # example [('GGCCCGGTGTGGTCGGCCTCATTTTGGAT.......', 88, 177),......]
    dataDict = dict(zip(list(range(len(seqHits))),info))
    df = pd.DataFrame.from_dict(dataDict,orient='index')
    #calculate log-likelihood ratio of k-mers in the + model vs - model
    df['kmerLLR'] = LLR_nol(seqHits,k,E)
    df['seqName'] = tHead
    df.columns = ['Sequence','Start','End','kmerLLR','seqName']
    df.sort_values(by='kmerLLR',inplace=True,ascending=False)
    df.reset_index(inplace=True)
    fa = df['Sequence']
    df = df[['Start','End','kmerLLR','seqName','Sequence']]

    return df



def hmmCalc_nol(tHead,tSeq,hmm,k):
    #tHead,tSeq = data
    O,oIdx,nBP = kmersWithAmbigIndex_nol(tSeq,k)
    A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
    bTrack = corefunctions.viterbi(O,A,E,states,pi)
    #Zip the indices of unambig k-mers with their viterbi derived HMM state labels
    coordBTrack = list(zip(oIdx,bTrack)) # [(1,'-'),(2,'+',...(n,'+'))]
    mergedTrack = coordBTrack + nBP # add back in ambig locations
    mergedTrack.sort(key=itemgetter(0)) # sort master list by index
    hmmTrack = [i[1] for i in mergedTrack] # fetch just state label from mergedTrack ['-','+',...,'+']
    # as each position is actually the kmer start position, we expand the state labels to the length of the kmer
    hmmTrack = [hmmT * k for hmmT in hmmTrack]
    groupedHits = corefunctions.groupHMM(hmmTrack) # ['-----','++++++++++','-','++++','------------']

    # Return sequences of HMM hits, and their start and end locations in the original sequence
    seqHits,starts,ends = formatHits_nol(groupedHits,k,tSeq)
    if (seqHits):
        df = hitOutput_nol(seqHits,starts,ends,k,E,tHead,tSeq)
        return tHead,df
    # Alternative output (transcript by transcript)

    else:
        return tHead,None





def findhits_nol(searchpool,modeldir,knum,outputname='hits',outputdir='./',alphabet='ATCG',fasta=True,progressbar=True):

    #Loop over values of k
    alphabet = alphabet.upper()

    hmm = pickle.load(open(modeldir,'rb'))
    A,E,pi,states = hmm['A'],hmm['E'],hmm['pi'],hmm['states']

    # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
    k = int(log(len(hmm['E']['+'].keys()),len(alphabet)))
    assert k == knum, 'Value of k provided does not match supplied hmm file'

    cookedFasta = corefunctions.getCookedFasta(searchpool)
    targetSeqs_ori,targetHeaders = cookedFasta[1::2],cookedFasta[::2]

    #Pool processes onto number of CPU cores specified by the user
    # with pool.Pool(args.n) as multiN:
    #     jobs = multiN.starmap(hmmCalc,product(*[list(zip(targetHeaders,targetSeqs))]))
    #     dataDict = dict(jobs)

    #dataDict = dict(starmap(hmmCalc,product(*[list(zip(targetHeaders,targetSeqs))])))

    # loop thru kmer size to generate targetSeqs and run findhits
    # each loop will remove 1nt from the start of each sequence in targetSeqs until the kmer size is reached
    # this will generate k different sequences of non-overlapping kmers
    for i in range(k):
        targetSeqs = [seq[i:] for seq in targetSeqs_ori]

        data_pairs = zip(targetHeaders, targetSeqs)

        total_items = len(targetHeaders)

        dataDict = {}

        # set iterable based on progress_bar
        iterable = tqdm(data_pairs,total=total_items) if progressbar else data_pairs

        # Loop through data_pairs with a progress bar
        for header, seq in iterable:
            # Call hmmCalc and get returned header and value
            theader, value = hmmCalc_nol(header, seq, hmm, k)
            
            # Assign the value to the corresponding header in your dictionary
            dataDict[theader] = value

        # dataDict = dict(hmmCalc(header, seq, hmm, k) for header, seq in data_pairs)


        #Check if no hits were found
        # if not all(v == None for v in dataDict.values()):
        dataFrames = pd.concat([df for df in dataDict.values() if not None])
        dataFrames['Start'] += i+1 #1-start coordinates and count in the shift
        dataFrames['End'] += i # count in the shift
        # this make sure the coordinates across the shifts are consistent with the original sequence
        dataFrames['Length'] = dataFrames['End'] - dataFrames['Start'] +1
        dataFrames = dataFrames[['Start','End','Length','kmerLLR','seqName','Sequence']]
        if not fasta:
            dataFrames = dataFrames[['Start','End','Length','kmerLLR','seqName']]
        dataFrames.sort_values(by='kmerLLR',ascending=False,inplace=True)
        dataFrames.reset_index(inplace=True,drop=True)

        mDir = outputdir
        if not mDir.endswith('/'):
            mDir+='/'

        dataFrames.to_csv(f'{mDir}{outputname}_{k}_shift{i}_viterbi.txt',sep='\t', index=False)

