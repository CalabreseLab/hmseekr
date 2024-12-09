###################################################################################################
### Description: 
# This program use precalculated model (emission matrix, prepared transition matrix, pi and states) from train function
# to find out HMM state path through sequences of interest
# therefore return sequences that has high similarity to the query sequence -- hits sequneces
# this function is different from the basic findhits function in that it calculates the emission probability of the next word given the current word

### Details:
# this function takes in a fasta file which defines the region to search for potential hits (highly similar regions) to the query seq
# also takes in the precalculated model (train function)
# along the searchpool fasta sequences, similarity scores to query seq will be calculated based on the model
# hits segments (highly similar regions) would be reported along with the sequence header from the input fasta file, start and end location of the hit segment
# kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment
# difference from the basic findhits function is that this function calculate the emission probability of the next word given the current word
# as the shift is by 1 nt, the next word has k-1 overlap with the current word
# for example, if the current word is 'TAGC', the next possible words are 'AGCA', 'AGCT', 'AGCC', 'AGCG'
# then the emission probability of the next word ('AGCA') given the current word as 'TAGC' is calculated as 
# np.log2(emission probability of 'AGCA' in the original E) - np.log2(sum of emission probability of all four possible worlds in the original E)


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
# a dataframe containing information about the hits regions: highly similar regions to query seq based on the precalculated model within the input fasta file
# information about the hits regions includes: the sequence header from the input fasta file, start and end location of the hit segment
# kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment if fasta=True

### Example:
# from hmseekr.findhits_condE import findhits_condE

# testhits = findhits_condE(searchpool='../fastaFiles/pool.fa',
#                           modeldir='../markovModels/hmm.dict',
#                           knum=4,
#                           outputname='hits',outputdir='./',
#                           alphabet='ATCG',fasta=True,
#                           progressbar=True)


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
import os


def find_new_starts(vector):

    starts = []  # don't need the first elemment as in viterbi_new it is already taken care of

    for i in range(1, len(vector)):
        # If current element is not consecutive to the previous one, it's the start of a new sequence
        if vector[i] != vector[i - 1] + 1:
            starts.append(i)

    return starts

'''
Viterbi: Calculate the most likely sequence of hidden states given observed sequence, transition matrix, and emission matrix

Inputs: O - list, observed sequence of k-mers
        A - dictionary, transition matrices of hidden states
        E - dictionary, emission matrices of hidden states
        states - list, hidden states (+,-)
        m: + to + transition probability
        n: - to - transition probability
        pi - Dictionary, initial probability of being in + or -
        cE - dictionary of dictionaries, normalized probability of the next word given the current word
Returns:    backTrack - list, sequence of hidden states

'''
def viterbi_new(O,A,E,states,pi,cE,oIdx):

    # Initialize list of dictionaries for the current step
    # and ukprev, which tracks the state transitions that maximize the 'viterbi function'
    uk=[{}]
    ukprev = [{}]
    N = len(O)
    # find new starts in the seq that is the start after N
    # as these position won't be dependent on the previous position
    # we need to use the original E to calculate these positions
    # need to return the position of the new starts in O
    new_starts = find_new_starts(oIdx)

    # calculate initial probabilities in each state given the first kmer
    # use E[state][O[0]] to get the emission probability of the first kmer in each state
    for state in states:
        uk[0][state]=pi[state]+E[state][O[0]]  # uk[0]['+'] = np.log2(.5) +  E['+']['AAAA']
        ukprev[0][state] = None # previous state does not exist, set to None
    # Loop through observed sequence
    # For each state, calculate the cumulative probability recursively
    # Store the state transition that maximizes this probability in ukprev for each current state
    # use cE to get the normalized probability of the current word given the previous word
    for n in range(1,N):
        uk.append({})
        ukprev.append({})
        for state in states:
            prevSelState = states[0] # this is just an arbitrary choice to start checking at the start of the list
            currMaxProb = A[prevSelState][state] + uk[n-1][prevSelState] # probability function
            for pState in states[1:]: # now check the other states...
                currProb = A[pState][state] + uk[n-1][pState]
                if currProb > currMaxProb: # if larger then the first one we checked, we have a new winner, store and continue loop and repeat
                    currMaxProb = currProb
                    prevSelState = pState
            # The emission probability is constant so add at the end rather than in the loop
            # use cE to get the normalized probability of the current word given the previous word if the position is not a new start
            # if it is a new start then calculate using the original E
            if n in new_starts:
                max_prob = currMaxProb + E[state][O[n]]
            else:
                max_prob = currMaxProb + cE[state][O[n-1]][O[n]]
            # save the cumalitive probability for each state
            uk[n][state] = max_prob
            # save the previous state that maximized this probability above
            ukprev[n][state] = prevSelState

    z = max(uk[-1],key=uk[-n].get) # retrieve the state associated with largest log probability
    prev = ukprev[-1][z] # get the state BEFORE "z" above that maximized z
    backtrack = [z,prev] # start our backtrack with knowns
    # Loop through BACKWARDS, getting the previous state that yielded the 'current' max state
    for n in range(N-2,-1,-1):
        backtrack.append(ukprev[n+1][prev]) # n+1 is the "current" state as we're going backwards, ukprev[n+1][prev] returns the state that maximized
        prev = ukprev[n+1][prev]
    backtrack = backtrack[::-1] # reverse the order
    return backtrack


''' LLR
Return log-likelihood ratio between two models in HMM for + k-mers
Input: sequnce of hits, value of k, k-mer frequencies in HMM emmission matrix and conditioned emmission matrix
Output: Array of LLRs for each hit

hits = seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
'''
def LLR_new(hits,k,E,cE):
    arr = np.zeros(len(hits))
    for i,hit in enumerate(hits):
        LLRPos,LLRNeg=0,0
        # calculate the first kmer based on E
        kmer=hit[0:k]
        LLRPos += E['+'][kmer]
        LLRNeg += E['-'][kmer]
        # calculate the rest of the kmers based on cE
        for j in range(1,len(hit)-k+1):
            prevkmer=hit[j-1:j-1+k]
            kmer=hit[j:j+k]
            LLRPos += cE['+'][prevkmer][kmer]
            LLRNeg += cE['-'][prevkmer][kmer]
        llr = LLRPos-LLRNeg
        arr[i] = llr
    return arr


'''
Combine all the data to create a dataframe
E: emission matrix
seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
'''

def hitOutput_new(seqHits,starts,ends,k,E,tHead,cE):
    info = list(zip(seqHits,starts,ends)) # example [('GGCCCGGTGTGGTCGGCCTCATTTTGGAT.......', 88, 177),......]
    dataDict = dict(zip(list(range(len(seqHits))),info))
    df = pd.DataFrame.from_dict(dataDict,orient='index')
    #calculate log-likelihood ratio of k-mers in the + model vs - model
    df['kmerLLR'] = LLR_new(seqHits,k,E,cE)
    df['seqName'] = tHead
    df.columns = ['Sequence','Start','End','kmerLLR','seqName']
    df.sort_values(by='kmerLLR',inplace=True,ascending=False)
    df.reset_index(inplace=True)
    fa = df['Sequence']
    df = df[['Start','End','kmerLLR','seqName','Sequence']]

    return df



# if in groupedHits, there is '+' string that is greater than streaklen, and '-' string that is less than gaplen
# for '+' string that is greater than streaklen, check the following '-' string
# if the '-' string is less than gaplen, change the '-' string to '+'
# merge the '+' streak, the converted '-' string, and the following '+' string
# repeat the process until no more '-' string can be converted

# def process_grouped_hits(groupedHits, streaklen, gaplen):
#     # Make a copy of the input list to avoid modifying the original
#     groupedHits_copy = groupedHits.copy()
    
#     # Check for any '+' string longer than streaklen
#     has_long_plus = any(s.startswith('+') and len(s) > streaklen for s in groupedHits_copy)
#     # Check for any '-' string shorter than gaplen
#     has_short_minus = any(s.startswith('-') and len(s) < gaplen for s in groupedHits_copy)

#     # If either condition is not met, return the copied list as is
#     if not (has_long_plus and has_short_minus):
#         return groupedHits_copy

#     changed = True
#     while changed:
#         changed = False
#         i = 0
#         while i < len(groupedHits_copy) - 1:
#             s_current = groupedHits_copy[i]
#             s_next = groupedHits_copy[i+1]

#             if s_current.startswith('+') and len(s_current) > streaklen:
#                 if s_next.startswith('-') and len(s_next) < gaplen:
#                     # Check if there is a next '+' string to merge
#                     if i+2 < len(groupedHits_copy):
#                         s_next_next = groupedHits_copy[i+2]
#                         if s_next_next.startswith('+'):
#                             # Convert the short '-' string to '+' and merge
#                             s_next_converted = '+' * len(s_next)
#                             new_plus = s_current + s_next_converted + s_next_next
#                             groupedHits_copy[i] = new_plus
#                             # Remove the next two entries
#                             del groupedHits_copy[i+1:i+3]
#                             changed = True
#                             continue
#                 # Move to the next index only if no changes were made
#                 i += 1
#             else:
#                 i += 1
#     return groupedHits_copy




def hmmCalc_new(tHead,tSeq,hmm,k,alphabet):
    #tHead,tSeq = data
    O,oIdx,nBP = corefunctions.kmersWithAmbigIndex(tSeq,k)
    A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
    cE = condition_E(E,alphabet)
    bTrack = viterbi_new(O,A,E,states,pi,cE,oIdx)
    #Zip the indices of unambig k-mers with their viterbi derived HMM state labels
    coordBTrack = list(zip(oIdx,bTrack)) # [(1,'-'),(2,'+',...(n,'+'))]
    mergedTrack = coordBTrack + nBP # add back in ambig locations
    mergedTrack.sort(key=itemgetter(0)) # sort master list by index
    hmmTrack = [i[1] for i in mergedTrack] # fetch just state label from mergedTrack ['-','+',...,'+']
    groupedHits = corefunctions.groupHMM(hmmTrack) # ['-----','++++++++++','-','++++','------------']
    # merge the gaps inbwteen the '+' streaks
    #groupedHits = process_grouped_hits(groupedHits_raw, streaklen, gaplen)

    # Return sequences of HMM hits, and their start and end locations in the original sequence
    seqHits,starts,ends = corefunctions.formatHits(groupedHits,k,tSeq)
    if (seqHits):
        df = hitOutput_new(seqHits,starts,ends,k,E,tHead,cE)
        return tHead,df
    # Alternative output (transcript by transcript)

    else:
        return tHead,None


def condition_E(E,alphabet):
    # as the shift is by 1 nt, the next word has k-1 overlap with the current word
    # so the next word (AGC(ATGC)) probability is dependent on the current word (TAGC)
    # generate a dictionary of a dictionary to store the normailized probability of the next word given the current word
    # E['+']['TAGC'] = {'AGCA':0.25,'AGCT':0.25,'AGCC':0.25,'AGCG':0.25}

    # initialize the dictionary
    cE = {state: {kmer: {} for kmer in E[state].keys()} for state in E.keys()}

    for state in E.keys():
        for kmer in E[state].keys():
            # get the last k-1 nt of the kmer
            # this is the overlapping part between the current kmer and the next kmer
            prefix = kmer[1:]
            # generate all possible next kmers based on alphabet
            # for example, if prefix is ATC and alphabet is ATCG, then the next kmers are ATCA,ATCC,ATCG,ATCT
            next_kmers = [prefix + suffix for suffix in alphabet]
            # get the sum of the probability of the next kmers
            # convert the probability from log2 back for summing up
            # then convert back to log2
            sum_prob = np.log2(sum([2**E[state][kmer] for kmer in next_kmers]))

            # normalize the probability of the next kmers as np.log2(prob)-np.log2(sum_prob)
            # adding the normalized probability to the dictionary
            for next_kmer in next_kmers:
                cE[state][kmer][next_kmer] = E[state][next_kmer] - sum_prob

    return cE




def findhits_condE(searchpool,modeldir,knum,outputname='hits',outputdir='./',alphabet='ATCG',fasta=True,progressbar=True):

    #Loop over values of k
    alphabet = alphabet.upper()

    hmm = pickle.load(open(modeldir,'rb'))
    A,E,pi,states = hmm['A'],hmm['E'],hmm['pi'],hmm['states']

    # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
    k = int(log(len(hmm['E']['+'].keys()),len(alphabet)))
    assert k == knum, 'Value of k provided does not match supplied hmm file'

    cookedFasta = corefunctions.getCookedFasta(searchpool)
    targetSeqs,targetHeaders = cookedFasta[1::2],cookedFasta[::2]


    #Pool processes onto number of CPU cores specified by the user
    # with pool.Pool(args.n) as multiN:
    #     jobs = multiN.starmap(hmmCalc,product(*[list(zip(targetHeaders,targetSeqs))]))
    #     dataDict = dict(jobs)

    #dataDict = dict(starmap(hmmCalc,product(*[list(zip(targetHeaders,targetSeqs))])))

    data_pairs = zip(targetHeaders, targetSeqs)

    total_items = len(targetHeaders)

    dataDict = {}

    # set iterable based on progress_bar
    iterable = tqdm(data_pairs,total=total_items) if progressbar else data_pairs

    # Loop through data_pairs with a progress bar
    for header, seq in iterable:
        # Call hmmCalc and get returned header and value
        theader, value = hmmCalc_new(header, seq, hmm, k, alphabet)
        
        # Assign the value to the corresponding header in your dictionary
        dataDict[theader] = value

    # dataDict = dict(hmmCalc(header, seq, hmm, k) for header, seq in data_pairs)


    #Check if no hits were found
    if not all(v is None for v in dataDict.values()):
        dataFrames = pd.concat([df for df in dataDict.values() if not None])
        dataFrames['Start']+=1 #1-start coordinates
        dataFrames['End']
        dataFrames['Length'] = dataFrames['End'] - dataFrames['Start'] +1
        dataFrames = dataFrames[['Start','End','Length','kmerLLR','seqName','Sequence']]
        if not fasta:
            dataFrames = dataFrames[['Start','End','Length','kmerLLR','seqName']]
        dataFrames.sort_values(by='kmerLLR',ascending=False,inplace=True)
        dataFrames.reset_index(inplace=True,drop=True)
    else:
        dataFrames = pd.DataFrame(columns=['Start', 'End', 'Length', 'kmerLLR', 'seqName', 'Sequence'])
        if not fasta:
            dataFrames = dataFrames[['Start', 'End', 'Length', 'kmerLLR', 'seqName']]

    mDir = outputdir
    if not mDir.endswith('/'):
        mDir+='/'
    
    if not os.path.exists(mDir):
        os.mkdir(mDir)

    dataFrames.to_csv(f'{mDir}{outputname}_{k}_viterbi.txt',sep='\t', index=False)

    return dataFrames
