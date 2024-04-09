###################################################################################################
### Description: 
# This function firstly filter the hits output from findhits.py 
# and then generate a bed file for the filtered hits
# this function only applies to the output from findhits.py with the regular fasta file as input
# for reversed fasta file, please use genbedrev.py

### Details:
# this function takes in the output from findhits.py and filter the hits based on
# the hit length and normalized kmerLLR score
# kmerLLR is the log likelihood ratio of of the probability 
# that the set of k-mers y within a hit derived from the QUERY versus the NULL state
# it is the sum of the log2(Q/N) ratio for each kmer within a hit
# the normalized kmerLLR (normLLR) is the kmerLLR divided by the hit length
# so the normLLR is the log2 of the Q/N ratio, i.e. if set normLLR > 0.5, the Q/N ratio is ~ 1.41
# then all the hits after the filtering will be converted to a bedfile

### Input:
# hitsdir: path to the input hits file, should be the output from findhits with the regular fasta file as input
# please use genbedrev.py for the reversed fasta file 
# outputdir: path and name to the output bed file, do not need to add .bed at the end
# lenfilter: the minimum length of the hit, only keep hits that has a length > lenfilter, default is 25
# llrfilter: the minimum normalized kmerLLR score, only keep hits that has a normLLR > llrfilter, default is 0.5
# progressbar: whether to show the progress bar, default is True

### Output:
# a bedfile with the following columns: 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'
# for all hits that passed the filtering
# score is the normLLR score, name is 'fwd' as it should only be applied for the regular fasta file
# for reversed fasta file, please use genbedrev.py

### Example:
# from hmseekr.genbed import genbed

# testbed = genbed(hitsdir='../mm10expmap_queryA_4_viterbi.txt', 
#                  outputdir='../mm10expmap_queryA_4_viterbi',
#                  lenfilter=25, llrfilter=0.5,progressbar=True)


########################################################################################################


import pandas as pd
from tqdm import tqdm

def genbed(hitsdir, outputdir, lenfilter=25, llrfilter=0.5, progressbar=True):

    print('Make sure the input is results of findhits using the regular fasta file!')
    # Read the table into a pandas DataFrame
    df = pd.read_csv(hitsdir, sep='\t', header=0)
    # Keep only the first 5 columns
    df = df.iloc[:, :5]

    # Filter rows based on Length
    df = df[df['Length'] > lenfilter]

    # check if there are still hits after filtering
    if len(df) == 0:
        print('No hits after length filtering')
        print('Please try a lower lenfilter value')
        return None
    
    else:
        # Filter based on normalized LLR 
        df['score'] = df['kmerLLR'] / df['Length']
        df = df[df['score'] > llrfilter]

        # check if there are still hits after filtering
        if len(df) == 0:
            print('No hits after normLLR filtering')
            print('Please try a lower llrfilter value')
            return None
        else: 

            # Process the seqName column
            temp = df['seqName'].str.replace('>', '', regex=False)
            temp = temp.str.replace(')', '', regex=False)
            temp = temp.str.split(':', expand=False)

            df['chrom'] = [x[0] for x in temp]
            temp = [x[1] for x in temp]

            # Split strand and update df
            # first split strand as it is the same as the - in the middle of the num
            temp = [x.split('(', 1) for x in temp]
            df['strand'] = [x[1] for x in temp]
            temp = [x[0] for x in temp]

            temp = [x.split('-', 1) for x in temp]
            df['seqStart'] = pd.to_numeric([x[0] for x in temp])
            df['seqEnd'] = pd.to_numeric([x[1] for x in temp])

            df['chromStart'] = None
            df['chromEnd'] = None

            # use progress bar


            iterable = tqdm(range(len(df))) if progressbar else range(len(df))

            for n in iterable:
                # integers in python are default to be standard notation (not scientific notation)

                if df['strand'].iloc[n] == '+':
                    df.iat[n,df.columns.get_loc('chromStart')] = int(df['seqStart'].iloc[n] + df['Start'].iloc[n])
                    df.iat[n,df.columns.get_loc('chromEnd')] = int(df['chromStart'].iloc[n] + df['Length'].iloc[n] - 1)
                else:
                    df.iat[n,df.columns.get_loc('chromStart')] = int(df['seqEnd'].iloc[n] - df['End'].iloc[n] + 1)
                    df.iat[n,df.columns.get_loc('chromEnd')] = int(df['chromStart'].iloc[n] + df['Length'].iloc[n] - 1)


            # Add 'name' columns
            df['name'] = 'fwd'

            # Reorder the columns
            df = df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]

            df.to_csv(f'{outputdir}.bed', sep='\t', index=False, header=False)

            return df

