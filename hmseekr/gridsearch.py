###################################################################################################
### Description: 
# This function performs a grid search to find the best trasnition probabilities for query to query and null to null states
# which is used in the train function

### Details:
# this function takes in query sequence, null squence and background sequence (see Inputs for details) fasta files
# ranges and steps for query to query transition rate (qT) and null to null transition rate (nT)
# a specific kmer number and performs the train function and findhits function for each combination of qT and nT
# within the hits sequences (findhits function results), only keep the sequence with length greater than lengthfilter 
# then calculates pearson correlation r score (seekr.pearson) between the filtered hit sequences (findhits results) and the query sequence
# it returns a dataframe (.csv file) containing the qT, nT, kmer number, the total number of hits sequences and the mean, median, standard deviation of the hits sequences' pearson correlation r score to the query sequence
# and the mean, median, standard deviation of the length of the hits sequences
# if query fasta contains more than one sequence, all the sequences in query fasta file will be merged to one sequence 
# for calculating kmer count files for hmseekr (hmseekr.kmers) and for calculating seekr.pearson 
# this function requires the seekr package to be installed
# as there are iterations of train and findhits functions, which could take long time, it is recommended to run this function on a high performance computing cluster


### Input:
# queryfadir: Path to the fasta file of query seq (e.g. functional regions of a ncRNA)
# if query fasta contains more than one sequence 
# all the sequences in query fasta file will be merged to one sequence for calculating seekr.pearson and for calculating kmer count files for hmseekr
# nullfadir: Path to the fasta file of null model sequences (e.g. transcriptome, genome, etc.)
# searchpool: Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
# bkgfadir: fasta file directory for background sequences, which serves as the normalizing factor for the input of seekr_norm_vectors and used by seekr_kmer_counts function
# this fasta file can be different from the nullfadir fasta file
# knum: a single integer value for kmer number
# queryTmin: minimal number of probability of query to query transition, this number should be greater than 0 and it is included in the iteration
# queryTmax: max number of probability of query to query transition, this number should be less than 1 and it is included in the iteration
# queryTstep: step width between queryTmin and queryTmax, numbers are limited to 6 decimal places
# nullTmin: minimal number of probability of null to null transition, this number should be greater than 0 and it is included in the iteration
# nullTmax: max number of probability of null to null transition, this number should be less than 1 and it is included in the iteration
# nullTstep: step width between nullTmin and nullTmax, numbers are limited to 6 decimal places
# lengthfilter: only keep hits sequences that have length > lengthfilter for calculating stats in the output, default=25. if no filter is needed, set to 0
# outputname: File name for output dataframe, default='gridsearch_results'
# outputdir: path of output directory to save outputs and intermediate files, default is a subfolder called gridsearch under current directory
# the intermediate fasta seq files, count files, trained models and hits files 
# are automatically saved under the outputdir into subfolders: seqs, counts, models, hits
# where qT and nT are included in the file names as the iterated transition probabilities
# alphabet: String, Alphabet to generate k-mers default='ATCG'
# progressbar: whether to show progress bar, default=True: show progress bar

### Output:
# a dataframe containing information about qT, nT, kmer number, the total number of hits sequences and the mean, median, standard deviation of the hits sequences' pearson correlation r score to the query sequence
# and the mean, median, standard deviation of the length of the hits sequences after filtering by lengthfilter

### Example:
# from hmseekr.gridsearch import gridsearch

# testsearch = gridsearch(queryfadir='../fastaFiles/mXist_rA.fa', 
#                         nullfadir='/Users/shuang/mSEEKR/fastaFiles/mm10_exp_map_200.fa', 
#                         searchpool='../fastaFiles/pool.fa',
#                         bkgfadir='/Users/shuang/mSEEKR/fastaFiles/vM25.lncRNA.can.500.nodup.fa',knum=4, 
#                         queryTmin=0.92, queryTmax=0.99, queryTstep=0.01,
#                         nullTmin=0.9, nullTmax=0.99, nullTstep=0.01,
#                         lengthfilter=25,
#                         outputname='gridsearch_results', 
#                         outputdir='/Users/shuang/gridsearch/', 
#                         alphabet='ATCG', progressbar=True)


########################################################################################################

from seekr.kmer_counts import BasicCounter as seekrBasicCounter 
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

from hmseekr import train
from hmseekr import findhits
from hmseekr import kmers

import os
import numpy as np
import pandas as pd
from tqdm import tqdm


def gridsearch(queryfadir, nullfadir, searchpool, bkgfadir, knum,
               queryTmin, queryTmax, queryTstep, nullTmin, nullTmax, nullTstep,
               lengthfilter=25, outputname='gridsearch_results',outputdir='./gridsearch/',  
               alphabet='ATCG', progressbar=True):

    # Check if specified directory exists
    # create new directory if not existing
    if not outputdir.endswith('/'):
        outputdir+='/'
    newDir = outputdir

    if not os.path.exists(newDir):
        os.mkdir(newDir)
    else:
        # check if the directory is empty
        if len(os.listdir(newDir)) != 0:
            print('The directory is not empty')
            print('Please specify a directory that is empty or does not exist for outputdir') 
            return None

    
    # based on input generate a list of queryT and nullT
    # include the min and max values in the list
    queryT_list= [round(i, 6) for i in list(np.arange(queryTmin, queryTmax+queryTstep, queryTstep))]
    nullT_list= [round(i, 6) for i in list(np.arange(nullTmin, nullTmax+nullTstep, nullTstep))]
    # filter queryT_list and nullT_list and only keep values that is between 0 and 1
    queryT_list = [i for i in queryT_list if i > 0 and i < 1]
    nullT_list = [i for i in nullT_list if i > 0 and i < 1]

    # check if queryT_list and nullT_list are empty
    if len(queryT_list) == 0 or len(nullT_list) == 0:
        print('No valid queryT or nullT values')
        print('queryT and nullT values should be between 0 and 1, but not equal to 0 or 1')
        print('Please check queryTmin, queryTmax, queryTstep, nullTmin, nullTmax, nullTstep')
        return None

    print('queryT_list:', queryT_list)
    print('nullT_list:', nullT_list)
    # load in the background sequences and calculate the seekr norm vectors
    print('Calculating background norm vectors')
    print('This could take a while if the background fasta file is large')
    bkg_norm = seekrBasicCounter(bkgfadir, k=knum, log2='Log2.post', binary=True, label=False, leave=True, silent=True, alphabet='AGTC') 
    bkg_norm.get_counts() 
    if not os.path.exists(f'{newDir}counts/'):
        os.mkdir(f'{newDir}counts/')
    mean_path = f'{newDir}counts/bkg_mean_{knum}mers.npy' 
    std_path = f'{newDir}counts/bkg_std_{knum}mers.npy' 
    np.save(mean_path, bkg_norm.mean)
    np.save(std_path, bkg_norm.std)
    print('Background norm vectors saved')



    # calculate the kmer counts of the query and null seq for hmseekr
    temp = kmers.kmers(queryfadir, str(knum), alphabet, outputdir=f'{newDir}counts/', outputname=f'hmquery{knum}')
    del temp
    print('Calculating null sequences kmer counts for hmseekr')
    print('This could take a while if the null fasta file is large')
    temp = kmers.kmers(nullfadir, str(knum), alphabet, outputdir=f'{newDir}counts/', outputname=f'hmnull{knum}')
    del temp
    print('Kmer counts for null sequences saved')
    querydir=f'{newDir}counts/hmquery{knum}.dict'
    nulldir=f'{newDir}counts/hmnull{knum}.dict'

    if not os.path.exists(f'{newDir}seqs/'):
                os.mkdir(f'{newDir}seqs/')

    # calculate the query kmer counts for seekr
    qseqs = seekrReader(queryfadir).get_seqs()
    
    if len(qseqs) > 1:
        print('More than one sequence in query fasta file')
        print('All the query sequences will be merged for calculating seekr.pearson')
        print('Merged fasta file is saved under seqs folder as seekrquery.fa')
        seekrqueryfadir = f'{newDir}seqs/seekrquery.fa'
        qseqs = '$'.join(qseqs)
        qseqfile = open(seekrqueryfadir, 'w')
        qseqfile.write('>concatenatedQuery' + '\n' + qseqs + '\n')
        qseqfile.close()

        query_count = seekrBasicCounter(infasta=seekrqueryfadir, outfile=f'{newDir}counts/seekrquery_counts.csv', k=knum, binary=False, label=True, mean=mean_path, std=std_path, log2='Log2.post', leave=True, silent=True, alphabet='ACGT') 
        query_count.make_count_file() 

    else:
        query_count = seekrBasicCounter(infasta=queryfadir, outfile=f'{newDir}counts/seekrquery_counts.csv', k=knum, binary=False, label=True, mean=mean_path, std=std_path, log2='Log2.post', leave=True, silent=True, alphabet='ACGT') 
        query_count.make_count_file() 

    # create an empty dataframe to store the results
    combstats = pd.DataFrame(columns=['qT','nT','knum','total_n','r_mean','r_median','r_std','len_mean','len_median','len_std'])

    # initiate a list to store all the pearson correlation r score
    # sim_all = pd.DataFrame(columns=['qT','nT','sim'])

    # loop through the queryT and nullT list
    total_iterations = len(queryT_list) * len(nullT_list)
    # qiter = tqdm(queryT_list,desc='queryT') if progressbar else queryT_list
    # niter = tqdm(nullT_list,desc='nullT', leave=True) if progressbar else nullT_list

    
    if progressbar:
        pbar = tqdm(total=total_iterations, desc='Overall Progress')
    else:
        pbar = None  # No progress bar

    for qT_idx, qT in enumerate(queryT_list, start=1):
        for nT_idx, nT in enumerate(nullT_list, start=1):
            if progressbar:
                # Update the progress bar description to reflect both loops' progress
                pbar.set_description(f'queryT {qT_idx}/{len(queryT_list)}, nullT {nT_idx}/{len(nullT_list)}')
                pbar.update(1)  # Update the progress bar by one step for each iteration of the inner loop

            # train the model
            if not os.path.exists(f'{newDir}models/'):
                os.mkdir(f'{newDir}models/')
            train.train(querydir, nulldir, str(knum), alphabet, qT, nT, queryPrefix=f'gsquery{qT}', nullPrefix=f'gsnull{nT}', outputdir=f'{newDir}models/')
            # ensemble the model directory
            modeldir = f'{newDir}models/gsquery{qT}_gsnull{nT}/{knum}/hmm.dict'
            hitsdir = f'{newDir}hits/'
            if not os.path.exists(f'{newDir}hits/'):
                os.mkdir(f'{newDir}hits/')
            # find the hits
            hits = findhits.findhits(searchpool=searchpool, modeldir=modeldir, knum=knum, outputname=f'hits_q{qT}_n{nT}', outputdir=hitsdir, alphabet=alphabet, fasta=True, progressbar=False)
            # only keep the rows in hits if Length col is greater than 25nt
            hits = hits[hits['Length']>lengthfilter]
            lenvec = np.array(hits['Length'])
            # save the hits sequences
            hits_seq = hits['Sequence']
            hits_header = hits['seqName']+'_'+hits['Start'].astype(str)+'_'+hits['End'].astype(str)
            hitseqdir = f'{newDir}seqs/hits_seq_q{qT}_n{nT}.fa'
            
            # save the hits sequences and its header to a fasta file
            with open(hitseqdir, 'w') as f:
                for i in range(len(hits_seq)):
                    f.write(f'{hits_header.iloc[i]}\n{hits_seq.iloc[i]}\n')

            # calculate kmer counts of hits sequences and the pearson correlation r score
            hits_count = seekrBasicCounter(infasta=hitseqdir, outfile=f'{newDir}counts/seekrhits_q{qT}_n{nT}_counts.csv', k=knum, binary=False, label=True, mean=mean_path, std=std_path, log2='Log2.post', leave=True, silent=True, alphabet='ACGT') 
            hits_count. make_count_file() 
            sim = seekrPearson(hits_count.counts,query_count.counts)
            # append the results to the dataframe
            newrow = {'qT':qT, 'nT':nT, 'knum':knum, 'total_n':len(hits_seq), 'r_mean':np.mean(sim), 'r_median':np.median(sim), 'r_std':np.std(sim), 'len_mean':np.mean(lenvec), 'len_median':np.median(lenvec), 'len_std':np.std(lenvec)}
            combstats.loc[len(combstats)] = newrow
            # append the results to sim_all
            # sim_all = sim_all.append({'qT':qT, 'nT':nT, 'sim':sim}, ignore_index=True)
        
    if progressbar:
        pbar.close()  # Close the progress bar when done

    # save the dataframe to a csv file
    combstats.to_csv(f'{newDir}{outputname}.csv', index=False)


