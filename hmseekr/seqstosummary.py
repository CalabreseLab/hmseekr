###################################################################################################
### Description: 
# This function takes in a fasta file of multiple query sequences, a transtion probabilty dataframe, a search pool fasta file, a null fasta file and a background fasta file
# and generates a summary dataframe where each row contains a sequence in the search pool fasta file
# and seven columns: seqName, feature, counts, len_sum, LLR_sum, LLR_meidan, pval_median (long format)
# or a wide format dataframe where each row corresponds to a sequence in the search pool fasta file
# and columns are eachfeature_counts, eachfeature_len_sum, eachfeature_LLR_sum, eachfeature_LLR_median, eachfeature_pval_median


### Details:
# this function is designed to get the overall likeliness of each search pool sequence to the query sequences
# the function takes in a fasta file of multiple query sequences, a transtion probabilty dataframe, a search pool fasta file, a null fasta file (for hmseekr) and a background fasta file (for seekr)
# here the transition probability dataframe must have the same rows as the query fasta file
# the columns should be '[qT,nT]' where qT is the probability of query to query transition, nT is the probability of null to null transition 
# the transition prbability for each query sequence can be different and can be optimized by the gridsearch function. 
# please include 'qT' and 'nT' as the first row (the column names) for the two columns in the csv file. 
# users can also choose to set the transition probability to be the same for all query sequences
# the function will run the kmers, train, findhits and hitseekr functions for each query sequence
# then the results can be filtered by the length of the hit regions, the kmer log likelihood ratio (kmerLLR) and the seekr pearson correlation p value
# finally for each search pool sequence, the function will calculate the counts of filtered hit regions with a specific query sequence
# and also the sum of kmerLLR and length, the median of kmerLLR and seekr pval for all the counts of a search pool sequence with each the query sequences
# for long format: each row of the output dataframe contains a sequence in the search pool fasta
# and has seven columns: seqName, feature, counts, len_sum, LLR_sum, LLR_meidan, pval_median
# seqName corresponds to the header in the search pool fasta file
# feature corresponds to the header in the query fasta file
# counts is the counts of filtered hit regions of the search pool sequences with the query sequences
# len_sum is the sum of the length of all counts of a search pool sequence with the query sequences
# LLR_sum is the sum of kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences
# LLR_median is the median of kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences
# pval_median is the median of seekr pearson correlation p value for each search pool sequence with the query sequences
# for wide format: each row of the output dataframe corresponds to a sequence in the search pool fasta
# and columns are eachfeature_counts, eachfeature_len_sum, eachfeature_LLR_sum, eachfeature_LLR_median, eachfeature_pval_median
# the output dataframe can then be used to generalize an overall likeliness of each search pool sequence to all the query sequences



### Input:
# queryfadir: Path to the fasta file of query seqs
# different from other functions such as kmers and gridsearch, if query fasta contains more than one sequence 
# each sequence will be treated as a separate query sequence
# transdf: Path to the transition probability dataframe in csv format, the dataframe should have the same rows as the query fasta file
# and the columns should be [qT,nT] where qT is the probability of query to query transition, nT is the probability of null to null transition
# please do not include the index column in the csv file
# nullfadir: Path to the fasta file of null model sequences (e.g. transcriptome, genome, etc.)
# searchpool: Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
# bkgfadir: fasta file directory for background sequences, which serves as the normalizing factor for the input of seekr_norm_vectors and used by seekr_kmer_counts function
# this fasta file can be different from the nullfadir fasta file
# knum: a single integer value for kmer number
# func: the function to use for finding hits, default='findhits_condE', other options include 'findhits'
# lenfilter: only keep hits sequences that have length > lengthfilter for calculating stats in the output, default=25. if no filter is needed, set to 0
# llrfilter: only keep hits sequences that have kmerLLR > llrfilter for calculating stats in the output, default=0. if no filter is needed, set to 0
# kmerLLR is the log likelihood ratio of of the probability 
# that the set of k-mers y within a hit derived from the QUERY versus the NULL state
# it is the sum of the log2(Q/N) ratio for each kmer within a hit
# pfilter: only keep hits sequences that have seekr pearson correlation p value < pfilter for calculating stats in the output, default=1. if no filter is needed, set to 1
# outputname: File name for output dataframe, default='seqstosummary_results'
# outputdir: path of output directory to save outputs and intermediate files, default is a subfolder called seqstosummary under current directory
# the intermediate fasta seq files, count files, trained models and hits files 
# are automatically saved under the outputdir into subfolders: seqs counts, models, hits
# outdfformat: the format of the output dataframe, default='long', other option is 'wide'
# alphabet: String, Alphabet to generate k-mers default='ATCG'
# progressbar: whether to show progress bar, default=True: show progress bar

### Output:
# a dataframe in long format: where each row contains a sequence in the search pool fasta file
# seven columns: seqName, feature, counts, len_sum, LLR_sum, LLR_meidan, pval_median
# seqName corresponds to the header in the search pool fasta file
# feature corresponds to the header in the query fasta file
# counts is the counts of filtered hit regions of the search pool sequences with the query sequences
# len_sum is the sum of the length of all counts of a search pool sequence with the query sequences
# LLR_sum is the sum of kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences
# LLR_median is the median of kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences
# pval_median is the median of seekr pearson correlation p value for each search pool sequence with the query sequences
# in wide format: each row of the output dataframe corresponds to a sequence in the search pool fasta
# and columns are eachfeature_counts, eachfeature_len_sum, eachfeature_LLR_sum, eachfeature_LLR_median, eachfeature_pval_median


### Example:
# from hmseekr.seqstosummary import seqstosummary

# testsum = seqstosummary(queryfadir='/Users/shuang/mSEEKR/fastaFiles/mXist_repeats.fa', 
#                         transdf='/Users/shuang/mSEEKR/fastaFiles/transdf.csv',
#                         nullfadir='/Users/shuang/mSEEKR/fastaFiles/mm10_exp_map_200.fa', 
#                         searchpool='/Users/shuang/mSEEKR/fastaFiles/pool.fa',
#                         bkgfadir='/Users/shuang/mSEEKR/fastaFiles/vM25.lncRNA.can.500.nodup.fa',
#                         knum=4, func='findhits_condE',
#                         lenfilter=25,llrfilter=0, pfilter=1,
#                         outputname='seqstosummary_results', 
#                         outputdir='/Users/shuang/seqstosummary/', 
#                         outdfformat='long',alphabet='ATCG', 
#                         progressbar=True)


########################################################################################################

from seekr.kmer_counts import BasicCounter as seekrBasicCounter 
from seekr.fasta_reader import Reader as seekrReader
from seekr.find_pval import find_pval
from seekr.find_dist import find_dist

from hmseekr.train import train
from hmseekr.kmers import kmers


import os
import numpy as np
import pandas as pd
from tqdm import tqdm


def seqstosummary(queryfadir, transdf, nullfadir, searchpool, bkgfadir, knum,
                  func='findhits_condE', lenfilter=25, llrfilter=0, pfilter=1,
                  outputname='seqstosummary_results',outputdir='./seqstosummary/', 
                  outdfformat='long',alphabet='ATCG', progressbar=True):
    
    # read in the query fasta file
    queryseqs = seekrReader(queryfadir).get_seqs()
    queryheaders = seekrReader(queryfadir).get_headers()
    # remove the '>' from the headers
    queryheaders = [x[1:] for x in queryheaders]

    # check if transdf is a csv file
    if not transdf.endswith('.csv'):
        print('The transition probability dataframe should be in csv format')
        print('Please check the input file')
        return None

    # read in the transition probability dataframe
    tdf = pd.read_csv(transdf)
    # check if tdf has two columns
    if len(tdf.columns) != 2:
        print('The transition probability dataframe should have two columns:qT and nT')
        print('Please check the input file')
        return None

    # test if the transition probability dataframe has the same rows as the query fasta file
    if len(queryseqs) != len(tdf): 
        print('The number of query sequences in the fasta file does not match the number of rows in the transition probability dataframe')
        print('Please check the input files')
        return None

    # rename the column names of tdf
    tdf.columns = ['qT','nT']
    print('make sure the transition probability dataframe 1st column is qT and 2nd column is nT')
    qTvec=tdf['qT']
    nTvec=tdf['nT']

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

    # load in corresponding findhits function
    if func == 'findhits':
        from hmseekr import findhits
        from hmseekr.findhits import findhits as findhits_cur
    elif func == 'findhits_condE':
        from hmseekr import findhits_condE
        from hmseekr.findhits_condE import findhits_condE as findhits_cur
    else:
        print('Please specify a valid function for finding hits')
        print('Options include: findhits, findhits_condE,')
        return None
    

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

    if not os.path.exists(f'{newDir}models/'):
        os.mkdir(f'{newDir}models/')

    # find the dist of background
    print('Fitting background model')
    fitres = find_dist(inputseq=bkgfadir, k_mer=knum, models='common10', 
                       subsetting=True, subset_size = 10000, 
                       fit_model=True, statsmethod='ks',progress_bar=progressbar, 
                       plotfit=f'{outputdir}models/modelfitplot', outputname=f'{outputdir}models/fitres')
    # delete the intermediate files generated by find_dist
    # if bkg_mean_{knum}mers.npy and bkg_std_{knum}mers.npy exist
    # remove them
    if os.path.exists(f'bkg_mean_{knum}mers.npy'):
        os.remove(f'bkg_mean_{knum}mers.npy')
    if os.path.exists(f'bkg_std_{knum}mers.npy'):
        os.remove(f'bkg_std_{knum}mers.npy')
        
    # calculate the kmer counts of the null seq for hmseekr
    print('Calculating null sequences kmer counts for hmseekr')
    print('This could take a while if the null fasta file is large')
    kmers(nullfadir, str(knum), alphabet, outputname=f'hmnull_{knum}',outputdir=f'{newDir}counts/')
    print('Kmer counts for null sequences saved')
    nulldir=f'{newDir}counts/hmnull_{knum}.dict'

    # initialize the output dataframe
    combstats = pd.DataFrame()


    # loop thru each query sequence

    iterable = tqdm(range(len(queryseqs))) if progressbar else range(len(queryseqs))

    for i in iterable:

        # create subdirectories for each query sequence
        if not os.path.exists(f'{newDir}seqs/'):
            os.mkdir(f'{newDir}seqs/')

        # write the query sequence to a fasta file
        queryfa=f'{newDir}seqs/query{i}.fa'
        qseqfile = open(queryfa, 'w')
        qseqfile.write(f'>{queryheaders[i]}\n{queryseqs[i]}\n')
        qseqfile.close()

        # calculate the kmer counts for the query sequence
        kmers(queryfa, str(knum), alphabet, outputname=f'hmquery{i}_{knum}',outputdir=f'{newDir}counts/')
        querydir=f'{newDir}counts/hmquery{i}_{knum}.dict'

        # train the model
        train(querydir, nulldir, str(knum), alphabet, qTvec[i], nTvec[i], queryPrefix=f'query{i}_q{qTvec[i]}', nullPrefix=f'n{nTvec[i]}', outputdir=f'{newDir}models/')
        # ensemble the model directory
        modeldir = f'{newDir}models/query{i}_q{qTvec[i]}_n{nTvec[i]}/{knum}/hmm.dict'
        hitsfd = f'{newDir}hits/'

        if not os.path.exists(f'{newDir}hits/'):
            os.mkdir(f'{newDir}hits/')

        # find hits
        print('Finding hits for query sequence',i+1)
        hits = findhits_cur(searchpool=searchpool, modeldir=modeldir, knum=knum, outputname=f'query{i}hits_q{qTvec[i]}_n{nTvec[i]}', outputdir=hitsfd, alphabet=alphabet, fasta=True, progressbar=progressbar)

        # hitseekr, do not use the function as we can recycle the bkg norm vecs and model fits

        # only keep the rows in hits if Length col is greater than lengthfilter
        hits = hits[hits['Length']>lenfilter]

        # check if there are still hits after filtering
        if len(hits) == 0:
            print('for query sequence',i+1,'in the query fasta file')
            print('No hits after length filtering')
            print('Please try a lower lenfilter value')
            # skip the rest of the loop and continue to the next query sequence
            print('continue to the next query sequence')
            continue
        
        else:
            # Filter based on normalized LLR 
            #hits['normLLR'] = hits['kmerLLR'] / hits['Length']
            hits = hits[hits['kmerLLR'] > llrfilter]

            # check if there are still hits after filtering
            if len(hits) == 0:
                print('for query sequence',i+1,'in the query fasta file')
                print('No hits after kmerLLR filtering')
                print('Please try a lower llrfilter value')
                # skip the rest of the loop and continue to the next query sequence
                print('continue to the next query sequence')
                continue
            else: 
                # proceed to calculate the seekr pearson correlation p value
                # save the hits sequences
                hits_seq = hits['Sequence']
                hits_header = hits['seqName']+'_'+hits['Start'].astype(str)+'_'+hits['End'].astype(str)
                hitseqdir = f'{newDir}seqs/query{i}hits_seqs.fa'
                
                # save the hits sequences and its header to a fasta file
                with open(hitseqdir, 'w') as f:
                    for n in range(len(hits_seq)):
                        f.write(f'{hits_header.iloc[n]}\n{hits_seq.iloc[n]}\n')
                
                # find the p value of the seekr r score
                print('Finding seekr p values for query sequence',i+1)
                pvals=find_pval(seq1file=hitseqdir, seq2file=queryfa, 
                                mean_path=mean_path, std_path=std_path,
                                k_mer=knum, fitres=fitres, log2='Log2.post', 
                                bestfit=1, outputname=f'{newDir}hits/query{i}hits_q{qTvec[i]}_n{nTvec[i]}_seekr_pval', progress_bar=progressbar)
                
                # change the first column name of pvals
                pvals.columns = ['seekr_pval']
                pvals = pvals.reset_index(drop=True)
                hits = hits.reset_index(drop=True)

                # add the pvals seekr_pval values to the hits dataframe as a new column seekr_pval
                hits['seekr_pval'] = pvals['seekr_pval']

                # filter with pvals
                hits = hits[hits['seekr_pval'] < pfilter]

                if len(hits) == 0:
                    print('for query sequence',i+1,'in the query fasta file')
                    print('No hits after seekr p val filtering')
                    print('Please try a lower pfilter value')
                    # skip the rest of the loop and continue to the next query sequence
                    print('continue to the next query sequence')
                    continue
                else: 
                    # save the hits dataframe
                    hits.to_csv(f'{newDir}hits/query{i}hits_q{qTvec[i]}_n{nTvec[i]}_filtered.csv',index=False)

                    # add the feature name to the hits dataframe
                    hits['feature'] = queryheaders[i]

                    # remove the '>' from the seqName
                    hits['seqName'] = hits['seqName'].str[1:]

                    # add the hits dataframe to the combstats dataframe
                    combstats = pd.concat([combstats,hits])

    # check if combstats is empty
    if len(combstats) == 0:
        print('No hits for any query sequence after all filtering')
        print('Please try different filtering values')
        return None
    else: 
        # group by the search pool sequence and the feature name
        combstats_summary = combstats.groupby(['seqName', 'feature']).agg({
            'Start': 'count',
            'Length': 'sum',
            'kmerLLR': ['sum', 'median'],  # Aggregating 'kmerLLR' with both sum and median
            'seekr_pval': 'median'   
        }).reset_index()

        # Flatten the MultiIndex columns and rename them
        combstats_summary.columns = ['seqName', 'feature', 'counts', 'len_sum', 'LLR_sum', 'LLR_median', 'pval_median']

        if outdfformat == 'long':
            # save the summary dataframe
            combstats_summary.to_csv(f'{newDir}{outputname}_long.csv',index=False)

            return combstats_summary
        elif outdfformat == 'wide':
            # pivot the dataframe
            combstats_summary_wide = combstats_summary.pivot(index='seqName',columns='feature',values=['counts', 'len_sum', 'LLR_sum', 'LLR_median', 'pval_median'])
            # Flatten the MultiIndex columns
            combstats_summary_wide.columns = [f'{col[1]}_{col[0]}' for col in combstats_summary_wide.columns]
            # Reset index to make 'seqName' a column again
            combstats_summary_wide = combstats_summary_wide.reset_index()
            # save the summary dataframe
            combstats_summary_wide.to_csv(f'{newDir}{outputname}_wide.csv',index=False)

            return combstats_summary_wide








        
         



    


