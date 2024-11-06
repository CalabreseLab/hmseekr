# hmseekr


This is a program for identifying regions of high similarity based on *k*-mer content to some set of query sequences, relative to a null background set of sequences.

<hr/>

## Dependencies

#### Python3.9
check with command 

```
python -V
```

if < 3.9, please install python 3.9 or above

<hr/>

## Installation

#### Install through Python Package Index (PyPI)
```
pip install hmseekr
```
This will make both the command line tool and the python module available.

#### Install through Github
```
pip install git+https://github.com/CalabreseLab/hmseekr.git
```
This will make both the command line tool and the python module available.

#### Install through Docker Hub
First you need to install Docker on your local computer. Then pull the Docker Image:
```
docker pull calabreselab/hmseekr:latest
```
This will install the Docker container which enables running hmseekr from the command line or Jupyter Notebook. See below the hmseekr Docker Image section for more details.

<hr/>

## Pipeline example

### kmers: counting kmer profile for query and background sequence 

This step counts k-mers for multiple specified values of k and saves them to a binary file.
The hmseekr_kmers function takes in fasta file such as the query sequence or the background sequences. If input fasta contains more than one sequence, all the sequences will be merged to one sequence before processing. It counts the kmer frequency for each kmer size in --kvec (-k), --kvec can also be a single kmer size. The output of the function is a dictionary of dictionaries with the outer dictionary keys as kmer size and the inner dictionary keys as kmer sequences and values as kmer counts.

This step should be performed for both query sequence and the background sequences, which then generate the .dict kmer profile for both query and background sequences.

#### Console Example:
generate kmer count files for kmer size 2, 3, 4 using mXist_rA.fa as input fasta file
```
hmseekr_kmers -fd './fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -name repeatA -dir './counts/' 
```

#### Python Example:
generate kmer count files for kmer size 4 using mXist_rA.fa as input fasta file
```python
from hmseekr.kmers import kmers

testdict = kmers(fadir='./fastaFiles/mXist_rA.fa',kvec='4',
                 alphabet='ATCG',outputname='repeatA',
                 outputdir='./counts/')
```

#### Inputs:

1. fadir (-fd) : Path to input fasta file, should be the query sequence or the background sequences.
2. kvec (-k) : Comma delimited string of possible k-mer values. For example, 3,4,5 or just 4.
3. alphabet (-a) : Alphabet to generate k-mers (e.g. ATCG). Default is 'ATCG'.
4. outputname (-name) : Desired output name for count file. Default is 'out'.
5. outputdir (-dir) : Directory to save output count file. Default is './', that is current directory.


#### Output:

Outputs binary .dict files containing count matrices


### train: training markov models

This step calculates the emission matrix based on kmer counts, prepare transition matrix, states and saves them to a binary file. This function takes in kmer count file for sequences of interest (e.g. query seq, functional regions of a lncRNA) and kmer count file for background sequences (e.g. null seq, transcriptome, genome). These kmer count files are generated using the kmers function with the k specified (which should be the same as calculted in the kmers function) and query to query transition rate (qT) and null to null transition rate (nT) specified.
It calculates the hidden state transition matrix, hidden state emission matrix, states and starting probability of each hidden state and saves them to a binary file.

#### Console Example:
train a model using previously generated kmer count files for repeatA and all lncRNA (kmers, or hmseekr_kmers function) with kmer size 4 and transition rates of 0.9999 for both query to query and null to null. Save the model to the current directory.
```
hmseekr_train -qd './counts/repeatA.dict' -nd './counts/all_lncRNA.dict' -k 4 -a ATCG -qT 0.9999 -nT 0.9999 -qPre repeatA -nPre lncRNA -dir './'
```

#### Python Example:
train a model using previously generated kmer count files for repeatA and all lncRNA (kmers, or hmseekr_kmers function) with kmer size 2,3,4 and transition rates of 0.9999 for both query to query and null to null, and save to the markovModels folder
```python
from hmseekr.train import train

testmodel = train(querydir='./counts/repeatA.dict', nulldir='./counts/all_lncRNA.dict',
                  kvec='2,3,4', alphabet='ATCG', queryT=0.9999, nullT=0.9999,
                  queryPrefix='repeatA', nullPrefix='lncRNA', outputdir='./markovModels/')
```

#### Inputs:

1. querydir (-qd): Path to kmer count file for sequence of interest or query seq (e.g. functional regions of a ncRNA)
2. nulldir (-nd): Path to kmer count file that compose null model/background sequences (e.g. transcriptome, genome, etc.)
3. kvec (-k): Comma delimited string of possible k-mer values. For example, '3,4,5' or just '4'. Numbers in kvec must be found in the k-mer count file (precalculated by kmers function)
4. alphabet (-a): Alphabet to generate k-mers, default=ATCG
5. queryT (-qT): Probability of query to query transition, default=0.9999, should be between 0 and 1 but not equal to 0 or 1
6. nullT (-nT): Probability of null to null transition, default=0.9999, should be between 0 and 1 but not equal to 0 or 1
7. queryPrefix (-qPre): prefix file name for query, defualt='query'
8. nullPrefix (-nPre): prefix file name for null, defualt='null'
9. outputdir (-dir): path of output directory to save output trained model file in .dict format, default is current directory


#### Output:

  Directory containing models (.dict file) for each value of k specified, the directory structure would look like:

    | markovModels
    |
    |--- repeatA_lncRNA
    |------- 2
    |------------- hmm.dict
    |------- 3
    |------------- hmm.dict
    |------- 4
    |------------- hmm.dict
    |--- repeatB_lncRNA
    .
    .
    .


### findhits: find high similar regions based on kmer profile within sequences of interest

This step uses precalculated model (emission matrix, prepared transition matrix, pi and states) from train (hmseekr_train) function to find out HMM state path through sequences of interest, therefore return sequences that has high similarity (based on kmer profile) to the query sequence -- hits sequneces. This function takes in a fasta file (searchpool) which defines the region to search for potential hits (highly similar regions to the query sequence), also takes in the precalculated model (train or hmseekr_train function). Along the searchpool fasta sequences, similarity scores to query sequence will be calculated based on the model, hits segments (highly similar regions) would be reported along with the sequence header from the input fasta file, start and end location of the hit segment, kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment. kmerLLR is defined as the sum of the log likelihood of each k-mer in the hit sequence being in the Q (query) state minus the log likelihood of them being in the N (null) state

#### Console Example:
use the previously trained model (hmm.dict by train/hmseekr_train function) to search for highly similar regions to query sequence (repeatA) within the pool.fa files (area of interest region to find sequences similar to repeatA, could be all lncRNAs or just chromosome 6) with kmer size 4 and save the hit sequences while showing progress bar
```
hmseekr_findhits -pool './fastaFiles/pool.fa' -m './markovModels/repeatA_lncRNA/4/hmm.dict' -k 4 -name 'hits' -dir './models/' -a 'ATCG' -fa -pb
```

#### Python Example:
use the previously trained model (hmm.dict by train/hmseekr_train function) to search for highly similar regions to query sequence (repeatA) within the pool.fa files (area of interest region to find sequences similar to repeatA, could be all lncRNAs or just chromosome 6) with kmer size 4 and save the hit sequences while showing progress bar
```python
from hmseekr.findhits import findhits

testhits = findhits(searchpool='./fastaFiles/pool.fa',
                    modeldir='./markovModels/repeatA_lncRNA/4/hmm.dict',
                    knum=4,outputname='hits',outputdir='./',
                    alphabet='ATCG',fasta=True,
                    progressbar=True)
```

#### Inputs:

1. searchpool (-pool): Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
2. modeldir (-m): Path to precalculated model .dict file output from train.py'
3. knum (-k): Value of k to use as an integer. Must be the same as the k value used in training (train function) that produced the model
4. outputname (-name): File name for output, useful to include information about the experiment, default='hits'
5. outputdir (-dir): Directory to save output dataframe. Default is './', that is current directory.
6. alphabet (-a): String, Alphabet to generate k-mers default='ATCG'
7. fasta (-fa): whether to save sequence of hit in the output dataframe, default=True: save the actual sequences
8. progressbar (-pb): whether to show progress bar

#### Output:
A dataframe containing information about the hit regions: highly similar regions to query sequence based on the precalculated model within the input fasta file. Information about the hits regions includes: the sequence header from the input fasta file, start and end location of the hit segment, kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment if fasta=True.


### findhits_condE: findhits with conditioned Emission probabilities 

this is a variant of findhits/hmseekr_findhits function: it calculates the emission probability of the next word given the current word
in findhits, as the shift is by 1 nt, the next word has k-1 overlap with the current word
for example, if the current word is 'TAGC', the next possible words are 'AGCA', 'AGCT', 'AGCC', 'AGCG'
then the emission probability of the next word ('AGCA') given the current word as 'TAGC' is calculated as 
np.log2(emission probability of 'AGCA' in the original E) - np.log2(sum of emission probability of all four possible worlds in the original E)


#### Console Example:
same example as hmseekr_findhits but with the conditioned emission probability
```
hmseekr_findhits_condE -pool './fastaFiles/pool.fa' -m './markovModels/repeatA_lncRNA/4/hmm.dict' -k 4 -name 'hits' -dir './models/' -a 'ATCG' -fa -pb
```

#### Python Example:
same example as findhits but with the conditioned emission probability
```python
from hmseekr.findhits_condE import findhits_condE

testhits = findhits_condE(searchpool='../fastaFiles/pool.fa',
                          modeldir='../markovModels/hmm.dict',
                          knum=4,outputname='hits',outputdir='./',
                          alphabet='ATCG',fasta=True,
                          progressbar=True)
```

#### Inputs and Output:

Inputs and Output are in the same format as findhits/hmseekr_findhits function


### findhits_nol: findhits with non-overlapping kmers

this is a variant of findhits/hmseekr_findhits function: it uses non-overlapping kmers and shifts the sequence by 1nt each time until the kmer size is reached. in this way, it reduces the kmer dependency of overlapping kmers.
if sequence is 'ATGCTTTTGCGC' the kmers would be 'ATGC','TTTT','GCGC'
then it chops off 1nt from the start of each sequence in the searchpool and re-scan with non-overlapping kmers
so the sequence would be 'TGCTTTTGCGC' and the kmers would be 'TGCT','TTTG'
the last bit of sequences that is less than the kmer size would be discarded
this is done until the kmer size is reached
in this way, it generates k result txt files each starts at different position of the sequence with non-overlapping kmers
this function reduces the kmer dependency of overlapping kmers in findhits, but also considered the different combinations of kmers at different start positions
the outcome txt files can be further process to find the best hit regions
for example, only use the regions that are hits in all k result txt files

#### Console Example:
same example as hmseekr_findhits but with non-overlapping kmers 
```
hmseekr_findhits_nol -pool './fastaFiles/pool.fa' -m './markovModels/repeatA_lncRNA/4/hmm.dict' -k 4 -name 'hits' -dir './models/' -a 'ATCG' -fa -pb
```

#### Python Example:
same example as findhits but with non-overlapping kmers 
```python
from hmseekr.findhits_nol import findhits_nol

testhits = findhits_nol(searchpool='../fastaFiles/pool.fa',
                        modeldir='../markovModels/hmm.dict',
                        knum=4,outputname='hits',outputdir='./',
                        alphabet='ATCG',fasta=True,
                        progressbar=True)
```

#### Inputs:

Inputs are in the same format as findhits/hmseekr_findhits function

#### Output:

Output has the same format as findhits/hmseekr_findhits function. But instead of one output txt file, it will generate k output txt files, each corresponds to a shift position. 



### gridsearch: search for best transition probabilities

This function performs a grid search to find the best trasnition probabilities for query to query and null to null states, which is used in the train function. This function takes in query sequence, null squence and background sequence (see Inputs for details) fasta files, ranges and steps for query to query transition rate (qT) and null to null transition rate (nT), a specific kmer number and performs the train function and findhits function for each combination of qT and nT. Within the hits sequences (findhits function results), only keep the sequence with length greater than 25nt. Then calculates pearson correlation r score (seekr.pearson) between the filtered hit sequences (findhits results) and the query sequence. It returns a dataframe (.csv file) containing the qT, nT, kmer number, the total number of hits sequences and the mean, median, standard deviation of the hits sequences' pearson correlation r score to the query sequence and the mean, median, standard deviation of the length of the hits sequences. If query fasta contains more than one sequence, all the sequences in query fasta file will be merged to one sequence, for calculating kmer count files for hmseekr (hmseekr.kmers) and for calculating seekr.pearson. This function requires the seekr package to be installed. As there are iterations of train and findhits functions, which could take long time, it is recommended to run this function on a high performance computing cluster. Variants of findhits functions can be specified to run.


#### Console Example:
perform a grid search to find the best transition probabilities for qT and nT each within the range of 0.9 to 0.99 with step of 0.01 with lengthfilter set to 25, which only keep the hit sequences with length greater than 25 for stats calculation. the regular 'findhits' function is used here.

```
hmseekr_gridsearch -qf './fastaFiles/repeatA.fa' -nf './fastaFiles/all_lncRNA.fa' -pool './fastaFiles/pool.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -qTmin 0.9 -qTmax 0.99 -qTstep 0.01 -nTmin 0.9 -nTmax 0.99 -nTstep 0.01 -fc 'findhits' -lf 25 -name 'gridsearch_results' -dir './gridsearch/' -a 'ATCG' -pb
```

#### Python Example:
perform a grid search to find the best transition probabilities for qT and nT each within the range of 0.9 to 0.99 with step of 0.01 with lengthfilter set to 25, which only keep the hit sequences with length greater than 25 for stats calculation. the regular 'findhits' function is used here.

```python
from hmseekr.gridsearch import gridsearch

testsearch = gridsearch(queryfadir='./fastaFiles/repeatA.fa', 
                        nullfadir='./fastaFiles/all_lncRNA.fa', 
                        searchpool='./fastaFiles/pool.fa',
                        bkgfadir='./fastaFiles/bkg.fa',knum=4, 
                        queryTmin=0.9, queryTmax=0.99, queryTstep=0.01,
                        nullTmin=0.9, nullTmax=0.99, nullTstep=0.01,
                        func='findhits_condE', lengthfilter=25,
                        outputname='gridsearch_results', 
                        outputdir='./gridsearch/', 
                        alphabet='ATCG', progressbar=True)
```

#### Inputs:

1. queryfadir (-qf): Path to the fasta file of query seq (e.g. functional regions of a ncRNA). If query fasta contains more than one sequence, all the sequences in query fasta file will be merged to one sequence for calculating seekr.pearson and for calculating kmer count files for hmseekr
2. nullfadir (-nf): Path to the fasta file of null model sequences (e.g. transcriptome, genome, etc.)
3. searchpool (-pool): Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
4. bkgfadir (-bkgf): fasta file directory for background sequences, which serves as the normalizing factor for the input of seekr_norm_vectors and used by seekr_kmer_counts function. This fasta file can be different from the nullfadir fasta file
5. knum (-k): a single integer value for kmer number
6. queryTmin (-qTmin): minimal number of probability of query to query transition, this number should be greater than 0 and it is included in the iteration
7. queryTmax (-qTmax): max number of probability of query to query transition, this number should be less than 1 and it is included in the iteration
8. queryTstep (-qTstep): step width between queryTmin and queryTmax, numbers are limited to 6 decimal places
9. nullTmin (-nTmin): minimal number of probability of null to null transition, this number should be greater than 0 and it is included in the iteration
10. nullTmax (-nTmax): max number of probability of null to null transition, this number should be less than 1 and it is included in the iteration
11. nullTstep (-nTstep): step width between nullTmin and nullTmax, numbers are limited to 6 decimal places
12. func (-fc): the function to use for finding hits, default='findhits_condE', other options include 'findhits'
13. lengthfilter (-lf): only keep hits sequences that have length > lengthfilter for calculating stats in the output, default=25. if no filter is needed, set to 0
14. outputname (-name): File name for output dataframe, default='gridsearch_results'
15. outputdir (-dir): path of output directory to save outputs and intermediate files, default is a subfolder called gridsearch under current directory. The intermediate fasta seq files, count files, trained models and hits files are automatically saved under the outputdir into subfolders: seqs, counts, models, hits, where qT and nT are included in the file names as the iterated transition probabilities
16. alphabet (-a): String, Alphabet to generate k-mers default='ATCG'
17. progressbar (-pb): whether to show progress bar, default=True: show progress bar

#### Output:
a dataframe containing information about qT, nT, kmer number, the total number of hits sequences and the mean, median, standard deviation of the hits sequences' pearson correlation r score to the query sequence and the mean, median, standard deviation of the length of the hits sequences after filtering by lengthfilter


### fastarev: reverse fasta sequences

This function reverses a fasta file and save it to a new file:'ATGC' will be reversed to 'CGTA'. This function takes in fasta file and reverse the sequence. If input fasta contains more than one sequence, each sequence will be reversed. Headers will be kept the same, only the sequences will be reversed. Keeping the headers the same will allow the user to match the reversed sequences to the original sequences easily. 

#### Console Example:
reverse the sequence of mouse Xist repeat A fasta file and save it to a new file while keeping the headers the same
```
hmseekr_fastarev -i '../fastaFiles/mXist_rA.fa' -o '../fastaFiles/mXist_rA_rev.fa'
```

#### Python Example:
reverse the sequence of mouse Xist repeat A fasta file and save it to a new file while keeping the headers the same
```python
from hmseekr.fastarev import fastarev

fastarev(input_file_path='../fastaFiles/mXist_rA.fa',
         output_file_path='../fastaFiles/mXist_rA_rev.fa')
```

#### Inputs:

1. input_file_path (-i): path to the input fasta file.
2. output_file_path (-o): path and name to the output fasta file.


#### Output:
A fasta file with the reversed sequences and the SAME headers as the input fasta file. If the input fasta file contains more than one sequence, each sequence will be reversed.


### genbed: filter hits sequences and generate bedfiles for regular/forward fasta

This function firstly filter the hits output from findhits/hmseekr_findhits and then generate a bed file for the filtered hits. This function only applies to the output from findhits/hmseekr_findhits with the regular fasta file as input. For reversed fasta file, please use genbedrev/hmseekr_genbedrev. This function takes in the output from hmseekr_findhits and filter the hits based on the hit length and normalized kmerLLR score. kmerLLR is the log likelihood ratio of of the probability that the set of k-mers y within a hit derived from the QUERY versus the NULL state. It is the sum of the log2(Q/N) ratio for each kmer within a hit. The normalized kmerLLR (normLLR) is the kmerLLR divided by the hit length. So the normLLR is the log2 of the Q/N ratio, i.e. if set normLLR > 0.5, the Q/N ratio is ~ 1.41. Then all the hits after the filtering will be converted to a bedfile

#### Console Example:
generate bedfile for filtered hits from the hmseekr_findhits output file mm10expmap_queryA_4_viterbi.txt
```
hmseekr_genbed -hd '../mm10expmap_queryA_4_viterbi.txt' -o '../mm10expmap_queryA_4_viterbi' -len 25 -llr 0.5 -pb
```

#### Python Example:
generate bedfile for filtered hits from the findhits output file mm10expmap_queryA_4_viterbi.txt
```python
from hmseekr.genbed import genbed

testbed = genbed(hitsdir='../mm10expmap_queryA_4_viterbi.txt', 
                 outputdir='../mm10expmap_queryA_4_viterbi',
                 lenfilter=25, llrfilter=0.5,progressbar=True)
```

#### Inputs:

1. hitsdir (-hd): path to the input hits file, should be the output from findhits with the regular fasta file as input. please use genbedrev/hmseekr_genbedrev for the reversed fasta file 
2. outputdir (-o): path and name to the output bedfile, do not need to add .bed at the end
3. lenfilter (-len): the minimum length of the hit, only keep hits that has a length > lenfilter, default is 25
4. llrfilter (-llr): the minimum normalized kmerLLR score, only keep hits that has a normLLR > llrfilter, default is 0.5
5. progressbar (-pb): whether to show the progress bar, default is True


#### Output:
A bedfile with the following columns: 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand' for all hits that passed the filtering. 'score' is the normLLR score, 'name' is 'fwd' as it should only be applied for the regular fasta file. For reversed fasta file, please use genbedrev/hmseekr_genbedrev


### genbedrev: filter hits sequences and generate bedfiles for reversed fasta

This function firstly filter the hits output from findhits/hmseekr_findhits with the reversed fasta file as input and then generate a bed file for the filtered hits. This function only applies to the output from findhits/hmseekr_findhits with the reversed fasta file as input for regular or forward fasta file, please use genbed/hmseekr_genbed. This function takes in the output from findhits/hmseekr_findhits with the reversed fasta file as input and filter the hits based on the hit length and normalized kmerLLR score. kmerLLR is the log likelihood ratio of of the probability that the set of k-mers y within a hit derived from the QUERY versus the NULL state. It is the sum of the log2(Q/N) ratio for each kmer within a hit. The normalized kmerLLR (normLLR) is the kmerLLR divided by the hit length. So the normLLR is the log2 of the Q/N ratio, i.e. if set normLLR > 0.5, the Q/N ratio is ~ 1.41. Then all the hits after the filtering will be converted to a bedfile.

#### Console Example:
generate bedfile for filtered hits from the hmseekr_findhits output file FLIPmm10expmap_queryA_4_viterbi.txt, which is generated using the reversed fasta file as input
```
hmseekr_genbedrev -hd '../FLIPmm10expmap_queryA_4_viterbi.txt' -o '../FLIPmm10expmap_queryA_4_viterbi' -len 25 -llr 0.5 -pb
```

#### Python Example:
generate bedfile for filtered hits from the findhits output file FLIPmm10expmap_queryA_4_viterbi.txt, which is generated using the reversed fasta file as input
```python
from hmseekr.genbedrev import genbedrev

testbed = genbedrev(hitsdir='../FLIPmm10expmap_queryA_4_viterbi.txt', 
                    outputdir='../FLIPmm10expmap_queryA_4_viterbi',
                    lenfilter=25, llrfilter=0.5,progressbar=True)
```

#### Inputs:

1. hitsdir (-hd): path to the input hits file, should be the output from findhits with the reversed fasta file as input. please use genbed/hmseekr_genbed for the regular/forward fasta file 
2. outputdir (-o): path and name to the output bedfile, do not need to add .bed at the end
3. lenfilter (-len): the minimum length of the hit, only keep hits that has a length > lenfilter, default is 25
4. llrfilter (-llr): the minimum normalized kmerLLR score, only keep hits that has a normLLR > llrfilter, default is 0.5
5. progressbar (-pb): whether to show the progress bar, default is True


#### Output:
A bedfile with the following columns: 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand' for all hits that passed the filtering. 'score' is the normLLR score, 'name' is 'rev' as it should only be applied for the reversed fasta file. For regular/forward fasta file, please use genbed/hmseekr_genbed.


### hitseekr: add seekr p values

This function takes in the output of findhits function and the query sequence fasta file with a background fasta file, and calculates the seekr pearson correlation p value between the hit sequences and the query sequence. It fit the background sequences to the common10 distributions and takes the best ranked distribution. It calculates the seekr pearson correlation p values between the hit sequences and the query sequence based on the best ranked distribution. It adds the seekr p value to the hits dataframe on top of the existing kmer log likelihood score (kmerLLR). The seekr p value could provide additional information about the similarity between the hit sequences and the query sequence.

#### Console Example:
add seekr p value to the hits dataframe generated by findhits function

```
hmseekr_hitseekr -hd './mm10_queryA_4_viterbi.txt' -qf './fastaFiles/repeatA.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -lf 25 -dir './' -pb
```

#### Python Example:
add seekr p value to the hits dataframe generated by findhits function

```python
from hmseekr.hitseekr import hitseekr

addpvals = hitseekr(hitsdir='./mm10_queryA_4_viterbi.txt',
                    queryfadir='./fastaFiles/mXist_rA.fa', 
                    bkgfadir='./fastaFiles/all_lncRNA.fa',
                    knum=4, lengthfilter=25, outputdir='./', progressbar=True)
```

#### Inputs:

1. hitsdir (-hd): Path to the directory of the .txt file generated by findhits function
2. queryfadir (-qf): Path to the fasta file of query seq (e.g. functional regions of a ncRNA). If query fasta contains more than one sequence, all the sequences in query fasta file will be merged to one sequence for calculating seekr.pearson and for calculating kmer count files for hmseekr
3. bkgfadir (-bkgf): fasta file directory for background sequences, which serves as the normalizing factor for the input of seekr_norm_vectors and used by seekr_kmer_counts function. This fasta file can be different from the nullfadir fasta file
4. knum (-k): a single integer value for kmer number, must be the same as the kmer number used in the findhits function
5. lengthfilter (-lf): only keep hits sequences that have length > lengthfilter for calculating stats in the output, default=25. if no filter is needed, set to 0
6. outputdir (-dir): path of output directory, default is current directory, save the final dataframe together with other intermediate files
7. progressbar (-pb): whether to show progress bar, default=True: show progress bar

#### Output:
a dataframe with seekr p value added to the findhits dataframe
outputname is automatically generated as the input findhits filename with '\_seekr' appended to it


### seqstosummary: search and quantify multiple features

This function is designed to get the overall likeliness of each search pool sequence to the query sequences. The function takes in a fasta file of multiple query sequences, a transtion probabilty dataframe, a search pool fasta file, a null fasta file (for hmseekr) and a background fasta file (for seekr). Here the transition probability dataframe must have the same rows as the query fasta file. The columns should be '\[qT,nT\]' where qT is the probability of query to query transition, nT is the probability of null to null transition The transition prbability for each query sequence can be different and can be optimized by the gridsearch function. Please include 'qT' and 'nT' as the first row (the column names) for the two columns in the csv file. Please check the transdf.csv file in the repository for an example. Users can also choose to set the transition probability to be the same for all query sequences. The function will run the kmers, train, findhits and hitseekr functions for each query sequence. Then the results can be filtered by the length of the hit regions, the normalized kmer log likelihood ratio (kmerLLR) and the seekr pearson correlation p value. Finally for each search pool sequence, the function will calculate the counts of filtered hit regions with a specific query sequence, and also the sum of length normalized kmerLLR (normLLR) and length for all the counts of a search pool sequence with each the query sequences. For long format each row of the output dataframe contains a sequence in the search pool fasta, and has five columns: seqName, feature, counts, sum_normLLR, sum_len: seqName corresponds to the header in the search pool fasta file; feature corresponds to the header in the query fasta file; counts is the counts of filtered hit regions of the search pool sequences with the query sequences; sum_normLLR is the sum of length normalized kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences; sum_len is the sum of the length of all counts of a search pool sequence with the query sequences. For wide format: each row of the output dataframe corresponds to a sequence in the search pool fasta and columns are eachfeature_counts, eachfeature_sum_normLLR, eachfeature_sum_len. The output dataframe can then be used to generalize an overall likeliness of each search pool sequence to all the query sequences



#### Console Example:
search all genes on chr16 for the potential hit counts and similarities to the query sequences include mXist repeat A, B, C and E

```
hmseekr_seqstosummary -qf './fastaFiles/mXist_repeats.fa' -td './transdf.csv' -nf './fastaFiles/all_lncRNA.fa' -pool './fastaFiles/chr16.fa' -bkgf './fastaFiles/bkg.fa' -k 4 -fc 'findhits_condE' -lenf 25 -llrf 0 -pf 1 -name 'seqstosummary_results' -dir './seqstosummary/' -format long -a 'ATCG' -pb
```

#### Python Example:
search all genes on chr16 for the potential hit counts and similarities to the query sequences include mXist repeat A, B, C and E

```python
from hmseekr.seqstosummary import seqstosummary

testsum = seqstosummary(queryfadir='./fastaFiles/mXist_repeats.fa', 
                        transdf='./trans.csv',
                        nullfadir='./fastaFiles/expressed_lncRNA.fa', 
                        searchpool='./fastaFiles/chr16.fa',
                        bkgfadir='./all_lncRNA.fa',
                        knum=4, func='findhits_condE',
                        lenfilter=25,llrfilter=0, pfilter=1,
                        outputname='seqstosummary_results', 
                        outputdir='/Users/shuang/seqstosummary/', 
                        outdfformat='long', alphabet='ATCG', progressbar=True)
```

#### Inputs:

1. queryfadir (-qf): Path to the fasta file of query seqs. Different from other functions such as kmers and gridsearch, if query fasta contains more than one sequence, each sequence will be treated as a separate query sequence
2. transdf (-td): Path to the transition probability dataframe in csv format, the dataframe should have the same rows as the query fasta file, and the columns should be '\[qT,nT\]' where qT is the probability of query to query transition, nT is the probability of null to null transition. Please include 'qT' and 'nT' as the first row (the column names) for the two columns in the csv file. Please do not include the index column in the csv file.
3. nullfadir (-nf): Path to the fasta file of null model sequences (e.g. transcriptome, genome, etc.)
4. searchpool (-pool): Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
5. bkgfadir (-bkgf): fasta file directory for background sequences, which serves as the normalizing factor for the input of seekr_norm_vectors and used by seekr_kmer_counts function. This fasta file can be different from the nullfadir fasta file
6. knum (-k): a single integer value for kmer number
7. func (-fc): the function to use for finding hits, default='findhits_condE', other options include 'findhits'
8. lenfilter (-lenf): only keep hits sequences that have length > lenfilter for calculating stats in the output, default=25. if no filter is needed, set to 0
9. llrfilter (-llrf): only keep hits sequences that have length normalized kmerLLR > llrfilter for calculating stats in the output, default=0. if no filter is needed, set to 0. kmerLLR is the log likelihood ratio of of the probability that the set of k-mers y within a hit derived from the QUERY versus the NULL state. It is the sum of the log2(Q/N) ratio for each kmer within a hit. The normalized kmerLLR (normLLR) is the kmerLLR divided by the hit length. So the normLLR is the log2 of the Q/N ratio, i.e. if set normLLR > 0.5, the Q/N ratio is ~ 1.41
10. pfilter (-pf): only keep hits sequences that have seekr pearson correlation p value < pfilter for calculating stats in the output, default=1. if no filter is needed, set to 1
11. outputname (-name): File name for output dataframe, default='seqstosummary_results'
12. outputdir (-dir): path of output directory to save outputs and intermediate files, default is a subfolder called seqstosummary under current directory. The intermediate fasta seq files, count files, trained models and hits files are automatically saved under the outputdir into subfolders: seqs, counts, models, hits, where qT and nT are included in the file names
13. outdfformat (-format): the format of the output dataframe, default='long', other option is 'wide'
14. alphabet (-a): String, Alphabet to generate k-mers default='ATCG'
15. progressbar (-pb): whether to show progress bar, default=True: show progress bar

#### Output:
a dataframe in long format: where each row contains a sequence in the search pool fasta file, and five columns: seqName, feature, counts, sum_normLLR, sum_len: seqName corresponds to the header in the search pool fasta file; feature corresponds to the header in the query fasta file; counts is the counts of filtered hit regions of the search pool sequences with the query sequences; sum_normLLR is the sum of length normalized kmer log likelihood ratio (kmerLLR) for each search pool sequence with the query sequences; sum_len is the sum of the length of all counts of a search pool sequence with the query sequences. In wide format: each row of the output dataframe corresponds to a sequence in the search pool fasta, and columns are eachfeature_counts, eachfeature_sum_normLLR, eachfeature_sum_len




<hr/>

## Example pipeline to run on both forward and reverse sequences to improve hits quality

Here we show an example of searching for similar regions of mouse Xist repeat A across mouse mm10 chromosome 16. Firstly we build and run the model for searching forward throughout chr16 and yield a bedfile after filtering for hit length and normalized LLR scores. Then we reverse the mouse Xist repeat A and the mm10 chromosome 16. We build and run the model on the reversed sequences and yield a bedfile in a similar way. In this way, we have the model searching for high similarity regions in both directions. If only one direction is run, all the hits are identified based on its previous sequences, which is how Hidden Markov Model works. When running in both directions, we can have hits identified based on preceding and succeeding sequences. By either merging or intersecting the bedfiles from both directions, we can improve the quality of the hits. All findhits/hmseekr_findhits function can be substituted by its variant: findhits_condE/hmseekr_findhits_condE or findhits_nol/hmseekr_findhits_nol as needed.

### run forwardly

query sequence: mXist_rA.fa
background sequence: mm10_all_lncRNA.fa
searchpool sequence: mm10_chr16.fa

#### Console Example:

count kmers for query and background sequences
```
hmseekr_kmers -fd '../fastaFiles/mXist_rA.fa' -k 4 -a ATCG -name repeatA_4 -dir '../counts/' 

hmseekr_kmers -fd '../fastaFiles/mm10_all_lncRNA.fa' -k 4 -a ATCG -name lncRNA_4 -dir '../counts/' 
```

train the model
```
hmseekr_train -qd '../counts/repeatA_4.dict' -nd '../counts/lncRNA_4.dict' -k 4 -a ATCG -qT 0.9999 -nT 0.9999 -qPre repeatA -nPre lncRNA -dir '../markovModels/'
```

find potential hits across chr16
```
hmseekr_findhits -pool '../fastaFiles/mm10_chr16.fa' -m '../markovModels/repeatA_lncRNA/4/hmm.dict' -k 4 -name 'chr16_queryA' -dir '../' -a 'ATCG' -fa -pb
```

generate bedfiles after filtering: as this is the forward run, we use hmseekr_genbed

```
hmseekr_genbed -hd '../chr16_queryA_4_viterbi.txt' -o '../chr16_queryA_4_viterbi' -len 25 -llr 0.5 -pb
```

#### Python Example:

```python
from hmseekr.kmers import kmers
from hmseekr.train import train
from hmseekr.findhits import findhits
from hmseekr.genbed import genbed

# count kmers for query and background sequences
qdict = kmers(fadir='../fastaFiles/mXist_rA.fa',kvec='4',
              alphabet='ATCG',outputname='repeatA_4',
              outputdir='../counts/')

bkgdict = kmers(fadir='../fastaFiles/mm10_all_lncRNA.fa',kvec='4',
                alphabet='ATCG',outputname='lncRNA_4',
                outputdir='../counts/')

# train the model
fwdmodel = train(querydir='../counts/repeatA_4.dict', nulldir='../counts/lncRNA_4.dict',
                 kvec='4', alphabet='ATCG', queryT=0.9999, nullT=0.9999,
                 queryPrefix='repeatA', nullPrefix='lncRNA', outputdir='../markovModels/')

# find potential hits across chr16
fwdhits = findhits(searchpool='../fastaFiles/mm10_chr16.fa',
                   modeldir='../markovModels/repeatA_lncRNA/4/hmm.dict',
                   knum=4,outputname='chr16_queryA',outputdir='../',
                   alphabet='ATCG',fasta=True,progressbar=True)

# generate bedfiles after filtering: as this is the forward run, we use genbed
fwdbed = genbed(hitsdir='../chr16_queryA_4_viterbi.txt', 
                outputdir='../chr16_queryA_4_viterbi',
                lenfilter=25, llrfilter=0.5,progressbar=True)

```


### run reversely

query sequence: reversed mXist_rA.fa
background sequence: reversed mm10_all_lncRNA.fa
searchpool sequence: reversed mm10_chr16.fa

#### Console Example:

reverse query, background and searchpool sequences. Please use the hmseekr_fastarev to reverse sequences to ensure proper formats for later steps
```
hmseekr_fastarev -i '../fastaFiles/mXist_rA.fa' -o '../fastaFiles/mXist_rA_rev.fa'

hmseekr_fastarev -i '../fastaFiles/mm10_all_lncRNA.fa' -o '../fastaFiles/mm10_all_lncRNA_rev.fa'

hmseekr_fastarev -i '../fastaFiles/mm10_chr16.fa' -o '../fastaFiles/mm10_chr16_rev.fa'
```

count kmers for reversed query and background sequences
```
hmseekr_kmers -fd '../fastaFiles/mXist_rA_rev.fa' -k 4 -a ATCG -name repeatA_4_rev -dir '../counts/' 

hmseekr_kmers -fd '../fastaFiles/mm10_all_lncRNA_rev.fa' -k 4 -a ATCG -name lncRNA_4_rev -dir '../counts/' 
```

train the model
```
hmseekr_train -qd '../counts/repeatA_4_rev.dict' -nd '../counts/lncRNA_4_rev.dict' -k 4 -a ATCG -qT 0.9999 -nT 0.9999 -qPre repeatA -nPre lncRNA_rev -dir '../markovModels/'
```

find potential hits reversely across chr16 
```
hmseekr_findhits -pool '../fastaFiles/mm10_chr16_rev.fa' -m '../markovModels/repeatA_lncRNA_rev/4/hmm.dict' -k 4 -name 'Rev_chr16_queryA' -dir '../' -a 'ATCG' -fa -pb
```

generate bedfiles after filtering: as this is the reverse run, we use hmseekr_genbedrev

```
hmseekr_genbedrev -hd '../Rev_chr16_queryA_4_viterbi.txt' -o '../Rev_chr16_queryA_4_viterbi' -len 25 -llr 0.5 -pb
```

#### Python Example:

```python
from hmseekr.fastarev import fastarev
from hmseekr.kmers import kmers
from hmseekr.train import train
from hmseekr.findhits import findhits
from hmseekr.genbedrev import genbedrev

# reverse query, background and searchpool sequences
# Please use the hmseekr_fastarev to reverse sequences to ensure proper formats for later steps
fastarev(input_file_path='../fastaFiles/mXist_rA.fa',
         output_file_path='../fastaFiles/mXist_rA_rev.fa')

fastarev(input_file_path='../fastaFiles/mm10_all_lncRNA.fa',
         output_file_path='../fastaFiles/mm10_all_lncRNA_rev.fa')

fastarev(input_file_path='../fastaFiles/mm10_chr16.fa',
         output_file_path='../fastaFiles/mm10_chr16_rev.fa')

# count kmers for reversed query and background sequences
revqdict = kmers(fadir='../fastaFiles/mXist_rA_rev.fa',kvec='4',
                 alphabet='ATCG',outputname='repeatA_4_rev',
                 outputdir='../counts/')

revbkgdict = kmers(fadir='../fastaFiles/mm10_all_lncRNA_rev.fa',kvec='4',
                   alphabet='ATCG',outputname='lncRNA_4_rev',
                   outputdir='../counts/')

# train the model
revmodel = train(querydir='../counts/repeatA_4_rev.dict', nulldir='../counts/lncRNA_4_rev.dict',
                 kvec='4', alphabet='ATCG', queryT=0.9999, nullT=0.9999,
                 queryPrefix='repeatA', nullPrefix='lncRNA_rev', outputdir='../markovModels/')

# find potential hits reversely across chr16
revhits = findhits(searchpool='../fastaFiles/mm10_chr16_rev.fa',
                   modeldir='../markovModels/repeatA_lncRNA_rev/4/hmm.dict',
                   knum=4,outputname='Rev_chr16_queryA',outputdir='../',
                   alphabet='ATCG',fasta=True,progressbar=True)

# generate bedfiles after filtering: as this is the reverse run, we use genbedrev
revbed = genbedrev(hitsdir='../Rev_chr16_queryA_4_viterbi.txt', 
                   outputdir='../Rev_chr16_queryA_4_viterbi',
                   lenfilter=25, llrfilter=0.5,progressbar=True)

```


### further analysis based on the bedfiles from both run

Each bedfiles would have the score column (5th), which is the normalized LLR score. We can treat this score as the quaility of each hit: the higher the score, the more similar the hit to the query. Therefore, if needed, hits region inside the bedfiles can also be sorted based on the score column. Users then can decided to apply further filters such as only taking the top 100 hits. This step can be easily done in Excel, R or Python. Please choose the method that works best for you.


For both forward and reverse run, the bedfiles have the correct start and end cooridnates. That is, although the reverse run searches the chr16 reversely, the hits region has been properly recognized and transformed to a bedfile with correct start and end position. These coordinates are not reversed. Therefore, we can use [bedtools](https://bedtools.readthedocs.io/en/latest/index.html) in a “command line” environment to either [merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) the hits or [intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) the hits. Please refer to the [bedtools](https://bedtools.readthedocs.io/en/latest/index.html) official website for installation guide.


#### merge the bedfiles

To merge the bedfiles, we firstly need to concatenate both bedfiles into one file
```
cat chr16_queryA_4_viterbi.bed Rev_chr16_queryA_4_viterbi.bed > Comb_chr16_queryA_4_viterbi.bed
```

Before merging, it is always good to sort the bedfiles
```
sort -k1,1 -k2,2n Comb_chr16_queryA_4_viterbi.bed > Comb_chr16_queryA_4_viterbi_sorted.bed
```

Merge the bedfiles. here we decide to take the mean score (5th column, normalized LLR) of the overlapping regions as the new score for the merged region. If this operation needs to be changed, please refer to the [merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) function website.

```
# load the module if working on HPC
module load bedtools

bedtools merge -i Comb_chr16_queryA_4_viterbi_sorted.bed -s -c 4,5,6 -o distinct,mean,distinct > chr16_queryA_4_viterbi_sorted_merged.bed
```
Users can choose to further filter the hits in the chr16_queryA_4_viterbi_sorted_merged.bed file based on the score column (5th), by either enforcing a threshold or order and choose the top 100 hits.


#### intersect the bedfiles

Before intersecting, it is always good to sort the bedfiles
```
sort -k1,1 -k2,2n chr16_queryA_4_viterbi.bed > chr16_queryA_4_viterbi_sorted.bed

sort -k1,1 -k2,2n Rev_chr16_queryA_4_viterbi.bed > Rev_chr16_queryA_4_viterbi_sorted.bed
```

Intersect the bedfiles. Here we use -s to enforce the intersection on the same strand. The name and score for the intersected region will be inherited from the -a input. For example, -a has the entry: chr1 0 10 fwd 0.9 + and -b has entry: chr1 5 12 rev 0.2 +. The intersected region will be reported as: chr1 5 10 fwd 0.9 +. If this operation needs to be changed, please refer to the [intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) function website.

```
# load the module if working on HPC
module load bedtools

bedtools intersect -a chr16_queryA_4_viterbi_sorted.bed -b Rev_chr16_queryA_4_viterbi_sorted.bed -s > chr16_queryA_4_viterbi_intsct.bed
```
The -wb argument can be added to report the original -b feature in the same row with the intersected feature. In this way, the first 6 columns are for the intersected feature, with the 4th (name), 5th (score) same as the -a feature, and the following 6 columns report the original -b feature. Users can recalcualte the score (5th column) for the intersected region by considering scores from both -a and -b feature. 
```
# load the module if working on HPC
module load bedtools

bedtools intersect -a chr16_queryA_4_viterbi_sorted.bed -b Rev_chr16_queryA_4_viterbi_sorted.bed -s -wb > chr16_queryA_4_viterbi_intsct.bed
```
Users can choose to further filter the hits in the chr16_queryA_4_viterbi_sorted_intsct.bed file based on the score column (5th), by either enforcing a threshold or order and choose the top 100 hits.


#### get fasta sequence from the bedfiles

Using bedtools [getfasta](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html) we can also easily extract the fasta sequence for each regions defined in the bedfile
```
# load the module if working on HPC
module load bedtools

bedtools getfasta -fi mm10_genome.fa -bed chr16_queryA_4_viterbi_intsct.bed -fo chr16_queryA_4.fa -s
```




<hr/>

## hmseekr Docker Image
We now provide a Docker image of hmseekr which hopefully avoids all the package dependencies and complications when install and run hmseekr through pip. 

### Docker Installation
Firstly, you should install Docker on your computer. Please refer to the [official website](https://www.docker.com/get-started/) for installation instructions.

After you have successfully installed Docker. Start/Run the application and make sure it is running properly – you should see the docker icon on your task bar with the status indicated.

### Pull Docker Image and Test Run
1. Start your command line tool: Terminal for MacOS and CMD for Windows. You can also use Powershell or Cygwin for Windows, but Cygwin might have interaction issues.

From the command line, pull the Docker Image:
```
docker pull calabreselab/hmseekr:latest
```
You can replace `latest` with a specific tag if needed.

2. Test Run the Docker Image
```
docker run -it --rm calabreselab/hmseekr:latest
```
The `-it` tag sets it to interactive mode. If you don't need to run the Docker container in interactive mode (i.e., you don't need a shell inside the container), you can simply omit the `-it` flag.
This will print the user manual out to the command line, which is basically the same as you run the command `hmseekr` directly from command line when you pip install seekr. 

### Run Docker Image from command line
You can run the seekr function from this Docker Image directly from command line with the following specified syntax.
```
docker run -v /path/to/your/files:/data calabreselab/hmseekr:latest hmseekr_kmers -fd '/data/fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -name repeatA -dir '/data/counts/'
```
In this command:
* `-v /path/to/your/files:/data`: This mounts the directory `/path/to/your/files` from your host machine (where /fastaFiles/mXist_rA.fa is located) to the `/data` directory inside the Docker container. Replace `/path/to/your/files` with the actual path to your files.
* `mseekr_kmers -fd '/data/fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -name repeatA -dir '/data/counts/'`: This is the command that gets executed inside the Docker container. Since we mounted our files to `/data` in the container, we reference them with `/data/fastaFiles/mXist_rA.fa`. 
* The `/data` folder is basically a mirror of the folder you specified in `/path/to/your/files`. So by specifying `-dir '/data/counts/'` (output into the counts folder under the path /data/) we can have the output files directly written in to the folder in `/path/to/your/files`.
* Please remember to **specify your output path to `/data/`** otherwise it will not be saved to your folder on local machine and it would be hard to locate it even inside the Docker Container Filesystem (in this case, when the Docker Container is removed, your data will be deleted as well). 

Examples of code mounts e:/test on Windows as the folder that contains the input and holds the output files:
```
docker run -v e:/test:/data calabreselab/hmseekr:latest hmseekr_kmers -fd '/data/fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -name repeatA -dir '/data/counts/'
```
Basically you need to add: `docker run -v /path/to/your/files:/data calabreselab/hmseekr:latest` before the command line code for hmseekr (see above for examples of all functions).

### Run Docker Image with Jupyter Notebook
If you want to work with python codes instead of directly calling from the command line, you can choose to run seekr with Jupyter Notebook inside the Docker Container.

1.  Run Docker Container with Interactive Terminal:
```
docker run -it -p 8888:8888 -v /path/on/host:/path/in/container calabreselab/hmseekr:latest /bin/bash
```
This command will start the container and give you a bash terminal inside it. The `-p 8888:8888` flag maps the port *8888* from the container to your host machine so that you can access the Jupyter Notebook server.

`/path/on/host` is the path to a directory on your local machine where you want the notebooks and other data to be saved. `/path/in/container` is the directory inside the container where you want to access and save the files from Jupyter Notebook.

When you use Jupyter Notebook now and create or save notebooks, they will be stored at `/path/in/container` inside the Docker container. Due to the volume mount, these files will also be accessible and stored at `/path/on/host` on your host machine so that you can later access the code and files even when the container is removed. 

Example of code:
```
docker run -it -p 8888:8888 -v e:/test:/data calabreselab/hmseekr:latest /bin/bash
```

2. Manually start Jupyter Notebook. From the bash terminal inside the Docker container:
```
jupyter notebook --ip=0.0.0.0 --port=8888 --NotebookApp.token='' --NotebookApp.password='' --allow-root
```
* `--ip=0.0.0.0`: Allow connections from any IP address. This is important for accessing Jupyter from outside the container.
* `--port=8888`: Run Jupyter on this port. You can change this if needed, but remember to adjust the `-p` flag when starting the Docker container.
* `--allow-root`: Allow Jupyter to be run as the root user inside the Docker container.
* `--NotebookApp.token=''`: This disables the token-based security measure.
* `--NotebookApp.password=''`: This ensures that there's no password required to access the Jupyter server. 

Disabling the token and password for Jupyter Notebook reduces security. It's generally okay for local development, but you should avoid doing this in production or any publicly accessible server.

3. Access the Jupyter Notebook from your host machine's browser by entering this address:
<http://localhost:8888/>

4. Run Python functions as demonstrated above. It would be convenient if all input files can be copied over to the folder you have mounted: `/path/on/host`. Pay attention that when you specify your input and output file route, use the `/path/in/container` route, as that is a mirror to your local folder. 

Example of code:
```python
from hmseekr.kmers import kmers

testdict = kmers(fadir='./fastaFiles/mXist_rA.fa',kvec='4',
                 alphabet='ATCG',outputname='repeatA',
                 outputdir='./counts/')
```
Once you are done, you can click the **Shut Down** button under the **File** tab of Jupyter Notebook to shut down the instance or you can just click `Ctrl+C` twice from the command line to kill the process. 
Then you need to exit the Docker Container from the interative session by typing `exit` and hit enter.

### Cleanup (Optional)
* Clean Up Docker Container:
    + List all containers, including the stopped ones: `docker ps -a`
    + To remove a specific container: `docker rm CONTAINER_ID_OR_NAME`
    + To remove all stopped containers: `docker container prune`

* Clean Up Docker Image. If you want to remove the Docker image you used:
    + List all images: `docker images`
    + Remove a specific image: `docker rmi IMAGE_ID_OR_NAME`

You'd typically only remove an image if you're sure you won't be using it again soon, or if you want to fetch a fresh version of it from a repository like Docker Hub.

* Additional Cleanup. Docker also maintains a cache of intermediate images and volumes. Over time, these can accumulate. To free up space:
    + Remove unused data: `docker system prune`
    + To also remove unused volumes (be careful, as this might remove data you want to keep): `docker system prune --volumes`

Remember to always be cautious when cleaning up, especially with commands that remove data. Make sure you have backups of any essential data, and always double-check what you're deleting.

<hr/>

## Issues and questions

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues 
Contact the authors at: shuang9@email.unc.edu; mauro_calabrese@med.unc.edu

