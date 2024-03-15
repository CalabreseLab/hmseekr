# hmseekr

This is a program for identifying regions of high similarity based on *k*-mer content to some set of query sequences, relative to a null background set of sequences.

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
The hmseekr_kmers function takes in fasta file such as the query sequence or the background sequences. It counts the kmer frequency for each kmer size in --kvec (-k), --kvec can also be a single kmer size. The output of the function is a dictionary of dictionaries with the outer dictionary keys as kmer size and the inner dictionary keys as kmer sequences and values as kmer counts.

This step should be performed for both query sequence and the background sequences, which then generate the .dict kmer profile for both query and background sequences.

#### Console Example:
generate kmer count files for kmer size 2, 3, 4 using mXist_rA.fa as input fasta file
```
  hmseekr_kmers -fd './fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -dir './counts/' -name repeatA
```

#### Python Example:
generate kmer count files for kmer size 4 using mXist_rA.fa as input fasta file
```python
from hmseekr import kmers

testdict = kmers(fadir='./fastaFiles/mXist_rA.fa',kvec='4',
                 alphabet='ATCG',outputdir='./counts/',
                 outputname='repeatA')
```

#### Inputs:

1. fadir (-fd) : Path to input fasta file, should be the query sequence or the background sequences.
2. kvec (-k) : Comma delimited string of possible k-mer values. For example, 3,4,5 or just 4.
3. alphabet (-a) : Alphabet to generate k-mers (e.g. ATCG). Default is 'ATCG'.
4. outputdir (-dir) : Directory to save output count file. Default is './', that is current directory.
5. outputname (-name) : Desired output name for count file. Default is 'out'.

#### Output:

Outputs binary .dict files containing count matrices

<hr/>

## train: training markov models

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
from hmseekr import train

testmodel = train(querydir='./counts/repeatA.dict', nulldir='./counts/all_lncRNA.dict',
                  kvec='2,3,4', alphabet='ATCG', queryT=0.9999, nullT=0.9999,
                  queryPrefix='repeatA', nullPrefix='lncNRA', outputdir='./markovModels/')
```

#### Inputs:

1. querydir (-qd): Path to kmer count file for sequence of interest or query seq (e.g. functional regions of a ncRNA)
2. nulldir (-nd): Path to kmer count file that compose null model/background sequences (e.g. transcriptome, genome, etc.)
3. kvec (-k): Comma delimited string of possible k-mer values. For example, '3,4,5' or just '4'. Numbers in kvec must be found in the k-mer count file (precalculated by kmers function)
4. alphabet (-a): Alphabet to generate k-mers, default=ATCG
5. queryT (-qT): Probability of query to query transition, default=0.9999
6. nullT (-nT): Probability of null to null transition, default=0.9999
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

<hr/>

## findhits: find high similar regions based on kmer profile within sequences of interest

This step uses precalculated model (emission matrix, prepared transition matrix, pi and states) from train (hmseekr_train) function to find out HMM state path through sequences of interest, therefore return sequences that has high similarity (based on kmer profile) to the query sequence -- hits sequneces. This function takes in a fasta file (searchpool) which defines the region to search for potential hits (highly similar regions to the query sequence), also takes in the precalculated model (train or hmseekr_train function). Along the searchpool fasta sequences, similarity scores to query sequence will be calculated based on the model, hits segments (highly similar regions) would be reported along with the sequence header from the input fasta file, start and end location of the hit segment, kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment.

#### Console Example:
use the previously trained model (hmm.dict by train/hmseekr_train function) to search for highly similar regions to query sequence (repeatA) within the pool.fa files (area of interest region to find sequences similar to repeatA, could be all lncRNAs or just chromosome 6) with kmer size 4 and save the hit sequences while showing progress bar
```
  hmseekr_findhits -pool './fastaFiles/pool.fa' -m './markovModels/repeatA_lncRNA/4/hmm.dict' -k 4 -name 'hits' -a 'ATCG' -fa -pb
```

#### Python Example:
use the previously trained model (hmm.dict by train/hmseekr_train function) to search for highly similar regions to query sequence (repeatA) within the pool.fa files (area of interest region to find sequences similar to repeatA, could be all lncRNAs or just chromosome 6) with kmer size 4 and save the hit sequences while showing progress bar
```python
from hmseekr import findhits

testhits = findhits(searchpool='./fastaFiles/pool.fa',
                    modeldir='./markovModels/repeatA_lncRNA/4/hmm.dict',
                    knum=4,outputname='hits',
                    alphabet='ATCG',fasta=True,
                    progressbar=True)
```

#### Inputs:

1. searchpool (-pool): Path to fasta file which defines the region to search for potential hits (highly similar regions) based on the precalculated model (train function)
2. modeldir (-m): Path to precalculated model .dict file output from train.py'
3. knum (-k): Value of k to use as an integer. Must be the same as the k value used in training (train function) that produced the model
4. outputname (-name): File name for output, useful to include information about the experiment, default='hits'
5. alphabet (-a): String, Alphabet to generate k-mers default='ATCG'
6. fasta (-fa): whether to save sequence of hit in the output dataframe, default=True: save the actual sequences
7. progressbar (-pb): whether to show progress bar

#### Output:
A dataframe containing information about the hit regions: highly similar regions to query sequence based on the precalculated model within the input fasta file. Information about the hits regions includes: the sequence header from the input fasta file, start and end location of the hit segment, kmer log likelihood score (kmerLLR), and the actual sequence of the hit segment if fasta=True.



## Issues and questions

For more details of the inputs and outputs, please refer to the manual listed under https://github.com/CalabreseLab/hmseekr/
Any issues can be reported to https://github.com/CalabreseLab/hmseekr/issues 
Contact the authors at: shuang9@email.unc.edu; mauro_calabrese@med.unc.edu

