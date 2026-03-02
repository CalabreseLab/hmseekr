###################################
# core functions from seekr package
# when updating these functions in seekr please also update these here

import numpy as np
from itertools import product
import pandas as pd
from tqdm import tnrange, trange, tqdm, tqdm_notebook

from collections import defaultdict
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import kstest

import warnings
import os
import sys


"""
Created on Mon Jul  4 13:10:49 2016

@author: jessime
"""

# based on seekr.fasta_reader function Reader as seekrReader

class seekrReader:
    """Fixes any compatibility issues a fasta file might have with this code.

    Parameters
    ----------
    infasta : str (default=None)
        Name of input fasta file to be manipulated
    outfasta : str (default=None)
        location to store extracted data from infasta
    names : iter (default=None)
        Common style names to use in header lines

    Attributes
    ----------
    data : list
        Raw lines of the infasta file
        Note: This is different than the data attribute in other classes

    Examples
    --------
    Putting the sequence on one line instead of breaking it every 80 chars.
    Making sure the whole sequence is capitalized.
    Restructuring the name line to work with GENCODE's naming.
    """

    def __init__(self, infasta=None, outfasta=None, names=None):
        self.infasta = infasta
        self.outfasta = outfasta
        self.names = names

        self.data = None

    def _read_data(self):
        """Sets data to stripped lines from the fasta file
        """
        with open(self.infasta) as infasta:
            self.data = [l.strip() for l in infasta]

    def _upper_seq_per_line(self):
        """Sets data to upper case, single line sequences for each header
        """
        new_data = []
        seq = ""
        for i, line in enumerate(self.data):
            if line[0] == ">":
                if seq:
                    new_data.append(seq.upper())
                    seq = ""
                else:
                    assert i == 0, "There may be a header without a sequence at line {}.".format(i)
                new_data.append(line)
            else:
                seq += line
        new_data.append(seq.upper())
        self.data = new_data

    def get_lines(self):
        self._read_data()
        self._upper_seq_per_line()
        return self.data

    def get_seqs(self):
        clean_data = self.get_lines()
        seqs = clean_data[1::2]
        return seqs

    def get_headers(self):
        clean_data = self.get_lines()
        headers = clean_data[::2]
        return headers

    def get_data(self, tuples_only=False):
        clean_data = self.get_lines()
        headers = clean_data[::2]
        seqs = clean_data[1::2]
        tuples = zip(headers, seqs)
        if tuples_only:
            return tuples
        else:
            return tuples, headers, seqs

    def supply_basic_header(self):
        """Convert headerlines to GENCODE format with only common name and length"""
        new_fasta = []

        if self.names is None:
            self.names = iter(self.get_headers())
        for i, line in enumerate(self.data):
            if line[0] == ">":
                name = next(self.names).strip(">")
                length = len(self.data[i + 1])
                new_fasta.append(">||||{}||{}|".format(name, length))
            else:
                new_fasta.append(line)
        return new_fasta

    def save(self):
        """Write self.data to a new fasta file"""
        with open(self.outfasta, "w") as outfasta:
            for line in self.data:
                outfasta.write(line + "\n")


"""
Created on Mon Sep 12 11:16:50 2016

@author: jessime

A couple of functions for autodetecting if code is being run in the notebook.
The appropriate tqdm progressbar will be returned.
"""

# from seekr.my_tqdm import my_tqdm

def _is_kernel():
    if "IPython" not in sys.modules:
        # IPython hasn't been imported, definitely not
        return False
    from IPython import get_ipython

    # check for `kernel` attribute on the IPython instance
    return getattr(get_ipython(), "kernel", None) is not None


def my_tqdm():
    return tqdm_notebook if _is_kernel() else tqdm


def my_trange():
    return tnrange if _is_kernel() else trange


"""
Description
-----------
Generates a kmer count matrix of m rows by n columns,
where m is the number of transcripts in a fasta file and n is 4^kmer.

Examples
--------
Generate a small, plain text .csv file:

```
from seekr.kmer_counts import BasicCounter

counter = BasicCounter(infasta='example.fa',
                       outfile='example_2mers.csv',
                       k=2,
                       binary=False)
counts = counter.make_count_file()
```

Notes
-----


Issues
------
Any issues can be reported to https://github.com/CalabreseLab

---
"""

class seekrBasicCounter:
    """Generates overlapping kmer counts for a fasta file

    Parameters
    ----------
    infasta: str (default=None)
        Full path to fasta file to be counted
    outfile: str (default=None)
        Full path to the counts file to be saved
    k: int (default=6)
        Size of kmer to be counted
    binary: bool (default=True)
        Saves as numpy array if True, else saves as csv
    mean: bool, np.array, str (default=True)
        Set the mean to 0 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated mean array.
    std: bool, np.array, str (default=True)
        Set the std. dev. to 1 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated std array.
    log2: Log2 (default='Log2.post')
        Log2 transformation can occur pre- or post-standardization, or not at all.
    leave: bool (default=True)
        Set to False if get_counts is used within another tqdm loop
    silent: bool (default=False)
        Set to True to turn off tqdm progress bar
    alphabet: str (default='AGTC')
        Valid letters to include in kmer.

    Attributes
    ----------
    counts: np.array
        Matrix of kmer counts. Dimensions are equal to # of transcripts by # of kmers.
    kmers: list
        Str elements of all kmers of size k
    map: dict
        Mapping of kmers to column values
    alpha_len: int
        Length of alphabet
    """

    def __init__(
        self,
        infasta=None,
        outfile=None,
        k=6,
        binary=True,
        mean=True,
        std=True,
        log2='Log2.post',
        leave=True,
        silent=False,
        label=False,
        alphabet="AGTC",
    ):
        self.infasta = infasta
        self.seqs = None
        if infasta is not None:
            self.seqs = seekrReader(infasta).get_seqs()
        self.outfile = outfile
        self.k = k
        self.binary = binary
        self.mean = mean
        if isinstance(mean, str):
            self.mean = np.load(mean)
        self.std = std
        if isinstance(std, str):
            self.std = np.load(std)
        self.log2 = log2
        self.leave = leave
        self.silent = silent
        self.label = label
        self.counts = None
        self.alpha_len = len(alphabet)
        self.kmers = ["".join(i) for i in product(alphabet, repeat=k)]
        self.map = {k: i for k, i in zip(self.kmers, range(self.alpha_len ** k))}

        if self.seqs is not None:
            if len(self.seqs) == 1 and self.std is True:
                err = (
                    "You cannot standardize a single sequence. "
                    "Please pass the path to an std. dev. array, "
                    "or use raw counts by setting std=False."
                )
                raise ValueError(err)

        # check if self.log2 is one of ['Log2.pre', 'Log2.post', 'Log2.none']
        if self.log2 not in ['Log2.pre', 'Log2.post', 'Log2.none']:
            raise ValueError("log2 must be one of ['Log2.pre', 'Log2.post', 'Log2.none']")

        # if not isinstance(self.log2, Log2):
        #     raise TypeError(f"log2 must be one of {list(Log2)}")

    def occurrences(self, row, seq):
        """Counts kmers on a per kilobase scale"""
        counts = defaultdict(int)
        length = len(seq)
        increment = 1000 / (length - self.k + 1)  # fix, original 1000/length, incorrect divisor, should be 1000/(l-k+1)
        for c in range(length - self.k + 1):
            kmer = seq[c : c + self.k]
            counts[kmer] += increment
        for kmer, n in counts.items():
            if kmer in self.map:
                row[self.map[kmer]] = n
        return row

    def _progress(self):
        """Determine which iterator to loop over for counting"""
        if self.silent:
            return self.seqs

        if not self.leave:
            tqdm_seqs = my_tqdm()(self.seqs, desc="Kmers", leave=False)
        else:
            tqdm_seqs = my_tqdm()(self.seqs)

        return tqdm_seqs

    def center(self):
        """Mean center counts by column"""
        if self.mean is True:
            self.mean = np.mean(self.counts, axis=0)
        self.counts -= self.mean

    def standardize(self):
        """Divide out the standard deviations from columns of the count matrix"""
        if self.std is True:
            self.std = np.std(self.counts, axis=0)
        self.counts /= self.std
        if np.isnan(self.counts).any():
            print(
                (
                    "\nWARNING: You have `np.nan` values in your counts "
                    "after standardization. This is likely due to "
                    "a kmer not appearing in any of your sequences. "
                    "Try: \n1) using a smaller kmer size, \n2) beginning "
                    "with a larger set of sequences, \n3) passing "
                    "precomputed normalization vectors from a larger "
                    "data set (e.g. GENCODE)."
                )
            )

    def log2_norm(self):
        """Apply a log2 transform to the count matrix"""
        self.counts += 1
        self.counts = np.log2(self.counts)

    def get_counts(self):
        """Generates kmer counts for a fasta file"""
        self.counts = np.zeros([len(self.seqs), self.alpha_len ** self.k], dtype=np.float32)
        seqs = self._progress()

        for i, seq in enumerate(seqs):
            self.counts[i] = self.occurrences(self.counts[i], seq)
        if self.log2 == 'Log2.pre':
            self.log2_norm()
        if self.mean is not False:
            self.center()
        if self.std is not False:
            self.standardize()
        if self.log2 == 'Log2.post':
            self.counts += np.abs(np.min(self.counts))
            self.log2_norm()

    def save(self, names=None):
        """Saves the counts appropriately based on current settings.

        There are four output methods for the counts:
        1. Binary. This saves just the counts as a binary numpy array.
        2. No labels. Saves in plain text, but without any labels.
        3. Default names. If no names are provided, fasta headers will be used as labels.
        4. Custom names. Provide a list of names if you want to label lncRNAs with your own names.

        Parameters
        ----------
        names : [str] (default=None)
            Unique names for rows of the Dataframe.
        """
        err_msg = (
            "You cannot label a binary file. "
            'Set only one of "binary" or "label" as True. '
            "If you used `-b` from the command line, "
            "try also using `-rl`."
        )
        assert not (self.binary and self.label), err_msg
        assert self.outfile is not None, "Please provide an outfile location."
        if self.binary:
            np.save(self.outfile, self.counts)
        elif self.label:
            if names is None:
                names = seekrReader(self.infasta).get_headers()
            df = pd.DataFrame(data=self.counts, index=names, columns=self.kmers)
            df.to_csv(self.outfile)
        else:
            np.savetxt(self.outfile, self.counts, delimiter=",", fmt="%1.6f")

    def make_count_file(self, names=None):
        """Wrapper function for the most common way to generate count files.

        Given a numpy file name, it will save a numpy file where counts have been:
        cast as a dense array, centered, and standardized.

        Parameters
        ----------
        names : [str] (default=None)
            lncRNA names to pass to self.save

        Returns
        -------
        counts: np.array
            Matrix of kmer counts. Dimensions are equal to # of transcripts by # of kmers.
        """
        self.get_counts()
        if self.outfile is not None:
            self.save(names)
        return self.counts


"""
Description
-----------
Generate a matrix of Pearson similarities from two kmer count files.

Examples
--------
Generate a small, plain text .csv file:

```
import pandas as pd
from seekr.pearson import pearson

dist = pearson(counts1=example_2mers,
               counts2=example_2mers)
pd.DataFrame(dist).to_csv('example_vs_example_2mers.csv')

Notes
-----

Issues
------
Any issues can be reported to https://github.com/CalabreseLab/seekr/issues

---
"""


def seekrPearson(counts1, counts2, row_standardize=True, outfile=None):
    """Calculates a column standardized Pearson correlation matrix"""
    if row_standardize:
        counts1 = (counts1.T - np.mean(counts1, axis=1)).T
        counts1 = (counts1.T / np.std(counts1, axis=1)).T
        counts2 = (counts2.T - np.mean(counts2, axis=1)).T
        counts2 = (counts2.T / np.std(counts2, axis=1)).T

    # Take the inner product and save
    dist = np.inner(counts1, counts2) / counts1.shape[1]
    if outfile:
        np.save(outfile, dist)
    return dist


####################################################################################################################
### Description: 
# Find the best fitted distribution to the input background sequences 

### Details:
# data to be fitted -- all possible pairwise kmer pearson correlation for the background sequences
# can also choose to return the actual data by itself
# the results (best fitted distribution or actual data) will be used to calculate p-values
####################################################################################################################

# Finds and returns distributions and respective parameters that best fit input sequences
def find_dist(inputseq='default', k_mer=4, log2='Log2.post', models='common10', subsetting=True, subset_size = 100000, fit_model=True, statsmethod='ks',progress_bar=False, plotfit=None, outputname=None):

    #Determine what background transcriptome to use
    if inputseq == 'default':
        print('Using default background sequences: mouse vM25 Long non-coding RNA transcript sequences from Gencode.')
        print('https://www.gencodegenes.org/mouse/release_M25.html')
        print('Only including unique sequences and have the full length Airn(chr17:12741398-12830151) appended.')
        # Get the directory where the current file is located
        current_dir = os.path.dirname(os.path.realpath(__file__))
        
        # Construct the path to the .fa file
        inputseq = os.path.join(current_dir, 'data', 'gencode.vM25.lncRNA_transcripts.unique.genesequence_withfullairn.fa')
    
    #preprare the list of distributions to be used
    if models == 'common10':
        distributions = ['cauchy', 'chi2', 'expon', 'exponpow', 'gamma', 
                         'lognorm', 'norm', 'pareto', 'rayleigh', 't']
        
    else: 
        #Create all distributions (continous and discrete) that would be used
        continuous_distributions = [d for d in dir(stats) if
                        isinstance(getattr(stats, d), stats.rv_continuous)]

        discrete_distributions = [d for d in dir(stats) if
                                isinstance(getattr(stats, d), stats.rv_discrete)]

        distributions = continuous_distributions + discrete_distributions 

        # it is a convention to start the name of private or special methods with an underscore
        # this line is used to filter out private or special methods from the scipy.stats module
        distributions = [d for d in distributions if d[0] != '_']
        # remove levy_stable and studentized_range, it has problem converging
        # and the shape of the distribution is not what we want
        distributions = [d for d in distributions if d != 'levy_stable']
        distributions = [d for d in distributions if d != 'studentized_range']
        # if 'gilbrat' is in distributions, replace 'gilbrat' with 'gibrat' as 'gilbrat' is a misspelling of 'gibrat'
        # if 'gilbrat' in distributions:
        #     distributions[distributions.index('gilbrat')] = 'gibrat'
        # use what it is from the stats package

        if models != 'all':
            # check if the user input distributions are included in the distributions list
            orilen1=len(models)
            distributions = [d for d in models if d in distributions]
            if len(distributions) < orilen1:
                print("Please enter valid distribution names available in scipy.stats. refer to https://docs.scipy.org/doc/scipy/reference/stats.html#continuous-distributions")
                # find out which distributions are excluded
                excluded1 = [d for d in models if d not in distributions]
                print(f"Excluding invalid distributions for fitting: {excluded1}")
        

    # Convert distribution names to distribution objects
    tempdist = distributions.copy()
    orilen2=len(distributions)
    distributions = [getattr(stats, d) for d in distributions]

    # Filter distributions - some distributions can throw errors during the fit step
    # filter out the distributions that do not have a 'fit' method
    distributions = [d for d in distributions if 'fit' in dir(d)]
    if (len(distributions) < orilen2) and (models != 'all' and models != 'common10') :
        # find out which distributions are excluded
        excluded2 = [d for d in tempdist if d not in distributions]
        print(f"Excluding distributions do not have a 'fit' method: {excluded2}")


    #Create normalization mean and std vector for the inputseq
    bkg_norm_counter = seekrBasicCounter(inputseq,k=k_mer,log2=log2,silent=True)
    bkg_norm_counter.get_counts()
    mean_path = f'bkg_mean_{k_mer}mers.npy'
    std_path = f'bkg_std_{k_mer}mers.npy'
    np.save(mean_path, bkg_norm_counter.mean)
    np.save(std_path, bkg_norm_counter.std)

    # Count k-mers
    bkg_counter = seekrBasicCounter(inputseq,mean=mean_path,std=std_path,k=k_mer,silent=True)
    bkg_counter.make_count_file()

    # Find similarities
    sim_counts = seekrPearson(bkg_counter.counts,bkg_counter.counts)

    # Extract the upper triangle and flatten it
    sim_triu = sim_counts[np.triu_indices(sim_counts.shape[0], k=1)]

    #If the user wants to subset the data
    if subsetting: 
        if len(sim_triu)>subset_size:
            # randomly sample the user defined number of the sim_triu
            sim_triu = np.random.choice(sim_triu, size=subset_size, replace=False)
        else:
            print("subset_size is larger than the actual data size, use the actual data size instead")
    
    

    #if fitting chosen, fit the data to the distributions
    if fit_model:
        
        if len(distributions)>50 and len(sim_triu)>5000000 and subsetting==False:
            print("The input sequence count and distribution number for fitting are both large, subsetting is recommended to save time")
    
        results = []
        # set iterable based on progress_bar
        iterable = tqdm(distributions) if progress_bar else distributions
        # Fit distributions to data
        for distribution in iterable:
            # Try to fit the distribution
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                try: 
                    params = distribution.fit(sim_triu)

                    if statsmethod == 'ks':
                        # Get the KS test result on fitted data
                        D, p = kstest(sim_triu, distribution.name, args=params)

                    elif statsmethod == 'mse':
        
                        if isinstance(distribution, stats.rv_continuous):
                            synthetic_data = distribution.rvs(*params, size=len(sim_triu))
                        else:
                            synthetic_data = distribution.rvs(*params[:-2], loc=params[-2], scale=params[-1], size=len(sim_triu))
                            
                        residuals = sim_triu - synthetic_data
                        D = np.mean(residuals**2)


                    elif statsmethod == 'aic' or statsmethod == 'bic':
                        # Calculate log-likelihood
                        if isinstance(distribution, stats.rv_continuous):
                            log_likelihood = np.sum(distribution.logpdf(sim_triu, *params))
                        else:
                            log_likelihood = np.sum(distribution.logpmf(sim_triu, *params[:-2], loc=params[-2], scale=params[-1]))
                        
                        # Number of parameters
                        k = len(params)
                        # Number of data points
                        n = len(sim_triu)
                        
                        # Calculate AIC and BIC
                        if statsmethod == 'aic':
                            D = 2 * k - 2 * log_likelihood

                        else: 
                            D = np.log(n) * k - 2 * log_likelihood

                    else:
                        print("Please enter a valid statsmethod: 'ks', 'mse', 'aic', or 'bic'. Use default 'ks' now.")
                        D, p = kstest(sim_triu, distribution.name, args=params)

                # catch the error if the distribution cannot be fit, exclude it from the results
                except Exception as e:
                    print(f"Could not fit {distribution.name} because {e}, excluding it from the results")
                    continue

                # save the name and test result
                results.append((distribution.name, D, params))

        # Sort by minimum D statistics
        results.sort(key=lambda x: x[1])
        # sort by p val
        # results.sort(key=lambda x: x[2], reverse=True)

        if plotfit:
            n = len(results)
            # get the num of cols as the minimal of 5 and n
            n_cols = min(5, n)
            n_rows = n // n_cols + (n % n_cols > 0)   # Calculate the number of rows needed

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3, n_rows*3))  # Increase figure size as necessary
            axes = axes.ravel()  # Flatten the axes array

            # Generate PDF for each distribution and plot it
            for idx, (ax, result) in enumerate(zip(axes, results)):
                dist_name, Dval,params = result
                # print(dist_name)
                distribution = getattr(stats, dist_name)
                
                # Generate PDF
                x = np.linspace(min(sim_triu), max(sim_triu), 1000)
                pdf = distribution.pdf(x, *params)

                # Plot the original data as histogram and the fitted model
                ax.hist(sim_triu, bins=100, density=True, alpha=0.6, color='skyblue')
                ax.plot(x, pdf, 'r--', linewidth=2)
                ax.set_title(f'{idx+1}: {dist_name} (Dev={Dval:.3f})')

            # Remove unused subplots
            for i in range(len(results), len(axes)):
                fig.delaxes(axes[i])

            plt.tight_layout()
            # save plot
            plt.savefig(f'{plotfit}.pdf',dpi=300)

        # results is a list of tuples like (distribution_name, D_statistics, params)
        if outputname:
            # convert to dataframe and adding row/column names
            results_df = pd.DataFrame(results, columns=['distribution_name', 'D_statistics', 'params'])
            # save the dataframe to csv file
            results_df.to_csv(f'{outputname}.csv', index=False)

        return results
    else:
        if plotfit:
            print('No plot will be produced as fit_model is set to False, please set fit_model=True to plot the fitted distributions vs the actual data')
        
        if outputname:
            # save sim_triu as a csv file
            np.savetxt(f'{outputname}.csv', sim_triu, delimiter=",")

        return sim_triu


###################################################################################################
### Description: 
# calculte p values of the seekr.pearson correlation values for input sequence 1 vs input sequence 2 
# p value is based on the output of find_dist, which is either a list of distributions or a npy array

### Details:
# this function connects the output of find_dist and the input sequnces of interests (2 input sequences)
# given the background sequencs, find_dist calculate all possible pairwise seekr.pearson values
# and then outputs either a list of fitted distributions or the npy array of the actual data 
# find_pval firstly calculates the seekr.pearson of the two input sequences, which produces a correlation matrix 
# find_pval then calculate the p value for each r value in the pearson correlation matrix based on the output of find_dist
# the output of find_pval is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences
########################################################################################################

# check if fitres is a list consisting of tuples of (string, number, tuple of numbers)
# th number here could be float64 or float32
def is_float_type(x):
    return isinstance(x, float) or np.isscalar(x)

def check_tuple_format(tup):
    if len(tup) != 3:
        return False
    return isinstance(tup[0], str) and \
           is_float_type(tup[1]) and \
           isinstance(tup[2], tuple) and \
           all(is_float_type(x) for x in tup[2])

def check_main_list(main_list):
    return all(check_tuple_format(tup) for tup in main_list)


def find_pval(seq1file, seq2file, mean_path, std_path, k_mer, fitres, log2='Log2.post', bestfit=1, outputname=None, progress_bar=True):

    # first check whether input k_mer is compatible with the mean and std files
    meanfile = np.load(mean_path)
    stdfile = np.load(std_path)

    if len(meanfile) != 4**(k_mer) | len(stdfile) != 4**(k_mer):
        print('k_mer size is not compatible with the normalization mean and/or std files.')
        print('Please make sure the normalization mean and std files are generated using the same kmer size as specified here in k_mer.')
        print('No p value is calculated. The output is None.')
        return None
    
    t1 = seekrBasicCounter(seq1file, mean=mean_path, std=std_path, k=k_mer, log2=log2, silent=True)
    t2 = seekrBasicCounter(seq2file, mean=mean_path, std=std_path, k=k_mer, log2=log2, silent=True)

    t1.make_count_file()
    t2.make_count_file()

    sim = seekrPearson(t1.counts,t2.counts)

    header1=seekrReader(seq1file).get_headers()
    header2=seekrReader(seq2file).get_headers()
    # remove the '>' in the header
    header1 = [i[1:] for i in header1]
    header2 = [i[1:] for i in header2]

    # check if header1 only contains unique elements
    if len(header1) != len(set(header1)):
        print('The headers of seq1file is not unique.')
        print('Be carefule during further analysis as there are potential indexing problems.')

    # check if header2 only contains unique elements
    if len(header2) != len(set(header2)):
        print('The headers of seq2file is not unique.')
        print('Be carefule during further analysis as there are potential indexing problems.')
        
    # check if fitres is a distribution or a npy array
    if isinstance(fitres, list):
        # check if fitres is a list consisting of (string, number, tuple of numbers)
        if not check_main_list(fitres):
            print('The format of fitres is wrong.')
            print('fitres should be a list consisting of tuples (string, number, tuple of numbers) corresponds to (distribution name, deviance, parameters)')
            print('fitres should be the output of find_dist.')
            print('No p value is calculated. The output is None.')
            return None
        else: 
            # get the best fit distribution by user input
            # for python the counting starts from 0
            bestdist = fitres[bestfit-1]
            distname = bestdist[0]
            params = bestdist[2]

            #Load distribution
            dist = getattr(stats, distname)
            distribution = dist(*params)

            # Initialize a matrix to hold the p-values
            p_values = np.zeros_like(sim)

            # Loop through each element in sim
            iterable = tqdm(range(sim.shape[0]), desc="Rows") if progress_bar else range(sim.shape[0])
            for i in iterable:
                for j in range(sim.shape[1]):
                    p_values[i, j] = 1 - distribution.cdf(sim[i, j])
            # convert to dataframe and adding row/column names
            pval_df = pd.DataFrame(p_values, index=header1, columns=header2)
            # save the dataframe to csv file
            if outputname:
                pval_df.to_csv(f'{outputname}.csv')

            return pval_df



    elif isinstance(fitres, np.ndarray):
        # check if fitres is a 1D array
        if len(fitres.shape) == 1:
            # Calculate p-values for entire sim matrix
            total_length_fitres = len(fitres)

            # alternative method to the loop below, faster but requires more memory
            # Use broadcasting to compare each element in sim with all elements in fitres
            # The result is a 3D array where the third dimension has size equal to total_length_fitres
            # comparison_result = fitres > sim[:,:,np.newaxis]
            # Sum along the third dimension and divide by total_length_fitres to get p-values
            # p_values = np.sum(comparison_result, axis=-1) / total_length_fitres

            # Initialize a matrix to hold the p-values
            p_values = np.zeros_like(sim)

            # Loop through each element in sim
            iterable = tqdm(range(sim.shape[0]), desc="Rows") if progress_bar else range(sim.shape[0])
            for i in iterable:
                for j in range(sim.shape[1]):
                    p_values[i, j] = np.sum(fitres > sim[i, j]) / total_length_fitres

            # convert to dataframe and adding row/column names
            pval_df = pd.DataFrame(p_values, index=header1, columns=header2)
            # save the dataframe to csv file
            if outputname:
                pval_df.to_csv(f'{outputname}.csv')
            return pval_df 
        
        else:
            print('The dimension of fitres as a numpy array is wrong. fitres should be a 1D numpy array.')
            print('fitres should be the output of find_dist.')
            print('No p value is calculated. The output is None.')
            return None


    else:
        print('fitres should be the output of find_dist. It should be either a list of distributions or a numpy array.')
        print('No p value is calculated. The output is None.')
        return None





