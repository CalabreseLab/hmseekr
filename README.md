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


### findhits: find high similar regions based on kmer profile within sequences of interest

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


<hr/>

## hmseekr Docker Image
We now provide a Docker image of hmseekr which hopefully avoids all the package dependencies and complications when install and run seekr through pip. 

### Docker Installation
Firstly, you should install Docker on your computer. Please refer to the [official website](https://www.docker.com/get-started/) for installation instructions.

After you have successfully installed Docker. Start/Run the application and make sure it is running properly – you should see the docker icon on your task bar with the status indicated.

### Pull Docker Image and Test Run
1. Start your command line tool: Terminal for MacOS and CMD for Windows. You can also use Powershell or Cygwin for Windows, but Cygwin might have interaction issues.

From the command line, pull the Docker Image:
```
$ docker pull calabreselab/hmseekr:latest
```
You can replace `latest` with a specific tag if needed.

2. Test Run the Docker Image
```
$ docker run -it --rm calabreselab/hmseekr:latest
```
The `-it` tag sets it to interactive mode. If you don't need to run the Docker container in interactive mode (i.e., you don't need a shell inside the container), you can simply omit the `-it` flag.
This will print the user manual out to the command line, which is basically the same as you run the command `seekr` directly from command line when you pip install seekr. 

### Run Docker Image from command line
You can run the seekr function from this Docker Image directly from command line with the following specified syntax.
```
$ docker run -v /path/to/your/files:/data calabreselab/hmseekr:latest hmseekr_kmers -fd '/data/fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -dir '/data/counts/' -name repeatA
```
In this command:
* `-v /path/to/your/files:/data`: This mounts the directory `/path/to/your/files` from your host machine (where /fastaFiles/mXist_rA.fa is located) to the `/data` directory inside the Docker container. Replace `/path/to/your/files` with the actual path to your files.
* `mseekr_kmers -fd '/data/fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -dir '/data/counts/' -name repeatA`: This is the command that gets executed inside the Docker container. Since we mounted our files to `/data` in the container, we reference them with `/data/fastaFiles/mXist_rA.fa`. 
* The `/data` folder is basically a mirror of the folder you specified in `/path/to/your/files`. So by specifying `-dir '/data/counts/'` (output into the counts folder under the path /data/) we can have the output files directly written in to the folder in `/path/to/your/files`.
* Please remember to **specify your output path to `/data/`** otherwise it will not be saved to your folder on local machine and it would be hard to locate it even inside the Docker Container Filesystem (in this case, when the Docker Container is removed, your data will be deleted as well). 

Examples of code mounts e:/test on Windows as the folder that contains the input and holds the output files:
```
$ docker run -v e:/test:/data calabreselab/hmseekr:latest hmseekr_kmers -fd '/data/fastaFiles/mXist_rA.fa' -k 2,3,4 -a ATCG -dir '/data/counts/' -name repeatA
```
Basically you need to add: `docker run -v /path/to/your/files:/data calabreselab/hmseekr:latest` before the command line code for hmseekr (see above for examples of all functions).

### Run Docker Image with Jupyter Notebook
If you want to work with python codes instead of directly calling from the command line, you can choose to run seekr with Jupyter Notebook inside the Docker Container.

1.  Run Docker Container with Interactive Terminal:
```
$ docker run -it -p 8888:8888 -v /path/on/host:/path/in/container calabreselab/hmseekr:latest /bin/bash
```
This command will start the container and give you a bash terminal inside it. The `-p 8888:8888` flag maps the port *8888* from the container to your host machine so that you can access the Jupyter Notebook server.

`/path/on/host` is the path to a directory on your local machine where you want the notebooks and other data to be saved. `/path/in/container` is the directory inside the container where you want to access and save the files from Jupyter Notebook.

When you use Jupyter Notebook now and create or save notebooks, they will be stored at `/path/in/container` inside the Docker container. Due to the volume mount, these files will also be accessible and stored at `/path/on/host` on your host machine so that you can later access the code and files even when the container is removed. 

Example of code:
```
$ docker run -it -p 8888:8888 -v e:/test:/data calabreselab/hmseekr:latest /bin/bash
```

2. Manually start Jupyter Notebook. From the bash terminal inside the Docker container:
```
$ jupyter notebook --ip=0.0.0.0 --port=8888 --NotebookApp.token='' --NotebookApp.password='' --allow-root
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
from hmseekr import kmers

testdict = kmers(fadir='./fastaFiles/mXist_rA.fa',kvec='4',
                 alphabet='ATCG',outputdir='./counts/',
                 outputname='repeatA')
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

