# installation

pip install hmseekr

# install from github repo
pip install git+https://github.com/CalabreseLab/hmseekr.git

####################################
# working example use XISTrB1 as the query and all v47 canonical lncRNA that has length > 500nt and depulicated as null and background sequences. use a user defined sequences as pool

# 1
# count kmer profiles
# need to count kmers for both query and null 

# query
hmseekr_kmers -fd 'XISTrB1.fa' -k 4 -a ATCG -name repeatB1 -dir './counts/' 

# null
hmseekr_kmers -fd 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -a ATCG -name v47lnc -dir './counts/' 


# 2
# use gridsearch to optimize transition probabilities

hmseekr_gridsearch -qf 'XISTrB1.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_user_defined_seqs.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 800 -name 'XISTrB1' -dir './gridsearch_XISTrB1/' -a 'ATCG'

# the output folder gridsearch_XISTrB1 will contain some intermediate outputs in subfolders such as hits for all combinations of qT and nT
# the XISTrB1.csv is the table that contains all the stats with each row being one combination of qT and nT, and listing the total number of filtered hits sequences and the median, 
# standard deviation of the filtered hits sequences' pearson correlation r score to the query sequence 
# and the median, standard deviation of the length of the filtered hits sequences
# and the same stats for the top 50 filtered hits sequences (ranked by seekr r score) if there are more than 50 hits
# filtered hits are filtered by lenmin (-li) and lenmax (-la)


# 3
# train Markov Models
# using optimized transition probability from gridsearch results: qT=0.9 and nT=0.93

hmseekr_train -qd './counts/repeatB1.dict' -nd './counts/v47lnc.dict' -k 4 -a ATCG -qT 0.9 -nT 0.93 -qPre XISTrB1 -nPre v47lnc -dir './markovModels/'


# 4
# search for hits in the user defined sequence pool

hmseekr_findhits_condE -pool 'v47_user_defined_seqs.fa' -m './markovModels/XISTrB1_v47lnc/4/hmm.dict' -k 4 -name XISTrB1_v47_search -dir './hits/' -a 'ATCG' -fa -pb

# 5 
# add SEEKR p value to the results but do no further filtering

hmseekr_hitseekr -hd './hits/XISTrB1_v47_search_4_viterbi.txt' -qf 'XISTrB1.fa.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -li 0 -la 1000000 -pf 1.1 -rf -1.1 -dir './hits/' -pb

####################################
# quanitfy the similarity to several queries (XIST all repeats each as a fasta entry, but saved in one file) at the same time

hmseekr_seqstosummary -qf 'XIST_allrepeats.fa' -td 'optimized_transit_p_XIST_allrepeats.csv' -lf 'lenfilter_XIST_allrepeats.csv' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_user_defined_seqs.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -fc 'findhits_condE' -pf 0.01 -name 'XISTallrepeats_v47' -dir './XISTallrepeats_seqstosummary/' -format wide -a 'ATCG' -pb

# please refer to this website for the details of each arguments
# https://github.com/CalabreseLab/hmseekr





