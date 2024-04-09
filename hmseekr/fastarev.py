###################################################################################################
### Description: 
# This function reverse a fasta file and save it to a new file
# 'ATGC' will be reversed to 'CGTA'

### Details:
# this function takes in fasta file and reverse the sequence
# if input fasta contains more than one sequence, each sequence will be reversed
# the headers will be kept the same, only the sequences will be reversed
# keeping the headers the same will allow the user to match the reversed sequences to the original sequences easily
# the output is a fasta file with the reversed sequences


### Input:
# input_file_path: path to the input fasta file. 
# output_file_path: path and name to the output fasta file.


### Output:
# a fasta file with the reversed sequences and the SAME headers as the input fasta file
# if the input fasta file contains more than one sequence, each sequence will be reversed

### Example:
# from hmseekr.fastarev import fastarev

# fastarev(input_file_path='../fastaFiles/mXist_rA.fa',
#          output_file_path='../fastaFiles/mXist_rA_rev.fa')


########################################################################################################


def fastarev(input_file_path, output_file_path):
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        for line in input_file:
            if line.startswith('>'):  # This is a header
                output_file.write(line)
            else:  # This is a sequence
                reversed_sequence = line.strip()[::-1]  # Reverse the sequence
                output_file.write(reversed_sequence + '\n')

