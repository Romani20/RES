# Relative Evolutionary Scoring (RES) algorithm

The RES algorithm uses a position weight matrix and background frequencies to score each
position of individual aligned sequences by taking the log odds ratio of observed frequency
to background frequency. It gives insights into several protein evolutionary features such as 
degree of conservation, divergence, and sites of functional importance. 

## Prerequisites

We provided files to test the results shown in our related paper, "Structural and mechanistic diversity
in p53-mediated regulation of organismal longevity across taxonomical orders." The algorithm can be
run on any set of closely related aligned sequences. The gap opening and extension penalties can be adjusted
based on sequence relatedness.

To set up the Python environment for run, PIP install:
 - argparse
 - math
 - pandas
 - Bio

## Running algorithm

Run: python res.py <path_to_aligned_.fasta_file> <name_of_new_output_.xlsx_file>



