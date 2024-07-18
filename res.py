import argparse
import math
import pandas as pd
from Bio import SeqIO

vals = 'ACDEFGHIKLMNPQRSTVWY-'

def fasta_to_array(fasta_file):
    """_summary_

    Args:
        fasta_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def scores_to_dict(scores):
    """_summary_

    Args:
        scores (_type_): _description_

    Returns:
        _type_: _description_
    """
    aa_score_dict = {val: score for val, score in zip(vals, scores)}
    return aa_score_dict

def calculate_pwm(sequences, pseudocount=1):
    """_summary_

    Args:
        sequences (_type_): _description_
        pseudocount (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
    counts = {nucleotide: [0] * len(sequences[0]) for nucleotide in vals}

    for sequence in sequences:
        for i, nucleotide in enumerate(sequence):
            if nucleotide in counts:
                counts[nucleotide][i] += 1
            else:
                print(f"Warning: Found invalid character '{nucleotide}' in sequence. Skipping.")

    for nucleotide in counts:
        for i in range(len(counts[nucleotide])):
            counts[nucleotide][i] += pseudocount

    pwm = {nucleotide: [count / (len(sequences) + pseudocount * len(counts)) for count in counts[nucleotide]] for nucleotide in counts}
    return pwm

def calculate_background_frequency(sequences):
    """_summary_

    Args:
        sequences (_type_): _description_

    Returns:
        _type_: _description_
    """
    dictbg = {val: 0 for val in vals}
    dictbg["total"] = 0

    for line in sequences:
        for i in line:
            if i in dictbg:
                dictbg[i] += 1
                dictbg["total"] += 1
            else:
                print(f"Warning: Found invalid character '{i}' in sequence. Skipping.")

    dictf = {i: dictbg[i] / dictbg["total"] for i in vals}
    return dictf

def calculate_position_score(sequence, pwm, background_freq, gap_open_penalty=14, gap_extend_penalty=3):
    """_summary_

    Args:
        sequence (_type_): _description_
        pwm (_type_): _description_
        background_freq (_type_): _description_
        gap_open_penalty (int, optional): _description_. Defaults to 14.
        gap_extend_penalty (int, optional): _description_. Defaults to 3.

    Returns:
        _type_: _description_
    """
    num_sequences = len(sequence)
    scores = []

    for i in range(num_sequences):
        score = 0.0
        gap_opened = False 

        nucleotide = sequence[i]
        if nucleotide == '-':
            if not gap_opened:
                score -= gap_open_penalty
                gap_opened = True  
            else:
                score -= gap_extend_penalty
        else:
            gap_opened = False
            if nucleotide in pwm and nucleotide in background_freq:
                if pwm[nucleotide][i] > 0 and background_freq[nucleotide] > 0:
                    score += math.log2(pwm[nucleotide][i] / background_freq[nucleotide])
        scores.append(score)
        print(f"Score for segment {i+1}: {score}")
    return scores


def calculate_and_save_scores_to_excel(sequences, output_file):
    """_summary_

    Args:
        sequences (_type_): _description_
        output_file (_type_): _description_
    """
    pwm = calculate_pwm(sequences)
    background_freq = calculate_background_frequency(sequences)

    df = pd.DataFrame()
    for i, sequence in enumerate(sequences):
        sequence_scores = calculate_position_score(sequence, pwm, background_freq)
        df[f'Sequence_{i+1}_scores'] = sequence_scores
    df.to_excel(output_file, index=False)
    print(f"Scores saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Calculate scores for sequences and save to an Excel file.')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file containing sequences')
    parser.add_argument('output_file', type=str, help='Path to the output Excel file')
    
    args = parser.parse_args()
    sequences = fasta_to_array(args.fasta_file)
    calculate_and_save_scores_to_excel(sequences, args.output_file)

if __name__ == '__main__':
    main()
