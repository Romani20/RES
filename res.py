import argparse
import math
import pandas as pd
from Bio import SeqIO

vals = 'ACDEFGHIKLMNPQRSTVWY-'

def fasta_to_array(fasta_file):
    """Converts a the sequences in an aligned fasta file 
    to an array of strings.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def calculate_pwm(sequences, pseudocount=1):
    """Calculates a position weight matrix from array of aligned 
    sequences.
    """
    counts = {residue: [0] * len(sequences[0]) for residue in vals}

    for sequence in sequences:
        for i, residue in enumerate(sequence):
            if residue in counts:
                counts[residue][i] += 1
            else:
                print(f"Invalid character '{residue}' in sequence. Skipping.")

    for residue in counts:
        for i in range(len(counts[residue])):
            counts[residue][i] += pseudocount

    pwm = {residue: [count / (len(sequences) + pseudocount * len(
        counts)) for count in counts[residue]] for residue in counts}
    return pwm

def calculate_background_frequency(sequences):
    """Calculates background frequency of each residue (or gap)
    over the aligned sequences.
    """
    dictbg = {val: 0 for val in vals}
    dictbg["total"] = 0

    for line in sequences:
        for i in line:
            if i in dictbg:
                dictbg[i] += 1
                dictbg["total"] += 1
            else:
                print(f"Invalid character '{i}' in sequence. Skipping.")

    dictf = {i: dictbg[i] / dictbg["total"] for i in vals}
    return dictf

def calculate_position_score(
        sequence, pwm, background_freq, gap_open_penalty=2.90, gap_extend_penalty=0):
    """Calculates the position-specific score for a sequence 
    using pwm and background frequency.
    """
    num_positions = len(sequence)
    scores = []

    gap_opened = False

    for i in range(num_positions):
        score = 0.0

        residue = sequence[i]
        if residue == '-':
            if not gap_opened:
                score -= gap_open_penalty
                gap_opened = True
            else:
                score -= gap_extend_penalty
        else:
            gap_opened = False
            if residue in pwm and i < len(pwm[residue]) and residue in background_freq:
                if pwm[residue][i] > 0 and background_freq[residue] > 0:
                    score += math.log2(pwm[residue][i] / background_freq[residue])

        scores.append(score)
        print(f"Score for position {i + 1}: {score}")

    return scores

def calculate_and_save_scores_to_excel(sequences, output_file):
    """Calculate and write scores to Excel file, where each column
    represents scores for an ortholog, and rows contain residue-wise 
    scores.
    """
    pwm = calculate_pwm(sequences)
    background_freq = calculate_background_frequency(sequences)

    df = pd.DataFrame()
    for i, sequence in enumerate(sequences):
        sequence_scores = calculate_position_score(sequence, pwm, background_freq)
        df[f'Sequence_{i+1}_scores'] = sequence_scores
    df.to_excel(output_file, index=False)
    print(f"Scoring completed.")

def main():
    parser = argparse.ArgumentParser(description='Calculate scores and write to Excel.')
    parser.add_argument('fasta_file', type=str, help='Path to fasta file of aligned sequences')
    parser.add_argument('output_file', type=str, help='Path to the output Excel file')
    
    args = parser.parse_args()
    sequences = fasta_to_array(args.fasta_file)
    calculate_and_save_scores_to_excel(sequences, args.output_file)

if __name__ == '__main__':
    main()
