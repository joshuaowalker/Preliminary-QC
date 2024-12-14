from Bio import SeqIO
from Bio.Align import PairwiseAligner
from collections import Counter
import os
import numpy as np

def find_fastq_file():
    """
    Find the FASTQ file in the current directory.
    """
    for file in os.listdir("."):
        if file.endswith(".fastq"):
            return file
    raise FileNotFoundError("No FASTQ file found in the current directory.")

def generate_consensus(sequences):
    """
    Generate a simple consensus sequence by taking the most common base at each position.
    """
    print("Generating consensus sequence...")
    length = max(len(record.seq) for record in sequences)
    consensus = []
    for i in range(length):
        bases = [record.seq[i] for record in sequences if i < len(record.seq)]
        most_common_base = Counter(bases).most_common(1)[0][0]
        consensus.append(most_common_base)
    return "".join(consensus)

def detect_length_outliers(sequences):
    """
    Detect length-based outliers using the interquartile range (IQR).
    """
    lengths = [len(record.seq) for record in sequences]
    q1, q3 = np.percentile(lengths, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    length_outliers = [i for i, record in enumerate(sequences) if len(record.seq) < lower_bound or len(record.seq) > upper_bound]
    return length_outliers, lower_bound, upper_bound

def calculate_alignment_scores(sequences, reference):
    """
    Calculate alignment scores for each sequence against the consensus.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    scores = [(i, record, aligner.score(record.seq, reference)) for i, record in enumerate(sequences)]
    return scores

def calculate_dynamic_threshold(scores, percentile=10):
    """
    Calculate a dynamic alignment score threshold based on the bottom percentile.
    """
    return np.percentile(scores, percentile)

def detect_alignment_outliers(sequences, alignment_scores, threshold):
    """
    Detect alignment-based outliers using a dynamic threshold.
    """
    outliers = [i for i, _, score in alignment_scores if score < threshold]
    return outliers

def align_and_reverse_complement(sequences, reference):
    """
    Align each sequence to a reference and reverse complement if it aligns better in reverse orientation.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"

    reverse_complemented = 0
    updated_sequences = []

    print("Aligning sequences and checking for reverse complements...")
    for record in sequences:
        forward_score = aligner.score(record.seq, reference)
        reverse_score = aligner.score(record.seq.reverse_complement(), reference)

        if reverse_score > forward_score:
            print(f"Sequence {record.id} reverse complemented.")
            record.seq = record.seq.reverse_complement()
            reverse_complemented += 1

        updated_sequences.append(record)

    return updated_sequences, reverse_complemented

def save_fastq(sequences, output_file):
    """
    Save sequences to a new FASTQ file.
    """
    print(f"Saving sequences to {output_file}...")
    with open(output_file, "w") as f:
        SeqIO.write(sequences, f, "fastq")

def main():
    try:
        # Step 1: Find and load FASTQ file
        fastq_file = find_fastq_file()
        print(f"Processing file: {fastq_file}")
        
        print("Parsing FASTQ sequences...")
        sequences = list(SeqIO.parse(fastq_file, "fastq"))

        # Step 2: Detect length-based outliers
        print("Detecting length-based outliers...")
        length_outliers, lower_bound, upper_bound = detect_length_outliers(sequences)
        outliers_length_removed = len(length_outliers)
        sequences = [seq for i, seq in enumerate(sequences) if i not in length_outliers]
        print(f"Length-based outliers removed: {outliers_length_removed}")
        print(f"Length thresholds: lower = {lower_bound}, upper = {upper_bound}")

        # Step 3: Generate consensus sequence
        consensus = generate_consensus(sequences)
        print(f"Consensus sequence generated (first 50 bases): {consensus[:50]}...")

        # Step 4: Calculate alignment scores
        print("Calculating alignment scores...")
        alignment_scores = calculate_alignment_scores(sequences, consensus)
        scores = [score for _, _, score in alignment_scores]  # Extract scores

        # Step 5: Calculate dynamic threshold and detect alignment-based outliers
        threshold = calculate_dynamic_threshold(scores, percentile=10)
        print(f"Dynamic alignment score threshold: {threshold:.2f}")
        alignment_outliers = detect_alignment_outliers(sequences, alignment_scores, threshold)
        outliers_alignment_removed = len(alignment_outliers)
        sequences = [seq for i, seq in enumerate(sequences) if i not in alignment_outliers]
        print(f"Alignment-based outliers removed: {outliers_alignment_removed}")

        # Step 6: Reverse complement alignment
        updated_sequences, reverse_complemented = align_and_reverse_complement(sequences, consensus)

        # Step 7: Save corrected sequences
        final_read_count = len(updated_sequences)
        corrected_file = f"Corrected_{os.path.basename(fastq_file).replace('.fastq', '')}-{final_read_count}seqs.fastq"
        save_fastq(updated_sequences, corrected_file)

        # Final report
        print("\nProcessing complete!")
        print(f"Total initial sequences: {len(alignment_scores) + outliers_length_removed}")
        print(f"Total length-based outliers removed: {outliers_length_removed}")
        print(f"Total alignment-based outliers removed: {outliers_alignment_removed}")
        print(f"Total reads remaining after preliminary QC: {final_read_count}")
        print(f"Total sequences reverse complemented: {reverse_complemented}")
        print(f"Corrected sequences saved to: {corrected_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
