from Bio import SeqIO
from collections import Counter, defaultdict
import os
import numpy as np
import edlib
from typing import List, Tuple
import subprocess
import random
import tempfile
import argparse

def find_fastq_file():
    """Find the FASTQ file in the current directory."""
    for file in os.listdir("."):
        if file.endswith(".fastq"):
            return file
    raise FileNotFoundError("No FASTQ file found in the current directory.")


def find_representative(sequences: List[SeqIO.SeqRecord], sample_size: int = 100, alignment_threshold: float = 0.4) -> \
        Tuple[SeqIO.SeqRecord, bool]:
    """
    Find a representative sequence by comparing a sample of sequences.
    Returns (sequence, should_reverse)
    """
    sample_size = min(sample_size, len(sequences))
    sample = sequences[:sample_size]

    # Calculate distance matrix including reverse complements
    distances = defaultdict(dict)

    for i, seq1 in enumerate(sample):
        seq1_str = str(seq1.seq)
        seq1_rc = str(seq1.seq.reverse_complement())

        for j, seq2 in enumerate(sample[i + 1:], i + 1):
            seq2_str = str(seq2.seq)

            # Compare all orientations
            d1 = edlib.align(seq1_str, seq2_str, task="distance")["editDistance"]
            d2 = edlib.align(seq1_rc, seq2_str, task="distance")["editDistance"]

            # Store normalized distances
            max_len = max(len(seq1_str), len(seq2_str))
            distances[i][j] = min(d1, d2) / max_len
            distances[j][i] = distances[i][j]

    best_count = 0
    best_idx = -1

    for i in range(sample_size):
        count = 0
        for j in range(sample_size):
            if j != i and distances[i][j] < alignment_threshold:
                count += 1

        if count > best_count:
            best_count = count
            best_idx = i

    return sample[best_idx], False


def detect_length_outliers(sequences):
    """
    Detect length-based outliers using the interquartile range (IQR).
    """
    lengths = [len(record.seq) for record in sequences]
    q1, q3 = np.percentile(lengths, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    length_outliers = [i for i, record in enumerate(sequences) if
                       len(record.seq) < lower_bound or len(record.seq) > upper_bound]
    return length_outliers, lower_bound, upper_bound


def generate_spoa_consensus(sequences, sample_size=100):
    """Run SPOA to generate consensus sequence."""
    if not sequences:
        return None

    try:
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n{str(seq.seq)}\n")
            temp_input = f.name

        # Construct SPOA command with parameters
        cmd = [
            "spoa",
            temp_input,
            "-r", "0",  # Result mode 0: consensus only
            "-l", "1",  # Global alignment mode
            "-m", "5",  # Match score
            "-n", "-4",  # Mismatch penalty
            "-g", "-8",  # Gap opening penalty
            "-e", "-6",  # Gap extension penalty
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        # Clean up input file
        os.unlink(temp_input)

        # Parse SPOA output for consensus sequence
        consensus = None
        for line in result.stdout.split('\n'):
            if not line.startswith('>') and line.strip():
                consensus = line.strip()
                break

        if not consensus:
            print("SPOA did not generate consensus sequence")
            return sequences[0]  # Fall back to first sequence

        # Create a SeqRecord object for the consensus
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        consensus_record = SeqRecord(Seq(consensus), id="consensus", description="spoa_consensus")

        return consensus_record

    except subprocess.CalledProcessError as e:
        print(f"SPOA failed with return code {e.returncode}")
        print(f"Command: {' '.join(cmd)}")
        print(f"Stderr: {e.stderr}")
        return sequences[0]  # Fall back to first sequence

    except Exception as e:
        print(f"Error running SPOA: {str(e)}")
        return sequences[0]  # Fall back to first sequence

def calculate_alignment_distances(sequences, reference):
    """
    Calculate alignment distances for each sequence against the reference using edlib.
    """
    ref_str = str(reference.seq)
    distances = []

    for i, record in enumerate(sequences):
        seq_str = str(record.seq)
        # Get edit distance
        result = edlib.align(seq_str, ref_str, task="distance")
        # Normalize by max length
        max_len = max(len(seq_str), len(ref_str))
        norm_distance = result["editDistance"] / max_len
        distances.append((i, record, norm_distance))

    return distances


def calculate_dynamic_threshold(distances, percentile=90):
    """
    Calculate a dynamic alignment distance threshold based on the top percentile.
    Note: For edit distances, higher values indicate worse alignment.
    """
    return np.percentile([d for _, _, d in distances], percentile)


def detect_alignment_outliers(sequences, alignment_distances, threshold):
    """
    Detect alignment-based outliers using a dynamic threshold.
    For edit distances, sequences with distance > threshold are outliers.
    """
    outliers = [i for i, _, distance in alignment_distances if distance > threshold]
    return outliers


def align_and_reverse_complement(sequences, reference):
    """
    Align each sequence to a reference and reverse complement if it aligns better in reverse orientation.
    Uses edlib for alignment.
    """
    ref_str = str(reference.seq)
    reverse_complemented = 0
    updated_sequences = []

    print("Aligning sequences and checking for reverse complements...")
    for record in sequences:
        seq_str = str(record.seq)
        seq_rc = str(record.seq.reverse_complement())

        forward_dist = edlib.align(seq_str, ref_str, task="distance")["editDistance"]
        reverse_dist = edlib.align(seq_rc, ref_str, task="distance")["editDistance"]

        if reverse_dist < forward_dist:
            print(f"Sequence {record.id} reverse complemented.")
            record.seq = record.seq.reverse_complement()
            reverse_complemented += 1

        updated_sequences.append(record)

    return updated_sequences, reverse_complemented


def save_fastq(sequences, output_file):
    """Save sequences to a new FASTQ file."""
    print(f"Saving sequences to {output_file}...")
    with open(output_file, "w") as f:
        SeqIO.write(sequences, f, "fastq")


def main():
    try:
        parser = argparse.ArgumentParser(description='Preliminary Quality Control Tool')
        parser.add_argument('-c', '--consensus', action='store_true',
                            help='Generate a consensus prior to calculating alignment outliers')
        args = parser.parse_args()

        # Step 1: Find and load FASTQ file
        fastq_file = find_fastq_file()
        print(f"Processing file: {fastq_file}")

        print("Parsing FASTQ sequences...")
        sequences = list(SeqIO.parse(fastq_file, "fastq"))

        iteration = 1
        convergence = False
        reverse_complemented_total = 0

        representative = None
        while not convergence:
            print(f"\nIteration {iteration}: Finding representative sequence and aligning sequences...")

            # Step 2: Find representative sequence
            representative, _ = find_representative(sequences)
            print(f"Representative sequence selected (first 50 bases): {str(representative.seq)[:50]}...")

            # Step 3: Reverse complement alignment
            updated_sequences, reverse_complemented = align_and_reverse_complement(sequences, representative)
            reverse_complemented_total += reverse_complemented

            # Check for convergence
            if reverse_complemented == 0:
                convergence = True
            else:
                sequences = updated_sequences

            iteration += 1

        # Step 4: Detect and remove length-based outliers
        print("Detecting length-based outliers...")
        length_outliers, lower_bound, upper_bound = detect_length_outliers(sequences)
        outliers_length_removed = len(length_outliers)
        sequences = [seq for i, seq in enumerate(sequences) if i not in length_outliers]
        print(f"Length-based outliers removed: {outliers_length_removed}")
        print(f"Length thresholds: lower = {lower_bound}, upper = {upper_bound}")

        if args.consensus:
            # Step 5: Generate consensus using spoa
            print("Generating consensus sequence using spoa...")
            consensus = generate_spoa_consensus(sequences, sample_size=100)
            print(f"Consensus sequence generated (first 50 bases): {str(consensus.seq)[:50]}...")
        else:
            consensus = representative

        # Step 6: Calculate alignment distances and detect alignment-based outliers using the spoa consensus
        print("Calculating alignment distances and detecting alignment-based outliers...")
        alignment_distances = calculate_alignment_distances(sequences, consensus)
        threshold = calculate_dynamic_threshold(alignment_distances, percentile=90)
        alignment_outliers = detect_alignment_outliers(sequences, alignment_distances, threshold)
        outliers_alignment_removed = len(alignment_outliers)
        sequences = [seq for i, seq in enumerate(sequences) if i not in alignment_outliers]
        print(f"Alignment-based outliers removed: {outliers_alignment_removed}")

        # Step 7: Save corrected sequences
        final_read_count = len(sequences)
        corrected_file = f"Corrected_{os.path.basename(fastq_file).replace('.fastq', '')}-{final_read_count}seqs.fastq"
        save_fastq(sequences, corrected_file)

        # Final report
        print("\nProcessing complete!")
        print(f"Total initial sequences: {len(sequences) + outliers_length_removed + outliers_alignment_removed}")
        print(f"Total length-based outliers removed: {outliers_length_removed}")
        print(f"Total alignment-based outliers removed: {outliers_alignment_removed}")
        print(f"Total reads remaining after preliminary QC: {final_read_count}")
        print(f"Total sequences reverse complemented: {reverse_complemented_total}")
        print(f"Corrected sequences saved to: {corrected_file}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()