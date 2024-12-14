# Preliminary Consensus Sequence QC - FASTQ

This Python script processes a FASTQ file by performing reverse complement correction, detecting outliers, and removing low-quality reads based on alignment scores. The final output includes a corrected FASTQ file with comprehensive quality control applied.

## Features

1. **Consensus Sequence Generation**  
   Generates a consensus sequence from the input FASTQ file by taking the most common base at each position.

2. **Outlier Detection**  
   Detects and removes:
   - **Length-based outliers**: Reads that are unusually short or long based on interquartile range (IQR) thresholds.
   - **Alignment-based outliers**: Reads with alignment scores below a threshold (default: 20).

3. **Reverse Complement Correction**  
   Aligns each read to the consensus sequence and applies reverse complement if the reverse alignment score is better.

4. **Worst 10% Removal**  
   Removes the bottom 10% of aligning reads based on alignment scores, but only if the input file contains 5 or more reads.

5. **Comprehensive Reporting**  
   Provides detailed on-screen output, including:
   - Total initial reads.
   - Number of length-based and alignment-based outliers removed.
   - Number of worst-aligning reads removed.
   - Total reads remaining after QC.

6. **Dynamic Output Filenames**  
   Corrected FASTQ file names include the total number of reads remaining after QC (e.g., `Corrected_InputFileName-XXseqs.fastq`).

---

## Usage

### Prerequisites

- Python 3.6 or higher
- Biopython library (`pip install biopython`)
- NumPy library (`pip install numpy`)

### How to Run

1. Place the script and your FASTQ file in the same directory.
2. Run the script:
   ```bash
   python script_name.py
