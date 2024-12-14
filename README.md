# Preliminary Consensus Sequence QC - FASTQ

This Python script processes a FASTQ file to:
1. Detect and remove outliers based on read length and alignment scores.
2. Correct the orientation of reads using reverse complement alignment.
3. Output a cleaned FASTQ file with comprehensive quality control applied.

The script dynamically calculates thresholds for alignment-based outliers using percentile-based scoring and for length-based outliers using interquartile range (IQR).

---

## Features

### 1. Reverse Complement Correction
- Aligns each read to the consensus sequence.
- Corrects the orientation of reads by applying reverse complement if the reverse alignment score is better.


### 2. Outlier Detection
- **Length-Based Outliers**:
  - Detects unusually short or long reads based on IQR thresholds.
  - Automatically removes these reads from the dataset.
- **Alignment-Based Outliers**:
  - Calculates dynamic thresholds based on the bottom 10% of alignment scores.
  - Flags and removes reads with alignment scores below this threshold.
  - Does not implement additional filtering at this step if there are fewer than 5 reads.

### 3. Comprehensive Reporting
- Provides detailed on-screen output for:
  - Total initial reads.
  - Length-based and alignment-based outliers removed.
  - Total reverse complement corrections applied.
  - Final number of reads after QC.

### 4. Dynamic Output Filenames
- The corrected FASTQ file includes the total number of reads remaining in the filename (e.g., `Corrected_InputFileName-XXseqs.fastq`).

---

## Usage

### Prerequisites

- Python 3.6 or higher
- Required Python libraries:
  - **Biopython**: Install using `pip install biopython`
  - **NumPy**: Install using `pip install numpy`

### How to Run

1. Place the script and your FASTQ file in the same directory.
2. Run the script:
   ```bash
   python script_name.py
