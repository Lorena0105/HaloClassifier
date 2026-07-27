# HaloClassifier

**Current release:** v1.0.1

HaloClassifier is a machine-learning tool for classifying haloarchaeal contigs as chromosomal or plasmidic.

It uses sequence-derived features and size-specific Random Forest models in a two-stage pipeline.

---

## Important

The pre-trained **small** Random Forest model is **not included** in this repository because it exceeds GitHub's 100 MB file size limit.

Download it from Zenodo:

https://doi.org/10.5281/zenodo.21620365

Then place it in:

```text
haloClassifier_pkg/models/small/
```

The directory should contain:

```text
random_forest_top265_model_small.pkl
top265_features.csv
```

HaloClassifier will then run normally.

---

## Quick start

```bash
conda create -n haloclassifier python=3.12
conda activate haloclassifier

git clone https://github.com/Lorena0105/HaloClassifier.git
cd HaloClassifier

pip install -e .

halo-generate-features path/to/fasta_folder/ [-o contigs_features.csv] [-t THREADS]
halo-classify contigs_features.csv [-o predictions_contigs.tsv] [-p THRESHOLD]
```

See the **Installation** section below for full details.


## Installation

### Requirements

- Python ≥ 3.10
- Prodigal (tested with v2.6.3)

### Python dependencies

The following packages will be installed automatically.
The versions shown below correspond to those used during development.

- pandas (2.3.3)
- scikit-learn (1.7.2)
- joblib (1.5.2)
- tqdm (4.67.1)

To provide a more reproducible installation, HaloClassifier restricts the versions of some dependencies to ranges that have been tested. This avoids potential compatibility issues caused by major changes in external libraries.

The installation therefore allows:
- `pandas >=2.3.3, <3.0`
- `scikit-learn >=1.7.2, <1.8`
- `joblib >=1.5.2, <2.0`

Installing different versions may also work, but these configurations have not been tested with the current release of HaloClassifier.

### Installation

```bash
conda create -n haloclassifier python=3.12
conda activate haloclassifier

git clone https://github.com/Lorena0105/HaloClassifier.git
cd HaloClassifier

pip install -e .
```

## Verify installation

```bash
haloclassifier
```

To test the example dataset:

```bash
halo-generate-features test/
halo-classify contigs_features.csv
```

## Pipeline Overview

## 1) Feature Generation

In this stage, input sequences are processed to generate a feature matrix in CSV format.

### What happens during feature generation?

- The program **scans the input folder** and processes only valid nucleotide sequence files (`.fasta`, `.fa`, `.fna`, `.fas`).
- Multi-contig FASTA files are supported.
- Contigs **>105 kb** are automatically fragmented into pieces between 40 kb and 105 kb, ensuring that all resulting fragments fall within the **extralarge** model range.
- Contigs **<1 kb** are not excluded, but are **flagged** and reported separately in `short_contigs_list.csv`, since they fall outside the training range.
- Each contig (or fragment) is first assigned to its corresponding size-specific model. HaloClassifier then computes **only the features required by that model**, avoiding unnecessary calculations.

The computed features may include:

#### **Non-coding features**
- GC content.  
- AT/GC skew.  
- Canonical k-mers (k = 2–7).

#### **Coding features** (ORF prediction via **Prodigal**)
- Complete and partial gene counts.
- Average complete gene length.
- Coding density.
- Codon usage (relative frequency).
- RSCU.
- ENC.

#### Parallel processing [-t]
The feature generation stage uses **all available CPU threads by default** for parallel processing.  
You can limit the number of threads by specifying the `-t` option.

Example usage:

```bash
halo-generate-features path/to/fasta_folder/ [-o contigs_features.csv] [-t THREADS]
```

This produces:
- `contigs_features.csv` (or a user-defined name) → main feature table for all contigs/fragments (CSV).
- `short_contigs_list.csv` (only if applicable) → list of sequences <1 kb.
- `haloClassifier_input_warnings.txt` (only if applicable) → files ignored and the reason.


## 2) Classification

The generated feature table is then classified using **pre-trained size-specific Random Forest models**.

- The pipeline maps each contig (or fragment) to a model based on its length:

| Model        | Range (bp)      |
|--------------|-----------------|
| Small        | 1,000–5,000     |
| Medium       | 5,001–10,000    |
| Large        | 10,001–40,000   |
| Extralarge   | 40,001–100,000  |

- Contigs <1 kb → classified with *small* model but flagged (`below_training_range = True`) in the final report.  
- Contigs 100–105 kb → treated as *extralarge* (no fragmentation).  
- Contigs >105 kb → automatically fragmented into pieces between 40 kb and 105 kb before classification. Fragment probabilities are averaged for the final prediction.

#### Probability Threshold [-p]
Users can optionally filter contigs based on a minimum predicted probability for either chromosome or plasmid.
A complementary file will be automatically generated → predictions_contigs_threshold_0.60.tsv (default for p = 0.6).

Example usage:

```bash
halo-classify contigs_features.csv [-o predictions_contigs.tsv] [-p THRESHOLD]
```

This produces:
- `predictions_contigs.tsv` (or a user-defined name) → final report file (TSV).
- `subfragments/` → predictions for each fragment generated from contigs >105 kb. The final report contains the mean probability for each contig, whereas this folder stores fragment-level predictions.
- `predictions_contigs_threshold_<THRESHOLD>.tsv` (e.g. `predictions_contigs_threshold_0.60.tsv`) → contains only contigs with a predicted probability ≥ the specified threshold (e.g. 0.6).

**Final report file**/**Threshold file** includes:

- `Assembly_ID`
- `Contig_ID` (`NA` for single-contig FASTA files)
- `Fragmentation` (`True` if the contig was fragmented before classification)
- `Length`
- `Prediction` (`Chromosome`/`Plasmid`)
- `Prob_Chromosome`
- `Prob_Plasmid`
- `Model_Used`
- `Below_Training_Range` (`True` if the contig is <1 kb)

## Citation

If you use HaloClassifier in your research, please cite the corresponding publication.

Citation information will be added upon publication.
