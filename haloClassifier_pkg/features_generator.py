#!/usr/bin/env python3
import os
import csv
import math
import uuid
import shutil
import time
import subprocess
import sys
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from importlib.resources import files

from . import input_reader
from . import noncoding
from . import coding

MAX_FRAGMENT_SIZE = 105000  # real fragmentation threshold

# -----------------------------------------------------------
# Check external dependencies
# -----------------------------------------------------------

def check_prodigal():
    """Exit if Prodigal is not installed or not available in PATH."""

    try:
        subprocess.run(
            ["prodigal", "-h"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

    except (FileNotFoundError, subprocess.CalledProcessError):

        sys.exit(
            "\nERROR: Prodigal is not installed or is not available in your PATH.\n\n"
            "Please install Prodigal before running HaloClassifier.\n"
        )


# Verify external programs before doing anything
check_prodigal()

# -----------------------------------------------------------
# Load feature lists used by each trained model
# -----------------------------------------------------------

def load_required_features():

    feature_sets = {}

    model_info = {
        "small": 265,
        "medium": 350,
        "large": 355,
        "extralarge": 595
    }

    for model, n in model_info.items():

        feature_file = (
            files(f"haloClassifier_pkg.models.{model}")
            / f"top{n}_features.csv"
        )

        with open(feature_file) as f:
            feature_sets[model] = set(
                line.strip() for line in f if line.strip()
            )

    return feature_sets


REQUIRED_FEATURES = load_required_features()

# -----------------------------------------------------------
# Split features only once
# -----------------------------------------------------------

REQUIRED_NONCODING = {}
REQUIRED_CODING = {}
REQUIRED_CODONS = {}
REQUIRED_RSCU = {}

for model, feats in REQUIRED_FEATURES.items():

    REQUIRED_NONCODING[model] = {
        f
        for f in feats
        if (
            f in {"GC_content", "AT_skew", "GC_skew"}
            or "-mer_" in f
        )
    }

    REQUIRED_CODING[model] = {
        f
        for f in feats
        if (
            f.startswith("codon_")
            or f.startswith("RSCU_")
            or f in {
                "total_genes",
                "complete_genes",
                "avg_length_complete",
                "coding_density",
                "ENC",
            }
        )
    }
    
    REQUIRED_CODONS[model] = [
        f[6:]
        for f in feats
        if f.startswith("codon_")
    ]

    REQUIRED_RSCU[model] = [
        f[5:]
        for f in feats
        if f.startswith("RSCU_")
    ]


def choose_model(length):

    if length <= 5000:
        return "small"
    elif length <= 10000:
        return "medium"
    elif length <= 40000:
        return "large"
    else:
        return "extralarge"

# -----------------------------------------------------------
# Fragmentation of contigs > 105 kb
# -----------------------------------------------------------
def fragment_contig(sequence, contig_id):
    """Split a sequence into fragments if it exceeds MAX_FRAGMENT_SIZE."""
    L = len(sequence)

    # No fragmentation needed
    if L <= MAX_FRAGMENT_SIZE:
        return [("full", sequence, 0, L)]

    fragments = []
    n_fragments = math.ceil(L / MAX_FRAGMENT_SIZE)
    frag_size = math.ceil(L / n_fragments)

    start = 0
    for i in range(n_fragments):
        end = start + frag_size
        if i == n_fragments - 1:
            end = L

        frag_id = f"{contig_id}_frag{i+1}_{start}_{end}"
        frag_seq = sequence[start:end]

        fragments.append((frag_id, sequence, start, end))
        start = end

    return fragments

# -----------------------------------------------------------
# Compute features for one fragment
# -----------------------------------------------------------
def compute_fragment_features(fragment_id, seq, start, end):
    """Compute coding and non-coding features for a fragment."""

    frag_seq = seq[start:end]
    length = len(frag_seq)

    model = choose_model(length)

    required_nc = REQUIRED_NONCODING[model]
    required_c = REQUIRED_CODING[model]
    required_codons = REQUIRED_CODONS[model]
    required_rscu = REQUIRED_RSCU[model]

    # -------------------------------------------------------
    # Non-coding features
    # -------------------------------------------------------

    res_nc = noncoding.process_sequence(
        frag_seq,
        required_nc,
	model
    )

    # -------------------------------------------------------
    # Temporary FASTA for Prodigal
    # -------------------------------------------------------

    unique_id = uuid.uuid4().hex
    temp_fasta = f"temp_{fragment_id}_{unique_id}.fasta"

    with open(temp_fasta, "w") as f:
        f.write(f">{fragment_id}\n{frag_seq}\n")

    frag_ffn = coding.run_prodigal(
        temp_fasta,
        "prodigal_fragments"
    )

    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)

    # -------------------------------------------------------
    # Coding features
    # -------------------------------------------------------

    total_genes, complete_lengths, intervals = (
        coding.read_ffn_gene_stats(frag_ffn)
    )

    codon_counts, codon_freq = coding.codon_usage_all_codons(
        frag_ffn
    )

    rscu = coding.calculate_rscu(
        codon_counts
    )

    enc = 0.0

    if "ENC" in required_c:
        enc = coding.calculate_enc(codon_counts)

    if os.path.exists(frag_ffn):
        os.remove(frag_ffn)

    num_complete = len(complete_lengths)

    avg_len_complete = (
        sum(complete_lengths) / num_complete
        if num_complete else 0
    )

    _, coding_density = coding.compute_coding_density(
        intervals,
        length
    )

    res_c = {}

    if "total_genes" in required_c:
        res_c["total_genes"] = total_genes

    if "complete_genes" in required_c:
        res_c["complete_genes"] = num_complete

    if "avg_length_complete" in required_c:
        res_c["avg_length_complete"] = avg_len_complete

    if "coding_density" in required_c:
        res_c["coding_density"] = coding_density

    if "ENC" in required_c:
        res_c["ENC"] = enc

    # -------------------------------------------------------
    # Only requested codon features
    # -------------------------------------------------------

    for cod in required_codons:
        res_c[f"codon_{cod}"] = codon_freq.get(cod, 0.0)

    for cod in required_rscu:
        res_c[f"RSCU_{cod}"] = rscu.get(cod, 0.0)

    return {
        "fragment_id": fragment_id,
        "length": length,
        **res_nc,
        **res_c
    }

# -----------------------------------------------------------
# Process a contig (fragmented if >105 kb)
# -----------------------------------------------------------
def process_contig(record):
    """Process a contig and compute features for all its fragments."""
    assembly_id = record["assembly_id"]
    contig_id = record["contig_id"]
    seq = record["sequence"]

    fragments = fragment_contig(seq, contig_id)

    results = []
    for frag_id, full_seq, start, end in fragments:
        frag_features = compute_fragment_features(
            frag_id, full_seq, start, end
        )
        frag_features["contig_id"] = contig_id
        frag_features["assembly_id"] = assembly_id
        results.append(frag_features)

    return results


# -----------------------------------------------------------
# MAIN
# -----------------------------------------------------------
def calc_features_for_contigs(folder_path, output_csv="contigs_features.csv",
                              short_contigs_csv="short_contigs_list.csv",
                              max_workers=None):

    records, ignored = input_reader.load_contigs(folder_path)
    input_reader.write_warning_file(ignored)

    # -------------------------------------------------------
    # Determine base directory (input can be a folder or a FASTA file)
    # -------------------------------------------------------

    if os.path.isdir(folder_path):
    	base_dir = folder_path
    else:
    	base_dir = os.path.dirname(folder_path)

    # Create Prodigal folder inside the base directory
    prodigal_dir_input = os.path.join(base_dir, "prodigal_fragments")
    os.makedirs(prodigal_dir_input, exist_ok=True)

    results = []
    short_contigs = []

    print("\nProcessing contigs (parallel mode)...")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_contig, rec): rec
            for rec in records
        }

        for fut in tqdm(as_completed(futures), total=len(futures)):
            frag_list = fut.result()
            for frag in frag_list:
                results.append(frag)
                if frag["length"] < 1000 and frag["fragment_id"] == "full":
                    short_contigs.append({
                        "assembly_id": frag["assembly_id"],
                        "contig_id": frag["contig_id"],
                        "length": frag["length"]
                    })

    # Save features CSV
    all_columns = sorted({k for row in results for k in row})
    ordered_cols = ["assembly_id", "contig_id", "fragment_id", "length"]
    all_columns = ordered_cols + [c for c in all_columns if c not in ordered_cols]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_columns)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\n📕 Feature file generated: {output_csv}")

    # Save short contigs (<1 kb)
    if short_contigs:
        with open(short_contigs_csv, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=["assembly_id", "contig_id", "length"])
            writer.writeheader()
            for row in short_contigs:
                writer.writerow(row)
        print(f"⚠️ Short contigs (<1 kb) saved in {short_contigs_csv}")

    # -------------------------------------------------------
    # Remove Prodigal folders
    # -------------------------------------------------------
    time.sleep(0.5)  # allow Prodigal to finish releasing file handles

    # Prodigal folder inside input FASTA directory
    if os.path.exists(prodigal_dir_input):
        try:
            shutil.rmtree(prodigal_dir_input)
        except Exception as e:
            print(f"\n⚠️ Could not delete {prodigal_dir_input}: {e}")

    # Prodigal folder inside current working directory
    prodigal_dir_cwd = os.path.join(os.getcwd(), "prodigal_fragments")
    if os.path.exists(prodigal_dir_cwd):
        try:
            shutil.rmtree(prodigal_dir_cwd)
        except Exception as e:
            print(f"\n⚠️ Could not delete {prodigal_dir_cwd}: {e}")
