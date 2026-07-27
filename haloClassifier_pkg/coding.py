#!/usr/bin/env python3
import os
import subprocess
import glob
import csv
import re
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed

# Codon table by amino acid
AMINO_CODONS = {
    'F': ['TTT','TTC'], 'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
    'I': ['ATT','ATC','ATA'], 'M': ['ATG'], 'V': ['GTT','GTC','GTA','GTG'],
    'S': ['TCT','TCC','TCA','TCG','AGT','AGC'], 'P': ['CCT','CCC','CCA','CCG'],
    'T': ['ACT','ACC','ACA','ACG'], 'A': ['GCT','GCC','GCA','GCG'],
    'Y': ['TAT','TAC'], 'H': ['CAT','CAC'], 'Q': ['CAA','CAG'], 'N': ['AAT','AAC'],
    'K': ['AAA','AAG'], 'D': ['GAT','GAC'], 'E': ['GAA','GAG'],
    'C': ['TGT','TGC'], 'W': ['TGG'], 'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
    'G': ['GGT','GGC','GGA','GGG'],
}
STOP_CODONS = {'TAA','TAG','TGA'}

# --------------------
# HELPER FUNCTIONS
# --------------------
def run_prodigal(fasta_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}.ffn")
    if os.path.isfile(output_file):
        return output_file
    cmd = ["prodigal", "-i", fasta_file, "-d", output_file, "-q", "-p", "meta"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return output_file

def read_fasta_length(fasta_file):
    return sum(len(line.strip()) for line in open(fasta_file) if not line.startswith(">"))

def read_ffn_gene_stats(ffn_file):
    total_genes, complete_lengths, intervals = 0, [], []
    current_seq, is_complete, start, end = "", False, None, None
    with open(ffn_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    if is_complete: complete_lengths.append(len(current_seq))
                    intervals.append((start, end))
                total_genes += 1
                current_seq = ""
                m_coord = re.search(r"#\s*(\d+)\s*#\s*(\d+)\s*#", line)
                if m_coord: start, end = int(m_coord.group(1)), int(m_coord.group(2))
                m_part = re.search(r"partial=(\d{2})", line)
                is_complete = (m_part.group(1) == "00") if m_part else False
            else:
                current_seq += line
        if current_seq:
            if is_complete: complete_lengths.append(len(current_seq))
            intervals.append((start, end))
    return total_genes, complete_lengths, intervals

def merge_intervals(intervals):
    if not intervals: return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for cs, ce in intervals[1:]:
        ps, pe = merged[-1]
        if cs <= pe: merged[-1] = (ps, max(pe, ce))
        else: merged.append((cs, ce))
    return merged

def compute_coding_density(intervals, fasta_len):
    merged = merge_intervals(intervals)
    coding_bases = sum(e - s + 1 for s, e in merged)
    return coding_bases, coding_bases / fasta_len if fasta_len else 0

# --------------------
# CODON PROCESSING
# --------------------
def clean_sequence(seq):
    return "".join(b for b in seq.upper() if b in ("A","C","G","T"))

def process_cds_sequence(seq, is_partial_01):
    seq = clean_sequence(seq)
    if is_partial_01 and len(seq) % 3: seq = seq[:-(len(seq) % 3)]
    return seq

def codon_usage_all_codons(ffn_file):
    bases = ["A","C","G","T"]
    all_codons = ["".join(p) for p in product(bases, repeat=3)]
    codon_counts = {c: 0 for c in all_codons}
    total_codons, current_seq, is_partial_01 = 0, [], False
    with open(ffn_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    seq = process_cds_sequence("".join(current_seq), is_partial_01)
                    for i in range(0, len(seq), 3):
                        cod = seq[i:i+3]
                        if len(cod) == 3:
                            codon_counts[cod] += 1
                            total_codons += 1
                current_seq = []
                is_partial_01 = "partial=01" in line
            else:
                current_seq.append(line)
        if current_seq:
            seq = process_cds_sequence("".join(current_seq), is_partial_01)
            for i in range(0, len(seq), 3):
                cod = seq[i:i+3]
                if len(cod) == 3:
                    codon_counts[cod] += 1
                    total_codons += 1
    codon_freq = {c: codon_counts[c]/total_codons if total_codons else 0 for c in all_codons}
    return codon_counts, codon_freq

def calculate_rscu(codon_counts):
    rscu = {}
    for aa, codons in AMINO_CODONS.items():
        total = sum(codon_counts.get(c, 0) for c in codons)
        n = len(codons)
        for c in codons: rscu[c] = codon_counts.get(c, 0) * n / total if total > 0 else 0.0
    return rscu

def calculate_enc(codon_counts):
    deg2, deg3, deg4, deg6 = [], [], [], []
    for aa, codons in AMINO_CODONS.items():
        k = len([c for c in codons if c not in STOP_CODONS])
        if k == 2: deg2.append(aa)
        elif k == 3: deg3.append(aa)
        elif k == 4: deg4.append(aa)
        elif k == 6: deg6.append(aa)
    def F_k(counts, codons):
        n_list = [counts.get(c, 0) for c in codons if c not in STOP_CODONS]
        n = sum(n_list)
        if n <= 1 or len(n_list) == 1: return 1.0
        sum_sq = sum(x**2 for x in n_list)
        return (sum_sq - n) / (n * (n - 1))
    def avg_F(aa_list):
        if not aa_list: return None
        return sum(F_k(codon_counts, AMINO_CODONS[aa]) for aa in aa_list)/len(aa_list)
    F2, F3, F4, F6 = avg_F(deg2), avg_F(deg3), avg_F(deg4), avg_F(deg6)
    ENC = 2.0
    if F2: ENC += 9.0 / F2
    if F3: ENC += 1.0 / F3
    if F4: ENC += 5.0 / F4
    if F6: ENC += 3.0 / F6
    return max(20.0, min(61.0, ENC))

# --------------------
# FRAGMENT CONTIGS > 105kb
# --------------------
def fragment_contig(fasta_file, max_len=105000, supplement_dir="subfragments"):
    os.makedirs(supplement_dir, exist_ok=True)
    fragments = []
    seq_id = None
    seq_lines = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    seq = "".join(seq_lines)
                    for i in range(0, len(seq), max_len):
                        frag_seq = seq[i:i+max_len]
                        frag_name = f"{seq_id}_frag{i//max_len+1}.fasta"
                        frag_path = os.path.join(supplement_dir, frag_name)
                        with open(frag_path, "w") as ff:
                            ff.write(f">{seq_id}_frag{i//max_len+1}\n{frag_seq}\n")
                        fragments.append(frag_path)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            seq = "".join(seq_lines)
            for i in range(0, len(seq), max_len):
                frag_seq = seq[i:i+max_len]
                frag_name = f"{seq_id}_frag{i//max_len+1}.fasta"
                frag_path = os.path.join(supplement_dir, frag_name)
                with open(frag_path, "w") as ff:
                    ff.write(f">{seq_id}_frag{i//max_len+1}\n{frag_seq}\n")
                fragments.append(frag_path)
    return fragments

# --------------------
# PROCESS EACH CONTIG/SUBFRAGMENT
# --------------------
def process_contig_or_fragment(fasta_file, ffn_folder, all_codons):
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    fasta_len = read_fasta_length(fasta_file)
    ffn_file = run_prodigal(fasta_file, ffn_folder)

    total_genes, complete_lengths, intervals = read_ffn_gene_stats(ffn_file)
    num_complete = len(complete_lengths)
    avg_len_complete = sum(complete_lengths)/num_complete if num_complete else 0
    _, coding_density = compute_coding_density(intervals, fasta_len)

    codon_counts, codon_freq = codon_usage_all_codons(ffn_file)
    rscu = calculate_rscu(codon_counts)
    enc = calculate_enc(codon_counts)

    row = {"fragment_id": base_name, "total_genes": total_genes, "complete_genes": num_complete,
           "avg_length_complete": avg_len_complete, "coding_density": coding_density, "ENC": enc}
    for c in all_codons:
        row[f"codon_{c}"] = codon_freq.get(c, 0.0)
        row[f"RSCU_{c}"] = rscu.get(c, 0.0)
    return row

# --------------------
# MAIN
# --------------------
def main(fasta_folder, ffn_folder="prodigal_ffn", output_csv="genes_stats.csv"):
    fasta_files = glob.glob(os.path.join(fasta_folder, "*.fasta"))
    bases = ["A","C","G","T"]
    all_codons = sorted(["".join(p) for p in product(bases, repeat=3)])
    results = []

    print(f"Processing {len(fasta_files)} files...\n")
    with ProcessPoolExecutor() as executor:
        futures = {}
        for f in fasta_files:
            fragments = fragment_contig(f)
            for frag in fragments:
                futures[executor.submit(process_contig_or_fragment, frag, ffn_folder, all_codons)] = frag
        for fut in as_completed(futures):
            results.append(fut.result())
            print(f"✓ Completed: {futures[fut]}")

    fieldnames = (["fragment_id","total_genes","complete_genes","avg_length_complete","coding_density","ENC"]
                  + [f"codon_{c}" for c in all_codons] + [f"RSCU_{c}" for c in all_codons])
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for r in results: writer.writerow(r)
    print("\nCSV saved to:", output_csv)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python coding.py <fasta_folder>")
        sys.exit(1)
    main(sys.argv[1])
