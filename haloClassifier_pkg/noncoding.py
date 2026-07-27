#!/usr/bin/env python3

import pickle
from importlib.resources import files

# -------------------------------------------------------
# Load lookup tables (loaded once when importing module)
# -------------------------------------------------------

with open(files("haloClassifier_pkg.models") / "kmer_lookup.pkl", "rb") as f:
    KMER_LOOKUPS = pickle.load(f)


# -------------------------------------------------------
# Basic sequence statistics
# -------------------------------------------------------

def gc_content(seq):

    g = seq.count("G")
    c = seq.count("C")

    return (g + c) * 100.0 / len(seq)


def at_gc_skew(seq):

    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")

    at = (a - t) / (a + t) if (a + t) else 0.0
    gc = (g - c) / (g + c) if (g + c) else 0.0

    return at, gc


# -------------------------------------------------------
# Selected k-mer frequencies
# -------------------------------------------------------

def selected_kmer_frequencies(seq, model):

    model_tables = KMER_LOOKUPS[model]

    tables = {}

    totals = {}

    counts = {}

    for k, data in model_tables.items():

        tables[k] = data["lookup"]

        totals[k] = 0

        for feature in data["features"]:
            counts[feature] = 0

    seq_len = len(seq)

    for i in range(seq_len):

        if 2 in tables and i + 2 <= seq_len:

            kmer = seq[i:i+2]

            if "N" not in kmer:

                totals[2] += 1

                feat = tables[2].get(kmer)

                if feat is not None:
                    counts[feat] += 1


        if 3 in tables and i + 3 <= seq_len:

            kmer = seq[i:i+3]

            if "N" not in kmer:

                totals[3] += 1

                feat = tables[3].get(kmer)

                if feat is not None:
                    counts[feat] += 1


        if 4 in tables and i + 4 <= seq_len:

            kmer = seq[i:i+4]

            if "N" not in kmer:

                totals[4] += 1

                feat = tables[4].get(kmer)

                if feat is not None:
                    counts[feat] += 1


        if 5 in tables and i + 5 <= seq_len:

            kmer = seq[i:i+5]

            if "N" not in kmer:

                totals[5] += 1

                feat = tables[5].get(kmer)

                if feat is not None:
                    counts[feat] += 1


        if 6 in tables and i + 6 <= seq_len:

            kmer = seq[i:i+6]

            if "N" not in kmer:

                totals[6] += 1

                feat = tables[6].get(kmer)

                if feat is not None:
                    counts[feat] += 1


        if 7 in tables and i + 7 <= seq_len:

            kmer = seq[i:i+7]

            if "N" not in kmer:

                totals[7] += 1

                feat = tables[7].get(kmer)

                if feat is not None:
                    counts[feat] += 1


    results = {}

    for k, data in model_tables.items():

        total = totals[k]

        if total == 0:

            for feature in data["features"]:
                results[feature] = 0.0

            continue

        inv_total = 1.0 / total

        for feature in data["features"]:
            results[feature] = counts[feature] * inv_total

    return results


# -------------------------------------------------------
# Main function
# -------------------------------------------------------

def process_sequence(seq, required_features, model):

    results = {}

    if "GC_content" in required_features:
        results["GC_content"] = gc_content(seq)

    if (
        "AT_skew" in required_features
        or "GC_skew" in required_features
    ):

        at, gc = at_gc_skew(seq)

        if "AT_skew" in required_features:
            results["AT_skew"] = at

        if "GC_skew" in required_features:
            results["GC_skew"] = gc

    results.update(
        selected_kmer_frequencies(
            seq,
            model
        )
    )

    return results
