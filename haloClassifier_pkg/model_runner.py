#!/usr/bin/env python3
import pandas as pd
import joblib
import os
import multiprocessing
from tqdm import tqdm
from importlib.resources import files

# -----------------------------------------------------
# Model size configuration
# -----------------------------------------------------
MODELS = {
    "small": {"min_bp": 1000, "max_bp": 5000},
    "medium": {"min_bp": 5001, "max_bp": 10000},
    "large": {"min_bp": 10001, "max_bp": 40000},
    "extralarge": {"min_bp": 40001, "max_bp": 100000}
}

def get_top_n(size_class):
    return {"small": 265, "medium": 350, "large": 355, "extralarge": 595}[size_class]

def get_model_paths(size_class):
    model_dir = files(f"haloClassifier_pkg.models.{size_class}")
    model_file = model_dir / f"random_forest_top{get_top_n(size_class)}_model_{size_class}.pkl"
    features_file = model_dir / f"top{get_top_n(size_class)}_features.csv"
    return model_file, features_file

def choose_model(length):
    if length < 5000:
        return "small"
    elif length <= 10000:
        return "medium"
    elif length <= 40000:
        return "large"
    else:
        return "extralarge"

# -----------------------------------------------------
# MAIN FUNCTION
# -----------------------------------------------------
def main(features_csv, short_contigs_csv=None,
         output_tsv="predictions_contigs.tsv",
         subfrag_dir="subfragments",
         n_jobs=1,
         probability_threshold=None):

    if n_jobs == -1:
        n_jobs = multiprocessing.cpu_count()

    df = pd.read_csv(features_csv)
    if "length" not in df.columns:
        raise ValueError("Feature CSV must contain a 'length' column.")

    df["model_group"] = df["length"].apply(choose_model)

    # -----------------------------------------------------
    # Load short contigs list (<1kb)
    # -----------------------------------------------------
    short_set = set()
    if short_contigs_csv is None:
        possible_path = os.path.join(os.path.dirname(features_csv),
                                     "short_contigs_list.csv")
        if os.path.exists(possible_path):
            short_contigs_csv = possible_path
    if short_contigs_csv and os.path.exists(short_contigs_csv):
        short_df = pd.read_csv(short_contigs_csv)
        short_set = set(short_df["contig_id"])

    final_predictions = []

    # -----------------------------------------------------
    # Process contigs
    # -----------------------------------------------------
    for (assembly_id, contig_id), group in tqdm(
            df.groupby(["assembly_id", "contig_id"]),
            desc="Processing contigs"):

        subfrag_list = []
        probs_chr = []
        probs_plasm = []
        models_used = []

        for _, row in group.iterrows():
            frag_id = row["fragment_id"]
            length = row["length"]
            model_key = row["model_group"]

            # Load model if needed
            if model_key not in globals():
                model_file, features_file = get_model_paths(model_key)
                clf = joblib.load(model_file)
                feature_list = pd.read_csv(features_file, header=None)[0].tolist()
                globals()[model_key] = (clf, feature_list)
            else:
                clf, feature_list = globals()[model_key]

            X_dict = {feat: row.get(feat, 0) for feat in feature_list}
            X = pd.DataFrame([X_dict])

            proba = clf.predict_proba(X)[0]
            label = clf.predict(X)[0]

            # -------------------------------
            # Parse fragment_id → clean + start/end
            # -------------------------------
            parts = frag_id.split("_")
            if len(parts) >= 4 and parts[-3].startswith("frag"):
                clean_fragment_id = "_".join(parts[:-2])
                try:
                    start = int(parts[-2])
                    end = int(parts[-1])
                except ValueError:
                    start, end = 0, length
            else:
                clean_fragment_id = frag_id
                start, end = 0, length

            subfrag_list.append({
                "contig_id": contig_id,
                "fragment_id": clean_fragment_id,
                "start": start,
                "end": end,
                "length": length,
                "prob_chromosome": proba[0],
                "prob_plasmid": proba[1],
                "model_used": model_key,
                "prediction": label if max(proba) > 0.5 else "Ambiguous"
            })

            probs_chr.append(proba[0])
            probs_plasm.append(proba[1])
            models_used.append(model_key)

        # Save subfragment file only if fragmentation occurred
        if len(subfrag_list) > 1:
            os.makedirs(subfrag_dir, exist_ok=True)
            out_path = os.path.join(subfrag_dir,
                                    f"{assembly_id}_{contig_id}_subfragments.tsv")
            cols = [
                "contig_id", "fragment_id", "start", "end",
                "length", "prob_chromosome", "prob_plasmid",
                "model_used", "prediction"
            ]
            df_sub = pd.DataFrame(subfrag_list)[cols].sort_values("start")
            df_sub.to_csv(out_path, sep="\t", index=False)

        # Average probabilities
        avg_chr = sum(probs_chr) / len(probs_chr)
        avg_plasm = sum(probs_plasm) / len(probs_plasm)
        model_final = max(set(models_used), key=models_used.count)

        final_pred = {
            "Assembly_ID": assembly_id,
            "Contig_ID": contig_id,
            "Fragmentation": len(subfrag_list) > 1,
            "Length": sum(sf["length"] for sf in subfrag_list),
            "Prediction": "Chromosome" if avg_chr > 0.5 else
                          "Plasmid" if avg_plasm > 0.5 else
                          "Ambiguous",
            "Prob_Chromosome": avg_chr,
            "Prob_Plasmid": avg_plasm,
            "Model_Used": model_final,
            "Below_Training_Range": contig_id in short_set
        }

        if final_pred["Assembly_ID"] == final_pred["Contig_ID"]:
            final_pred["Contig_ID"] = "NA"

        final_predictions.append(final_pred)

    # Save full predictions
    df_final = pd.DataFrame(final_predictions)
    df_final.to_csv(output_tsv, sep="\t", index=False)
    print(f"\n✅ Predictions saved to: {output_tsv}")
    if os.path.exists(subfrag_dir):
    	print(f"Subfragment files stored in: {subfrag_dir}/")

    # -----------------------------------------------------
    # Apply probability threshold filter if requested
    # -----------------------------------------------------
    if probability_threshold is not None:
        df_filtered = df_final[df_final[["Prob_Chromosome","Prob_Plasmid"]].max(axis=1) >= probability_threshold]
        base, ext = os.path.splitext(output_tsv)
        threshold_file = f"{base}_threshold_{probability_threshold:.2f}.tsv"
        df_filtered.to_csv(threshold_file, sep="\t", index=False)
        print(f"✅ Filtered predictions saved to: {threshold_file}")


# -----------------------------------------------------
# CLI mode
# -----------------------------------------------------
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python model_runner.py <features_csv> [<short_contigs_csv>] [-p probability_threshold]")
        sys.exit(1)

    features_csv = sys.argv[1]
    short_csv = sys.argv[2] if len(sys.argv) > 2 else None
    main(features_csv, short_csv)
