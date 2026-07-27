import argparse
from . import model_runner, features_generator

# ---------------------------
#  GENERAL COMMAND (MENU)
# ---------------------------
def menu():

    parser = argparse.ArgumentParser(
        prog="haloClassifier",
        description=(
            "haloClassifier: Classification of Haloarchaea contigs as chromosomal or plasmidic\n"
            "Version 1.0.1"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.epilog = (
        "Available commands:\n"
        "  halo-generate-features   Generate a CSV feature table from FASTA files\n"
        "  halo-classify            Classify contigs using pre-trained models\n\n"
        "Examples:\n"
        "  halo-generate-features fasta_folder/ [-o contigs_features.csv] [-t THREADS]\n"
        "      # Default output file: contigs_features.csv\n"
        "      # Optional: use -t to select number of threads (default = all cores)\n\n"
        "  halo-classify contigs_features.csv [-o predictions_contigs.tsv] [-p THRESHOLD]\n"
        "      # Default output file: predictions_contigs.tsv\n"
        "      # Optional: filter contigs with minimum probability threshold (default = 0.5)\n"
    )

    parser.print_help()


# ---------------------------
#  ENTRY POINT 1 — Feature generation
# ---------------------------
def run_generate():
    parser = argparse.ArgumentParser(
        description="Generate a CSV feature table from FASTA input files"
    )
    parser.add_argument("fasta_folder", help="Folder containing FASTA files")
    parser.add_argument(
        "-o", "--output",
        default="contigs_features.csv",
        help="Name of the output CSV file (default: contigs_features.csv)"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=None,
        help="Number of CPU threads to use (default = all available cores)"
    )
    args = parser.parse_args()

    features_generator.calc_features_for_contigs(
        folder_path=args.fasta_folder,
        output_csv=args.output,
        max_workers=args.threads
    )


# ---------------------------
#  ENTRY POINT 2 — Classification
# ---------------------------
def run_classify():
    parser = argparse.ArgumentParser(
        description="Classify contigs using the trained Random Forest models"
    )
    parser.add_argument("features_csv", help="CSV file containing contig feature values")
    parser.add_argument(
        "-o", "--output",
        default="predictions_contigs.tsv",
        help="Name of the output TSV file (default: predictions_contigs.tsv)"
    )
    parser.add_argument(
        "-p", "--threshold", type=float, default=None,
        help="Optional probability threshold to filter final contigs (0.0-1.0, default = 0.5)"
    )
    args = parser.parse_args()

    model_runner.main(
        features_csv=args.features_csv,
        output_tsv=args.output,
        probability_threshold=args.threshold
    )


# ---------------------------
#  LEGACY COMPATIBILITY
# ---------------------------
def main():
    menu()
