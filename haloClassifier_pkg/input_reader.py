import os

# Valid nucleotide FASTA extensions
VALID_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas")

# Valid IUPAC nucleotide codes
VALID_BASES = set("ATGCNRYSWKMBVDH")

# Minimum length used in training
TRAINING_MIN_LENGTH = 1000

def scan_input_folder(path):
    """
    Accept either

      - a directory containing FASTA files
      - a single FASTA file

    Returns
        valid_files
        ignored
    """

    valid_files = []
    ignored = []

    # ----------------------------------------
    # Single FASTA file
    # ----------------------------------------

    if os.path.isfile(path):

        filename = os.path.basename(path)

        if filename.lower().endswith(VALID_EXTENSIONS):
            valid_files.append(path)
        else:
            ignored.append((filename, "unsupported extension"))

        return valid_files, ignored

    # ----------------------------------------
    # Directory
    # ----------------------------------------

    if not os.path.isdir(path):
        raise ValueError(f"{path} is not a valid file or directory")

    for filename in os.listdir(path):

        file_path = os.path.join(path, filename)

        if not os.path.isfile(file_path):
            continue

        if filename.lower().endswith(VALID_EXTENSIONS):
            valid_files.append(file_path)
        else:
            ignored.append((filename, "unsupported extension"))

    return valid_files, ignored


def read_multifasta(file_path):
    """
    Read a FASTA or multi-FASTA file and return a list of
    (header_id, sequence) tuples.
    """
    contigs = []
    seq_id = None
    seq_lines = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if seq_id:
                    contigs.append((seq_id, "".join(seq_lines)))
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.upper())

    if seq_id:
        contigs.append((seq_id, "".join(seq_lines)))

    return contigs


def sequence_is_valid(seq):
    """
    Returns True if the sequence contains only valid IUPAC nucleotide characters.
    """
    return all(base in VALID_BASES for base in seq)


def load_contigs(folder):
    """
    Load all contigs from valid FASTA files.
    
    Additional rule:
    - If any contig inside a file contains invalid characters → ignore the whole file.
    
    Returns:
        records → list of contig dictionaries
        ignored → list of ignored files with reasons
    """

    valid_files, ignored = scan_input_folder(folder)
    records = []

    for file_path in valid_files:
        filename = os.path.basename(file_path)
        assembly_id = os.path.splitext(filename)[0]

        contigs = read_multifasta(file_path)

        # Check validity of ALL contigs. If any is invalid → ignore file.
        invalid_count = 0
        for contig_id, seq in contigs:
            if not sequence_is_valid(seq):
                invalid_count += 1

        if invalid_count > 0:
            if len(contigs) == 1:
                reason = "invalid characters in sequence"
            else:
                reason = "invalid characters in one or more contigs"

            ignored.append((filename, reason))
            continue  # skip entire file

        # If file is valid → add all its contigs
        for contig_id, seq in contigs:
            seq_len = len(seq)
            records.append({
                "assembly_id": assembly_id,
                "contig_id": contig_id,
                "sequence": seq,
                "length": seq_len,
                "below_training_range": seq_len < TRAINING_MIN_LENGTH
            })

    return records, ignored


def write_warning_file(ignored, outfile="haloClassifier_input_warnings.txt"):
    """
    Write a report listing files ignored by HaloClassifier and the reason.
    """
    if not ignored:
        return

    with open(outfile, "w") as f:
        f.write("Files ignored by haloClassifier\n\n")

        for filename, reason in ignored:
            f.write(f"{filename}\n")
            f.write(f"   reason: {reason}\n\n")
