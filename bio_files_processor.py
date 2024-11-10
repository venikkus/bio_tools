def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    Converts a multi-line FASTA file to a single-line FASTA format.

    This function reads a FASTA file where sequences might be spread across
    multiple lines and outputs a new FASTA file or prints the sequences, where
    each sequence is written on a single line.

    Parameters
    ----------
    input_fasta : str
        Path to the input multi-line FASTA file.
    output_fasta : str, optional
        Path to the output FASTA file. If not provided, the results will be
        printed to the console.

    Returns
    -------
    None
    """
    with open(input_fasta, "r") as infile:
        sequences = {}
        header = None
        sequence = []

        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)

        sequences[header] = "".join(sequence) if header else None

    if output_fasta:
        with open(output_fasta, "w") as outfile:
            for header, seq in sequences.items():
                outfile.write(f"{header}\n{seq}\n")
    else:
        for header, seq in sequences.items():
            print(f"{header}\n{seq}")


def parse_blast_output(input_file, output_file=None):
    """
    Parses the output of a BLAST search and extracts the first description line
    for each query.

    This function reads a BLAST results file, extracts descriptions
    for each query, and outputs the sorted results to a file or prints
    them to the console.

    Parameters
    ----------
    input_file : str
        Path to the input BLAST results file.
    output_file : str, optional
        Path to the output file. If not provided, the results will be printed
        to the console.

    Returns
    -------
    None
    """
    descriptions = []
    query_found = False

    with open(input_file, "r") as infile:
        for line in infile:
            line = line.strip()

            if line.startswith("Query #"):
                query_found = True
                continue

            if query_found:
                if ("]" in line) and not line.startswith("Description"):
                    description = line.split("]")[0].strip() + "]"
                    descriptions.append(description)
                    query_found = False

                elif ("..." in line) and not line.startswith("Description"):
                    if "...." in line:
                        description = line.split("....")[0].strip() + "...."
                    else:
                        description = line.split("...")[0].strip() + "..."
                    descriptions.append(description)
                    query_found = False
                continue

    descriptions = sorted(descriptions)

    if output_file:
        with open(output_file, "w") as outfile:
            for desc in descriptions:
                outfile.write(f"{desc}\n")
    else:
        for desc in descriptions:
            print(desc)


def select_genes_from_gbk_to_fasta(
    input_gbk,
    genes,
    n_before=1,
    n_after=1,
    output_fasta=None
):
    """
    Extracts specified gene sequences along with n_before upstream
    and n_after downstream genes from a GenBank file.

    Parameters
    ----------
    input_gbk : str
        Path to the input GenBank (.gbk) file containing genomic data.
    genes : list of str
        List of gene names to extract from the GenBank file.
    n_before : int, optional, default: 1
        Number of upstream genes (before the selected gene) to include in
        the output.
    n_after : int, optional, default: 1
        Number of downstream genes (after the selected gene) to include in
        the output.
    output_fasta : str, optional
        Path to the output FASTA file. If not provided, the sequences will
        be printed.
    """
    gene_records = []
    current_record = []
    translation_lines = []
    query_found = False

    with open(input_gbk, "r") as infile:
        for line in infile:
            line = line.strip()

            if "/gene=" in line:
                start, end = line.find('"'), line.rfind('"')
                gene_name = line[start + 1:end]
                current_record = [gene_name]
                translation_lines = []
                query_found = True
                continue

            if query_found and "/translation=" in line:
                start = line.find('"') + 1
                translation_lines.append(line[start:].replace('"', ""))
                continue

            if query_found and translation_lines and not line.endswith('"'):
                translation_lines.append(line.strip())
                continue

            if query_found and translation_lines and line.endswith('"'):
                full_translation = "".join(translation_lines)
                current_record.append(full_translation)
                gene_records.append(tuple(current_record))
                query_found = False
                continue

    matching_genes = []
    for idx, (gene_name, translation) in enumerate(gene_records):
        if any(gene in gene_name for gene in genes):
            start_idx = max(0, idx - n_before)
            end_idx = min(len(gene_records), idx + n_after + 1)
            matching_genes.extend(gene_records[start_idx:end_idx])

    if output_fasta:
        with open(output_fasta, "w") as outfile:
            for gene_name, translation in matching_genes:
                outfile.write(f">{gene_name}\n{translation}\n")
    else:
        for gene_name, translation in matching_genes:
            print(f">{gene_name}\n{translation}\n")
