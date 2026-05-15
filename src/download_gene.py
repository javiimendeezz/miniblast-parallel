import argparse
import csv
import os
import re
from io import StringIO

from Bio import Entrez, SeqIO


def download_fasta(accession):

    handle = Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="fasta",
        retmode="text"
    )

    fasta_text = handle.read()

    handle.close()

    return fasta_text


def fasta_to_record(fasta_text):

    record = next(
        SeqIO.parse(
            StringIO(fasta_text),
            "fasta"
        )
    )

    return record


def extract_gene_name(description):

    match = re.search(
        r"\(([^()]+)\),\s*mRNA",
        description
    )

    if match:
        return match.group(1)

    return "NA"


def extract_organism(description):
    parts = description.split()

    if len(parts) >= 3:
        return f"{parts[1]} {parts[2]}"

    return "NA"


def main():

    parser = argparse.ArgumentParser(
        description="Descarga un transcrito RefSeq desde NCBI."
    )

    parser.add_argument(
        "--accession",
        required=True,
        help="Accession RefSeq (ej: NM_079109.3)"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="Email requerido por NCBI"
    )

    parser.add_argument(
        "--outdir",
        default="data"
    )

    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.tool = "MiniBlastParalelo"

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Descargando {args.accession}...")

    fasta_text = download_fasta(args.accession)

    record = fasta_to_record(fasta_text)

    sequence = str(record.seq)

    description = record.description

    organism = extract_organism(description)

    gene_name = extract_gene_name(description)

    output_path = os.path.join(
        args.outdir,
        f"{args.accession}.tsv"
    )

    with open(output_path, "w", encoding="utf-8", newline="") as out_file:

        writer = csv.writer(
            out_file,
            delimiter="\t"
        )

        writer.writerow([
            "accession",
            "organism",
            "gene_name",
            "description",
            "sequence"
        ])

        writer.writerow([
            args.accession,
            organism,
            gene_name,
            description,
            sequence
        ])

    print("Descarga completada.")
    print(f"Archivo guardado en: {output_path}")


if __name__ == "__main__":
    main()