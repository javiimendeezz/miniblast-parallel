import argparse
import csv
import os
import re
import time
from io import StringIO

from Bio import Entrez, SeqIO


def search_transcripts(organism, min_length, max_length, max_results):
    query = (
        f'"{organism}"[Organism] '
        f'AND refseq[filter] '
        f'AND biomol_mrna[PROP] '
        f'AND {min_length}:{max_length}[SLEN]'
    )

    handle = Entrez.esearch(
        db="nuccore",
        term=query,
        retmax=max_results
    )
    record = Entrez.read(handle)
    handle.close()

    return record["IdList"]


def fetch_summaries(ids):
    if not ids:
        return []

    handle = Entrez.esummary(
        db="nuccore",
        id=",".join(ids)
    )
    summaries = Entrez.read(handle)
    handle.close()

    return summaries


def download_fasta_by_id(ncbi_id):
    handle = Entrez.efetch(
        db="nuccore",
        id=ncbi_id,
        rettype="fasta",
        retmode="text"
    )
    fasta_text = handle.read()
    handle.close()

    return fasta_text


def fasta_to_sequence(fasta_text):
    record = next(SeqIO.parse(StringIO(fasta_text), "fasta"))
    return str(record.seq)


def extract_gene_name(title):
    """
    Extrae el símbolo del gen desde títulos tipo:
    Drosophila melanogaster doublesex (dsx), transcript variant G, mRNA
    Drosophila melanogaster oxoglutarate dehydrogenase, transcript variant K (Ogdh), mRNA
    """

    matches = re.findall(r"\(([^()]+)\)", title)

    if matches:
        return matches[-1]

    return title


def load_organisms(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        organisms = [
            line.strip()
            for line in file
            if line.strip()
        ]

    return organisms


def main():
    parser = argparse.ArgumentParser(
        description="Busca transcritos RefSeq mRNA y genera una base de datos TSV."
    )

    parser.add_argument(
        "--organisms-file",
        required=True,
        help="Fichero TXT con un organismo por línea"
    )

    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--max-length", type=int, required=True)
    parser.add_argument("--max-results", type=int, default=20)
    parser.add_argument("--email", required=True)
    parser.add_argument("--outdir", default="data")
    parser.add_argument("--output", default="sequence_database.tsv")

    args = parser.parse_args()

    organisms = load_organisms(args.organisms_file)

    Entrez.email = args.email
    Entrez.tool = "MiniBlastParalelo"

    os.makedirs(args.outdir, exist_ok=True)

    output_path = os.path.join(args.outdir, args.output)

    with open(output_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")

        writer.writerow([
            "accession",
            "organism",
            "gene_name",
            "length",
            "description",
            "sequence"
        ])

        for organism in organisms:
            print(f"\nBuscando transcritos en: {organism}")

            ids = search_transcripts(
                organism=organism,
                min_length=args.min_length,
                max_length=args.max_length,
                max_results=args.max_results
            )

            summaries = fetch_summaries(ids)

            print(f"Transcritos encontrados: {len(summaries)}")

            for item in summaries:
                ncbi_id = item["Id"]
                accession = item.get("AccessionVersion", "NA")
                title = item.get("Title", "NA")
                length = int(item.get("Length", 0))
                gene_name = extract_gene_name(title)

                print(
                    f"Descargando {accession} | "
                    f"{organism} | "
                    f"{gene_name} | "
                    f"{length} bp"
                )

                fasta_text = download_fasta_by_id(ncbi_id)
                sequence = fasta_to_sequence(fasta_text)

                writer.writerow([
                    accession,
                    organism,
                    gene_name,
                    length,
                    title,
                    sequence
                ])

                time.sleep(0.34)

    print(f"\nBase de datos generada correctamente: {output_path}")


if __name__ == "__main__":
    main()