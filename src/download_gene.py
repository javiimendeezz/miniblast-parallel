# src/download_gene_transcript.py

import argparse
import os
from Bio import Entrez


def search_refseq_transcript(gene, organism):
    query = (
        f'{gene}[Gene Name] AND "{organism}"[Organism] '
        f'AND (refseq[filter]) AND biomol_mrna[PROP]'
    )

    handle = Entrez.esearch(
        db="nuccore",
        term=query,
        retmax=5
    )
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError("No se encontraron transcritos RefSeq para ese gen.")

    return record["IdList"][0]


def download_fasta(nucleotide_id, output_path):
    handle = Entrez.efetch(
        db="nuccore",
        id=nucleotide_id,
        rettype="fasta",
        retmode="text"
    )

    fasta = handle.read()
    handle.close()

    with open(output_path, "w", encoding="utf-8") as file:
        file.write(fasta)


def main():
    parser = argparse.ArgumentParser(
        description="Descarga un transcrito RefSeq en FASTA desde NCBI."
    )

    parser.add_argument("--gene", required=True)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--email", required=True)
    parser.add_argument("--outdir", default="data")

    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.tool = "MiniBlastParalelo"

    os.makedirs(args.outdir, exist_ok=True)

    print("Buscando transcrito RefSeq...")
    nucleotide_id = search_refseq_transcript(args.gene, args.organism)

    output_path = os.path.join(
        args.outdir,
        f"{args.gene}_{args.organism}_transcript.fasta".replace(" ", "_")
    )

    print(f"ID encontrado: {nucleotide_id}")
    print("Descargando FASTA...")

    download_fasta(nucleotide_id, output_path)

    print(f"Archivo guardado en: {output_path}")


if __name__ == "__main__":
    main()