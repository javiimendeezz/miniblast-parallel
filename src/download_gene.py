# src/download_gene.py

import argparse
import os

from Bio import Entrez


def get_nucleotide_id_from_gene(gene_id):
    """
    Obtiene un identificador de nucleótido asociado al gen.
    """

    handle = Entrez.elink(
        dbfrom="gene",
        db="nuccore",
        id=gene_id
    )

    record = Entrez.read(handle)
    handle.close()

    links = record[0]["LinkSetDb"]

    if not links:
        raise ValueError("No se encontraron secuencias asociadas al gen.")

    nucleotide_id = links[0]["Link"][0]["Id"]

    return nucleotide_id


def download_fasta(nucleotide_id, output_path):
    """
    Descarga la secuencia FASTA.
    """

    handle = Entrez.efetch(
        db="nuccore",
        id=nucleotide_id,
        rettype="fasta",
        retmode="text"
    )

    fasta_data = handle.read()

    handle.close()

    with open(output_path, "w", encoding="utf-8") as file:
        file.write(fasta_data)


def main():

    parser = argparse.ArgumentParser(
        description="Descarga secuencias FASTA desde NCBI usando Gene ID."
    )

    parser.add_argument(
        "--gene-id",
        required=True,
        help="Gene ID de NCBI"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="Email requerido por NCBI"
    )

    parser.add_argument(
        "--outdir",
        default="data",
        help="Directorio de salida"
    )

    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.tool = "MiniBlastParallel"

    os.makedirs(args.outdir, exist_ok=True)

    print("Buscando secuencia asociada al gen...")

    nucleotide_id = get_nucleotide_id_from_gene(args.gene_id)

    print(f"Nucleotide ID encontrado: {nucleotide_id}")

    output_file = os.path.join(
        args.outdir,
        f"gene_{args.gene_id}.fasta"
    )

    print("Descargando FASTA...")

    download_fasta(nucleotide_id, output_file)

    print(f"Archivo guardado en: {output_file}")


if __name__ == "__main__":
    main()