import argparse
import csv
import os
import time 
from concurrent.futures import ProcessPoolExecutor, as_completed  # Para paralelizar procesos

# Cargo el Script que generamos con la referencia
def load_reference(reference_path):
    with open(reference_path, "r", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        row = next(reader)

    return {
        "accession": row["accession"],
        "organism": row["organism"],
        "gene_name": row["gene_name"],
        "description": row["description"],
        "sequence": row["sequence"]
    }

# Carga todas los genes de Drosophila que tenemos en la base de datos
def load_database(database_path):
    sequences = []

    with open(database_path, "r", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:
            sequences.append(row)

    return sequences

# Generamos los Kmeros y le pasamos: long Kmer, cómo avanza y el numero máximo de Kmers a generar. La referencia también.
def generate_reference_kmers(reference_sequence, kmer_size, step, max_kmers):
    # Guardo en mi lista los kmers
    kmers = []

    for start in range(0, len(reference_sequence) - kmer_size + 1, step):
        kmer = reference_sequence[start:start + kmer_size]

        kmers.append({
            "kmer": kmer,
            "reference_start": start
        })

        # En cuanto supera mi secuencia de referencia, lo descarto
        if max_kmers is not None and len(kmers) >= max_kmers:
            break

    return kmers

# Para ver cuántos matches hay entre un kmer y una ventana de la secuencia objetivo
def score_window(kmer, window):
    matches = 0

    for a, b in zip(kmer, window):
        if a == b:
            matches += 1

    return matches

# Usando el score_window vemos el mejor = más matches.
def best_local_match(reference_kmers, target_sequence, window_step):
    best_matches = -1
    best_score_percent = 0.0
    best_reference_start = None
    best_target_start = None

    kmer_size = len(reference_kmers[0]["kmer"])

    if len(target_sequence) < kmer_size:
        return {
            "matches": 0,
            "score_percent": 0.0,
            "reference_start": None,
            "target_start": None
        }

    for ref_item in reference_kmers:
        kmer = ref_item["kmer"]

        for target_start in range(0, len(target_sequence) - kmer_size + 1, window_step):
            window = target_sequence[target_start:target_start + kmer_size]

            matches = score_window(kmer, window)

            if matches > best_matches:
                best_matches = matches
                best_score_percent = (matches / kmer_size) * 100
                best_reference_start = ref_item["reference_start"]
                best_target_start = target_start

    return {
        "matches": best_matches,
        "score_percent": round(best_score_percent, 4),
        "reference_start": best_reference_start,
        "target_start": best_target_start
    }

# Función que se ejecutará en paralelo para procesar cada secuencia de la base de datos
def _process_sequence(args):
    row, reference_kmers, window_step = args

    result = best_local_match(
        reference_kmers=reference_kmers,
        target_sequence=row["sequence"],
        window_step=window_step
    )

    kmer_size = len(reference_kmers[0]["kmer"]) if reference_kmers else 0

    return {
        "accession": row["accession"],
        "organism": row["organism"],
        "gene_name": row["gene_name"],
        "length": row["length"],
        "matches": result["matches"],
        "kmer_size": kmer_size,
        "score_percent": result["score_percent"],
        "reference_start": result["reference_start"],
        "target_start": result["target_start"],
        "description": row["description"]
    }

# Es la función que emplea el POOL=PARALLELO para las secuencias.
def compare_sequences_parallel(reference, database, kmer_size, reference_step, window_step, max_kmers, workers):
    reference_sequence = reference["sequence"]

    reference_kmers = generate_reference_kmers(
        reference_sequence=reference_sequence,
        kmer_size=kmer_size,
        step=reference_step,
        max_kmers=max_kmers
    )

    print(f"K-mers generados desde la referencia: {len(reference_kmers)}")

    args_list = [(row, reference_kmers, window_step) for row in database]

    results = []
    completed = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_process_sequence, args): i for i, args in enumerate(args_list)}

        for future in as_completed(futures):
            results.append(future.result())
            completed += 1

            if completed % 100 == 0:
                print(f"Procesadas {completed}/{len(database)} secuencias...", flush=True)

    return results


def save_results(results, output_path):
    sorted_results = sorted(
        results,
        key=lambda row: row["score_percent"],
        reverse=True
    )

    with open(output_path, "w", encoding="utf-8", newline="") as file:
        fieldnames = [
            "accession",
            "organism",
            "gene_name",
            "length",
            "matches",
            "kmer_size",
            "score_percent",
            "reference_start",
            "target_start",
            "description"
        ]

        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted_results)


def print_top_results(results, top_n):
    sorted_results = sorted(
        results,
        key=lambda row: row["score_percent"],
        reverse=True
    )

    print(f"\nTop {top_n} mejores coincidencias locales:\n")

    for i, row in enumerate(sorted_results[:top_n], start=1):
        print(
            f"{i}. {row['gene_name']} | "
            f"{row['organism']} | "
            f"{row['accession']} | "
            f"Score: {row['score_percent']}% | "
            f"Matches: {row['matches']}/{row['kmer_size']} | "
            f"Ref pos: {row['reference_start']} | "
            f"Target pos: {row['target_start']}"
        )

# ===================== MAIN =====================
'''
El orden que tenemos que llamar para poder comparar los datos de forma paralela sería:

1) Cargar la referencia y la base de datos.
2) Generar los kmers de la referencia. --> compare_sequences_parallel.
  2.1) Genero los kmers de la referencia.
  2.2) Para cada secuencia de la base de datos, en paralelo. LA PARALELIZACIÓN SE HACE AQUÍ.
  2.3) LLamo a process_sequence, que se encarga de llamar a best_local_match para cada secuencia de la base de datos.
3) Guardar los resultados en un archivo de salida.
4) Imprimir los mejores resultados.
'''
def main():
    parser = argparse.ArgumentParser(
        description="Mini-BLAST paralelo con búsqueda local por ventanas móviles."
    )

    parser.add_argument("--reference", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--output", default="miniblast_window_parallel_results.tsv")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--workers", type=int, default=os.cpu_count())

    parser.add_argument("--kmer-size", type=int, default=50)
    parser.add_argument("--reference-step", type=int, default=25)
    parser.add_argument("--window-step", type=int, default=5)
    parser.add_argument("--max-kmers", type=int, default=50)

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    output_path = os.path.join(args.outdir, args.output)

    reference = load_reference(args.reference)
    database = load_database(args.database)

    print("Iniciando MiniBLAST paralelo por ventanas...")
    print(f"Referencia: {reference['gene_name']} | {reference['organism']} | {reference['accession']}")
    print(f"Longitud referencia: {len(reference['sequence'])} bp")
    print(f"Secuencias en base de datos: {len(database)}")
    print(f"k-mer size: {args.kmer_size}")
    print(f"reference step: {args.reference_step}")
    print(f"window step: {args.window_step}")
    print(f"max k-mers: {args.max_kmers}")
    print(f"workers: {args.workers}")

    start_time = time.perf_counter()

    results = compare_sequences_parallel(
        reference=reference,
        database=database,
        kmer_size=args.kmer_size,
        reference_step=args.reference_step,
        window_step=args.window_step,
        max_kmers=args.max_kmers,
        workers=args.workers
    )

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time

    save_results(results, output_path)
    print_top_results(results, args.top)

    print(f"\nTiempo de ejecución paralela: {elapsed_time:.6f} segundos")
    print(f"Workers utilizados: {args.workers}")
    print(f"Resultados guardados en: {output_path}")


if __name__ == "__main__":
    main()
