import argparse
import csv
import os
import sys
import time

# Los workers de PySpark se lanzan como subprocesos Python independientes.
# Sin esto, en Windows buscan "python" en PATH y no encuentran el venv,
# fallando con "Python worker exited unexpectedly (crashed)".
os.environ["PYSPARK_PYTHON"] = sys.executable
os.environ["PYSPARK_DRIVER_PYTHON"] = sys.executable

from pyspark.sql import SparkSession


# ===================== CARGA DE DATOS =====================

def load_reference(reference_path):
    with open(reference_path, "r", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        row = next(reader)

    return {
        "accession": row["accession"],
        "organism":  row["organism"],
        "gene_name": row["gene_name"],
        "description": row["description"],
        "sequence":  row["sequence"]
    }


# ===================== ALGORITMO =====================

def generate_reference_kmers(reference_sequence, kmer_size, step, max_kmers):
    """Fragmenta la referencia en k-mers con desplazamiento fijo."""
    kmers = []

    for start in range(0, len(reference_sequence) - kmer_size + 1, step):
        kmer = reference_sequence[start:start + kmer_size]
        kmers.append({"kmer": kmer, "reference_start": start})

        if max_kmers is not None and len(kmers) >= max_kmers:
            break

    return kmers


def score_window(kmer, window):
    """Cuenta coincidencias posición a posición entre un k-mer y una ventana."""
    matches = 0
    for a, b in zip(kmer, window):
        if a == b:
            matches += 1
    return matches


def best_local_match(reference_kmers, target_sequence, window_step):
    """
    Desliza cada k-mer sobre target_sequence y devuelve el mejor
    alineamiento local encontrado (mayor número de matches).
    """
    best_matches      = -1
    best_score_pct    = 0.0
    best_ref_start    = None
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
        for target_start in range(
            0, len(target_sequence) - kmer_size + 1, window_step
        ):
            window  = target_sequence[target_start:target_start + kmer_size]
            matches = score_window(kmer, window)

            if matches > best_matches:
                best_matches      = matches
                best_score_pct    = (matches / kmer_size) * 100
                best_ref_start    = ref_item["reference_start"]
                best_target_start = target_start

    return {
        "matches":         best_matches,
        "score_percent":   round(best_score_pct, 4),
        "reference_start": best_ref_start,
        "target_start":    best_target_start
    }


# ===================== MOTOR PYSPARK =====================

_JVM_OPTIONS = (
    "--add-opens=java.base/sun.nio.ch=ALL-UNNAMED "
    "--add-opens=java.base/java.lang=ALL-UNNAMED "
    "--add-opens=java.base/java.lang.invoke=ALL-UNNAMED "
    "--add-opens=java.base/java.io=ALL-UNNAMED "
    "--add-opens=java.base/java.net=ALL-UNNAMED "
    "--add-opens=java.base/java.nio=ALL-UNNAMED "
    "--add-opens=java.base/java.util=ALL-UNNAMED "
    "--add-opens=java.base/java.util.concurrent=ALL-UNNAMED "
    "--add-opens=java.base/java.util.concurrent.atomic=ALL-UNNAMED "
    "-Djava.security.manager=allow"
)


def _get_or_create_spark(workers):
    spark = (
        SparkSession.builder
        .appName("MiniBLAST-Paralelo")
        .master(f"local[{workers}]")
        .config("spark.sql.shuffle.partitions", str(workers))
        .config("spark.ui.showConsoleProgress", "false")
        .config("spark.driver.extraJavaOptions", _JVM_OPTIONS)
        .config("spark.executor.extraJavaOptions", _JVM_OPTIONS)
        .getOrCreate()
    )
    spark.sparkContext.setLogLevel("ERROR")
    return spark
'''
Flujo de comparación paralela con PySpark:

1) Inicialización: SparkSession en modo local (usa todos los cores disponibles).
2) Broadcast: los k-mers de referencia se envían una sola vez a todos los
   executors, evitando tráfico de red redundante.
3) UDF: una función de usuario (UDF) calcula el mejor alineamiento local
   para cada secuencia del DataFrame.
4) Transformación: el DataFrame de la BD se enriquece con los resultados
   de la UDF de forma distribuida.
5) Acción (collect): los resultados se ordenan por score y se recogen
   en el driver para guardarlos en TSV.
'''


def compare_sequences_spark(
    spark, reference, database_path,
    kmer_size, reference_step, window_step, max_kmers, workers
):
    sc = spark.sparkContext

    # Funciones locales: cloudpickle las serializa por valor (bytecode),
    # no por referencia al módulo, evitando 'No module named miniblast_parallel'
    # en los workers.
    def _score_window(kmer, window):
        matches = 0
        for a, b in zip(kmer, window):
            if a == b:
                matches += 1
        return matches

    def _best_local_match(reference_kmers, target_sequence, window_step):
        best_matches = -1
        best_score_pct = 0.0
        best_ref_start = None
        best_target_start = None
        kmer_size = len(reference_kmers[0]["kmer"])
        if len(target_sequence) < kmer_size:
            return {"matches": 0, "score_percent": 0.0,
                    "reference_start": None, "target_start": None}
        for ref_item in reference_kmers:
            kmer = ref_item["kmer"]
            for target_start in range(0, len(target_sequence) - kmer_size + 1, window_step):
                window = target_sequence[target_start:target_start + kmer_size]
                matches = _score_window(kmer, window)
                if matches > best_matches:
                    best_matches = matches
                    best_score_pct = (matches / kmer_size) * 100
                    best_ref_start = ref_item["reference_start"]
                    best_target_start = target_start
        return {
            "matches": best_matches,
            "score_percent": round(best_score_pct, 4),
            "reference_start": best_ref_start,
            "target_start": best_target_start,
        }

    # --- Generar k-mers ---
    reference_kmers = generate_reference_kmers(
        reference_sequence=reference["sequence"],
        kmer_size=kmer_size,
        step=reference_step,
        max_kmers=max_kmers
    )
    print(f"K-mers generados desde la referencia: {len(reference_kmers)}")

    kmer_size_val = len(reference_kmers[0]["kmer"]) if reference_kmers else 0

    # --- Broadcast: los k-mers viajan una sola vez a cada executor ---
    broadcast_kmers = sc.broadcast(reference_kmers)
    broadcast_window_step = sc.broadcast(window_step)
    broadcast_kmer_size_val = sc.broadcast(kmer_size_val)

    # --- Cargar la base de datos en el driver y paralelizar como RDD ---
    # Se evita el framework de UDFs de PySpark (inestable en Windows) usando
    # RDD.map directamente, que serializa funciones Python puras sin overhead
    # de la capa SQL/Arrow.
    with open(database_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        db_rows = list(reader)

    # --- Transformación distribuida con RDD.map ---
    def process_row(row):
        kmers = broadcast_kmers.value
        wstep = broadcast_window_step.value
        ksz   = broadcast_kmer_size_val.value
        sequence = row.get("sequence") or ""

        if sequence:
            result = _best_local_match(kmers, sequence, wstep)
        else:
            result = {
                "matches": 0,
                "score_percent": 0.0,
                "reference_start": None,
                "target_start": None,
            }

        return {
            "accession":       row.get("accession"),
            "organism":        row.get("organism"),
            "gene_name":       row.get("gene_name"),
            "length":          row.get("length"),
            "description":     row.get("description"),
            "matches":         result["matches"],
            "kmer_size":       ksz,
            "score_percent":   result["score_percent"],
            "reference_start": result["reference_start"],
            "target_start":    result["target_start"],
        }

    rdd = sc.parallelize(db_rows, workers)

    # --- Acción: recoger en el driver y ordenar localmente ---
    results = sorted(
        rdd.map(process_row).collect(),
        key=lambda r: r["score_percent"] if r["score_percent"] is not None else 0.0,
        reverse=True
    )

    return results


def compare_sequences_parallel(reference, database_path, kmer_size, reference_step, window_step, max_kmers, workers):
    """Wrapper con la misma firma que compare_sequences (secuencial) para el benchmark."""
    spark = _get_or_create_spark(workers)
    return compare_sequences_spark(
        spark=spark,
        reference=reference,
        database_path=os.path.abspath(database_path),
        kmer_size=kmer_size,
        reference_step=reference_step,
        window_step=window_step,
        max_kmers=max_kmers,
        workers=workers,
    )


# ===================== SALIDA =====================

def save_results(results, output_path):
    sorted_results = sorted(
        results, key=lambda r: r["score_percent"], reverse=True
    )

    with open(output_path, "w", encoding="utf-8", newline="") as file:
        fieldnames = ["accession", "organism", "gene_name", "length", "matches", "kmer_size", "score_percent", "reference_start", "target_start", "description"]
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted_results)


def print_top_results(results, top_n):
    print(f"\nTop {top_n} mejores coincidencias locales:\n")
    for i, row in enumerate(results[:top_n], start=1):
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

def main():
    parser = argparse.ArgumentParser(
        description="Mini-BLAST paralelo con PySpark y búsqueda local por ventanas móviles."
    )
    parser.add_argument("--reference", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--output", default="miniblast_window_parallel_results.tsv")
    parser.add_argument("--top",  type=int, default=10)
    parser.add_argument("--workers", type=int, default=os.cpu_count())

    parser.add_argument("--kmer-size", type=int, default=50)
    parser.add_argument("--reference-step", type=int, default=25)
    parser.add_argument("--window-step", type=int, default=5)
    parser.add_argument("--max-kmers", type=int, default=50)

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    output_path = os.path.join(args.outdir, args.output)

    spark = _get_or_create_spark(args.workers)

    reference = load_reference(args.reference)

    print("Iniciando MiniBLAST paralelo con PySpark...")
    print(f"Referencia: {reference['gene_name']} | {reference['organism']} | {reference['accession']}")
    print(f"Longitud referencia: {len(reference['sequence'])} bp")
    print(f"k-mer size: {args.kmer_size}  |  reference step: {args.reference_step}")
    print(f"window step: {args.window_step}  |  max k-mers: {args.max_kmers}")
    print(f"workers: {args.workers}")

    start_time = time.perf_counter()

    results = compare_sequences_spark(
        spark=spark,
        reference=reference,
        database_path=os.path.abspath(args.database),
        kmer_size=args.kmer_size,
        reference_step=args.reference_step,
        window_step=args.window_step,
        max_kmers=args.max_kmers,
        workers=args.workers,
    )

    elapsed = time.perf_counter() - start_time

    save_results(results, output_path)
    print_top_results(results, args.top)

    print(f"\nTiempo de ejecución paralela (PySpark): {elapsed:.6f} segundos")
    print(f"Workers utilizados: {args.workers}")
    print(f"Resultados guardados en: {output_path}")

    spark.stop()


if __name__ == "__main__":
    main()
