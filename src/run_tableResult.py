"""
Ejecuta MiniBLAST secuencial y paralelo (con varios numeros de workers),
REPITIENDO cada configuracion N veces y quedandose con la MEDIANA para
calcular el speedup. Muestra una tabla con tiempos, speedup y eficiencia
y la guarda en CSV.

Uso de ejemplo:
    python src/run_tableResult.py \
        --reference data/NM_079109.3.tsv \
        --database data/sequence_database.tsv \
        --workers-list 2,4,8,16 \
        --repeats 5
"""

import argparse
import os

from tableResult import BenchmarkTable
from miniblast_sequential import (
    load_reference,
    load_database,
    compare_sequences
)
from miniblast_parallel import compare_sequences_parallel


def parse_workers(workers_str):
    """Convierte '1,2,4,8' en [1, 2, 4, 8]."""
    return [int(w.strip()) for w in workers_str.split(",") if w.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark MiniBLAST secuencial vs paralelo (con mediana de N repeticiones)."
    )

    parser.add_argument("--reference", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--outdir", default="results")
    parser.add_argument(
        "--workers-list",
        default="1,2,4,8",
        help="Lista de workers a probar separados por coma (ej: '1,2,4,8')."
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=5,
        help="Numero de repeticiones por configuracion (se toma la mediana). Por defecto 5."
    )
    parser.add_argument(
        "--csv",
        default="benchmark_results.csv",
        help="Nombre del archivo CSV donde guardar la tabla."
    )

    parser.add_argument("--kmer-size", type=int, default=50)
    parser.add_argument("--reference-step", type=int, default=25)
    parser.add_argument("--window-step", type=int, default=5)
    parser.add_argument("--max-kmers", type=int, default=50)

    args = parser.parse_args()

    workers_list = parse_workers(args.workers_list)

    os.makedirs(args.outdir, exist_ok=True)
    csv_path = os.path.join(args.outdir, args.csv)

    print("Cargando referencia y base de datos...")
    reference = load_reference(args.reference)
    database = load_database(args.database)

    print(f"Referencia: {reference['gene_name']} | {reference['accession']}")
    print(f"Base de datos: {len(database)} secuencias")
    print(f"Configuraciones a probar: secuencial + paralelo con workers={workers_list}")
    print(f"Repeticiones por configuracion: {args.repeats} (se toma la mediana)\n")

    bench = BenchmarkTable(label="MiniBLAST - Secuencial vs Paralelo")

    total = 1 + len(workers_list)

    # 1) Secuencial (linea base)
    print(f"[1/{total}] Ejecutando version SECUENCIAL ({args.repeats} reps)...")
    _, t_seq = bench.measure_sequential(
        compare_sequences,
        reference=reference,
        database=database,
        kmer_size=args.kmer_size,
        reference_step=args.reference_step,
        window_step=args.window_step,
        max_kmers=args.max_kmers,
        repeats=args.repeats
    )
    print(f"  Mediana secuencial: {t_seq:.4f} s\n")

    # 2) Paralelo con distintos workers
    for i, workers in enumerate(workers_list, start=2):
        print(f"[{i}/{total}] Ejecutando version PARALELA con {workers} workers ({args.repeats} reps)...")
        _, t_par = bench.measure_parallel(
            workers,
            compare_sequences_parallel,
            reference=reference,
            database=database,
            kmer_size=args.kmer_size,
            reference_step=args.reference_step,
            window_step=args.window_step,
            max_kmers=args.max_kmers,
            workers=workers,
            repeats=args.repeats
        )
        print(f"  Mediana paralelo ({workers} workers): {t_par:.4f} s")
        print(f"  Speedup vs secuencial: {t_seq / t_par:.2f}x\n")

    # 3) Mostrar y guardar
    bench.display()
    bench.to_csv(csv_path)
    print(f"\nTabla guardada en: {csv_path}")


if __name__ == "__main__":
    main()