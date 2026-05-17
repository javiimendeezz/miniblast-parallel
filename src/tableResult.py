"""
Utilidad de benchmarking para comparar ejecuciones secuenciales y paralelas.

Define la clase BenchmarkTable, que almacena el tiempo de la version
secuencial (linea base) y de una o varias ejecuciones paralelas con
distinto numero de workers. Soporta REPETIR cada configuracion N veces
y quedarse con la MEDIANA para reducir el ruido del sistema, ademas de
calcular la desviacion estandar como medida de variabilidad.

Calcula automaticamente el speedup y la eficiencia, y permite mostrar
los resultados como tabla en consola o exportarlos a CSV.

Speedup    = T_secuencial_mediana / T_paralelo_mediana
Eficiencia = Speedup / numero_de_workers
"""

import csv
import os
import statistics
import time


class BenchmarkTable:
    """
    Tabla de benchmark que compara una ejecucion secuencial frente
    a varias ejecuciones paralelas con distinto numero de workers,
    con soporte para multiples repeticiones (mediana + desviacion).
    """

    def __init__(self, label=""):
        self.label = label
        self.sequential_time = None       # mediana de los tiempos secuenciales
        self.sequential_std = 0.0         # desviacion estandar
        self.sequential_runs = []         # lista de tiempos individuales
        self.parallel_runs = []           # lista de dicts (uno por config de workers)

    # ----------------- Registro manual -----------------

    def set_sequential_time(self, elapsed_seconds, all_times=None):
        if elapsed_seconds <= 0:
            raise ValueError("El tiempo secuencial debe ser positivo.")

        self.sequential_time = elapsed_seconds
        self.sequential_runs = all_times if all_times is not None else [elapsed_seconds]
        self.sequential_std = (
            statistics.stdev(self.sequential_runs)
            if len(self.sequential_runs) > 1 else 0.0
        )

    def add_parallel_run(self, workers, elapsed_seconds, all_times=None):
        if self.sequential_time is None:
            raise ValueError(
                "Primero hay que registrar el tiempo secuencial "
                "con set_sequential_time() o measure_sequential()."
            )
        if workers <= 0:
            raise ValueError("El numero de workers debe ser positivo.")
        if elapsed_seconds <= 0:
            raise ValueError("El tiempo paralelo debe ser positivo.")

        all_times = all_times if all_times is not None else [elapsed_seconds]
        std = statistics.stdev(all_times) if len(all_times) > 1 else 0.0

        speedup = self.sequential_time / elapsed_seconds
        efficiency = speedup / workers

        self.parallel_runs.append({
            "workers": workers,
            "time": elapsed_seconds,
            "time_std": std,
            "all_times": all_times,
            "repeats": len(all_times),
            "speedup": speedup,
            "efficiency": efficiency
        })

    # ----------------- Registro con N repeticiones automaticas -----------------

    def _run_with_repeats(self, repeats, func, args, kwargs):
        """Ejecuta func 'repeats' veces y devuelve (resultado_ultimo, lista_de_tiempos)."""
        times = []
        result = None
        for i in range(repeats):
            print(f"    Repeticion {i + 1}/{repeats}...", flush=True)
            start = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start
            times.append(elapsed)
            print(f"    -> {elapsed:.4f} s")
        return result, times

    def measure_sequential(self, func, *args, repeats=1, **kwargs):
        """
        Ejecuta func 'repeats' veces y registra la mediana como linea base.
        Devuelve (resultado_de_la_ultima_ejecucion, mediana_en_segundos).
        """
        result, times = self._run_with_repeats(repeats, func, args, kwargs)
        median_time = statistics.median(times)
        self.set_sequential_time(median_time, all_times=times)
        return result, median_time

    def measure_parallel(self, n_workers, func, *args, repeats=1, **kwargs):
        """
        Ejecuta func 'repeats' veces y registra la mediana como ejecucion paralela
        con 'n_workers' workers.
        Devuelve (resultado_de_la_ultima_ejecucion, mediana_en_segundos).
        """
        result, times = self._run_with_repeats(repeats, func, args, kwargs)
        median_time = statistics.median(times)
        self.add_parallel_run(n_workers, median_time, all_times=times)
        return result, median_time

    # ----------------- Conversion / salida -----------------

    def to_rows(self):
        rows = []

        if self.sequential_time is not None:
            rows.append({
                "mode": "Secuencial",
                "workers": 1,
                "repeats": len(self.sequential_runs),
                "time_seconds": round(self.sequential_time, 6),
                "time_std": round(self.sequential_std, 6),
                "speedup": 1.0,
                "efficiency": 1.0
            })

        for run in self.parallel_runs:
            rows.append({
                "mode": "Paralelo",
                "workers": run["workers"],
                "repeats": run["repeats"],
                "time_seconds": round(run["time"], 6),
                "time_std": round(run["time_std"], 6),
                "speedup": round(run["speedup"], 4),
                "efficiency": round(run["efficiency"], 4)
            })

        return rows

    def display(self):
        rows = self.to_rows()

        if not rows:
            print("No hay datos de benchmark que mostrar.")
            return

        if self.label:
            print(f"\n=== Benchmark: {self.label} ===")

        # Tabla final compacta: solo workers, mediana y speedup
        simple_rows = [
            {
                "workers": row["workers"],
                "time_seconds": row["time_seconds"],
                "speedup": row["speedup"]
            }
            for row in rows
        ]

        try:
            from tabulate import tabulate
            print(tabulate(
                simple_rows,
                headers={
                    "workers": "Workers",
                    "time_seconds": "Mediana (s)",
                    "speedup": "Speedup"
                },
                tablefmt="grid",
                floatfmt=".4f"
            ))
        except ImportError:
            self._display_plain(simple_rows)

    def _display_plain(self, rows):
        header = (
            f"{'Workers':>8} "
            f"{'Mediana (s)':>14} "
            f"{'Speedup':>10}"
        )
        sep = "-" * len(header)

        print(sep)
        print(header)
        print(sep)
        for row in rows:
            print(
                f"{row['workers']:>8} "
                f"{row['time_seconds']:>14.4f} "
                f"{row['speedup']:>10.4f}"
            )
        print(sep)

    def to_csv(self, output_path):
        rows = self.to_rows()

        if not rows:
            raise ValueError("No hay datos de benchmark que guardar.")

        outdir = os.path.dirname(output_path)
        if outdir:
            os.makedirs(outdir, exist_ok=True)

        with open(output_path, "w", encoding="utf-8", newline="") as file:
            fieldnames = [
                "mode", "workers", "repeats",
                "time_seconds", "time_std",
                "speedup", "efficiency"
            ]
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

# Usamos MEDIANA por: Refleja mejor cuanto tarda porque no tiene en cuenta los valores anómalos. (Ordena y coge el del medio)
# Para llamarlo: python src/run_tableResult.py --reference data/NM_079109.3.tsv --database data/sequence_database.tsv --workers-list 2,4,8,16 --repeats 5