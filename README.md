# MiniBLAST Parallel

Implementación secuencial y paralela (PySpark) de un algoritmo de alineamiento local por ventanas deslizantes, similar a BLAST.

## Requisitos previos

- **Python 3.9+**
- **Java JDK 11 o superior** — necesario para PySpark (pip no lo instala)
  - Descarga: https://adoptium.net/
  - Verificar instalación: `java -version`
  - En Windows, asegúrate de que `JAVA_HOME` apunta al JDK instalado

## Instalación

```bash
# Crear y activar entorno virtual
python -m venv .venv

# Windows
.venv\Scripts\activate

# Linux / macOS
source .venv/bin/activate

# Instalar dependencias
pip install -r requirements.txt
```

## Uso

### Versión secuencial

```bash
python src/miniblast_sequential.py \
    --reference data/NM_079109.3.tsv \
    --database  data/sequence_database.tsv
```

### Versión paralela (PySpark)

```bash
python src/miniblast_parallel.py \
    --reference data/NM_079109.3.tsv \
    --database  data/sequence_database.tsv \
    --workers   4
```

### Benchmark secuencial vs paralelo

```bash
python src/run_tableResult.py \
    --reference   data/NM_079109.3.tsv \
    --database    data/sequence_database.tsv \
    --workers-list 2,4,8,16 \
    --repeats     5
```

Genera `results/benchmark_results.csv` con tiempos, speedup y eficiencia.

### Gráficas

```bash
python src/plot_benchmark.py
```

Genera en `results/`:
- `benchmark_speedup.png` — curva de speedup vs número de workers
- `benchmark_times.png` — tiempos de ejecución por configuración

## Estructura

```
data/          Referencia y base de datos en formato TSV
src/
  miniblast_sequential.py   Algoritmo secuencial
  miniblast_parallel.py     Algoritmo paralelo con PySpark
  run_tableResult.py        Benchmark comparativo
  plot_benchmark.py         Generación de gráficas
  download_gene.py          Descarga de secuencias desde NCBI
results/       Salida de resultados y gráficas
```
