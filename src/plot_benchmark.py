"""
Genera dos graficas PNG con seaborn a partir del CSV de benchmark:

  1. Grafica LINEAL: Speedup vs numero de workers (con linea ideal de
     referencia y_= x). Muestra como escala el paralelismo.

  2. Grafica de BARRAS: Tiempo de ejecucion de cada configuracion
     (secuencial en color distinto al paralelo). Visualiza el ahorro
     real de tiempo.

Uso de ejemplo:
    python src/plot_benchmark.py
    python src/plot_benchmark.py --csv results/benchmark_results.csv \
                                 --line-output results/speedup.png \
                                 --bar-output  results/times.png
"""

import argparse
import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# --------------------------------------------------------------------
# 1) GRAFICA LINEAL: speedup vs workers
# --------------------------------------------------------------------
def plot_speedup_line(df, output_path, title):
    sns.set_theme(style="whitegrid", palette="deep", font_scale=1.1)

    fig, ax = plt.subplots(figsize=(9, 6))

    # Curva real
    sns.lineplot(
        data=df,
        x="workers",
        y="speedup",
        marker="o",
        markersize=10,
        linewidth=2.5,
        label="Speedup real",
        ax=ax
    )

    # Curva ideal y = x
    max_workers = int(df["workers"].max())
    ideal_x = list(range(1, max_workers + 1))
    ax.plot(
        ideal_x,
        ideal_x,
        linestyle="--",
        color="gray",
        linewidth=1.5,
        label="Speedup ideal (lineal)"
    )

    # Anotaciones sobre cada punto: speedup y tiempo
    for _, row in df.iterrows():
        label = f"{row['speedup']:.2f}x\n({row['time_seconds']:.2f}s)"
        ax.annotate(
            label,
            xy=(row["workers"], row["speedup"]),
            xytext=(8, 8),
            textcoords="offset points",
            fontsize=9,
            ha="left"
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Numero de workers (procesos)", fontsize=12)
    ax.set_ylabel("Speedup (T_secuencial / T_paralelo)", fontsize=12)
    ax.set_xticks(df["workers"].tolist())
    ax.legend(loc="upper left", frameon=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Grafica lineal guardada en: {output_path}")


# --------------------------------------------------------------------
# 2) GRAFICA DE BARRAS: tiempo de ejecucion por configuracion
# --------------------------------------------------------------------
def plot_times_bar(df, output_path, title):
    sns.set_theme(style="whitegrid", font_scale=1.1)

    fig, ax = plt.subplots(figsize=(9, 6))

    # Color distinto para secuencial vs paralelo
    colors = []
    labels = []
    for _, row in df.iterrows():
        is_sequential = str(row["mode"]).lower().startswith("sec")
        if is_sequential:
            colors.append("#c95136")          # naranja para secuencial
            labels.append("Secuencial\n(1 worker)")
        else:
            colors.append("#5dc6db")          # azul para paralelo
            labels.append(f"Paralelo\n({int(row['workers'])} workers)")

    bars = ax.bar(
        labels,
        df["time_seconds"],
        color=colors,
        edgecolor="black",
        linewidth=0.6
    )

    # Etiqueta encima de cada barra con el tiempo en segundos
    for bar, secs in zip(bars, df["time_seconds"]):
        ax.annotate(
            f"{secs:.2f}s",
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            fontsize=10,
            fontweight="bold"
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Configuracion", fontsize=12)
    ax.set_ylabel("Tiempo de ejecucion (s)", fontsize=12)

    # Leyenda manual con los dos colores
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#c95136", edgecolor="black", label="Secuencial"),
        Patch(facecolor="#5dc6db", edgecolor="black", label="Paralelo"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Grafica de barras guardada en: {output_path}")


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Genera dos graficas PNG (lineal y de barras) a partir del CSV de benchmark."
    )
    parser.add_argument(
        "--csv",
        default="results/benchmark_results.csv",
        help="Ruta al CSV con los resultados del benchmark."
    )
    parser.add_argument(
        "--line-output",
        default="results/benchmark_speedup.png",
        help="Ruta del PNG con la grafica lineal de speedup."
    )
    parser.add_argument(
        "--bar-output",
        default="results/benchmark_times.png",
        help="Ruta del PNG con la grafica de barras de tiempos."
    )
    parser.add_argument(
        "--line-title",
        default= "Speedup vs Numero de workers",
        help="Titulo de la grafica lineal."
    )
    parser.add_argument(
        "--bar-title",
        default="Tiempos de ejecucion por configuracion de hebras",
        help="Titulo de la grafica de barras."
    )
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        raise FileNotFoundError(
            f"No se encontro el CSV: {args.csv}. "
            "Ejecuta primero run_tableResult.py para generarlo."
        )

    # Cargar y ordenar por workers (asi la curva queda bien)
    df = pd.read_csv(args.csv)
    df = df.sort_values("workers").reset_index(drop=True)

    # Asegurar carpetas de salida
    for path in (args.line_output, args.bar_output):
        outdir = os.path.dirname(path)
        if outdir:
            os.makedirs(outdir, exist_ok=True)

    # Generar las dos graficas
    plot_speedup_line(df, args.line_output, args.line_title)
    plot_times_bar(df, args.bar_output, args.bar_title)


if __name__ == "__main__":
    main()
    
# Para llamarlo: python src/plot_benchmark.py