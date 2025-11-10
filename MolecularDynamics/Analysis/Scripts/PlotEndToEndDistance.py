import re
import pandas as pd
import matplotlib
import numpy as np

from pathlib import Path
from PlottingPreferences import *
from matplotlib.lines import Line2D


def load_long(csv_path):
    df = pd.read_csv(csv_path)
    val_cols = [c for c in df.columns if c.startswith("v")]

    time_w = df[df["kind"] == "time"].set_index("series")[val_cols]
    dist_w = df[df["kind"] == "distance"].set_index("series")[val_cols]

    time_l = time_w.stack().reset_index(name="time").rename(columns={"level_1": "step"})
    dist_l = dist_w.stack().reset_index(name="distance").rename(columns={"level_1": "step"})

    time_l["step"] = time_l["step"].str[1:].astype(int)
    dist_l["step"] = dist_l["step"].str[1:].astype(int)
    long = pd.merge(time_l, dist_l, on=["series", "step"], how="outer")

    long = long.dropna(subset=["time", "distance"]).reset_index(drop=True)

    m = long["series"].str.extract(r'^(?P<sample_id>HP\d+N?)_(?P<force>[\d.]+)nN_R(?P<run>\d+)$', flags=re.I)
    long["sample_id"] = m["sample_id"]
    long["force_nN"] = m["force"].astype(float)
    long["run"] = m["run"].astype(int)

    long = long[["series", "sample_id", "force_nN", "run", "time", "distance"]]
    return long


def available_series(long: pd.DataFrame) -> pd.DataFrame:
    return (
        long[["sample_id", "force_nN", "run", "series"]]
        .drop_duplicates()
        .sort_values(["sample_id", "force_nN", "run"])
        .reset_index(drop=True)
    )


def plot_subset(long: pd.DataFrame, file_name, graphics_path,
                samples=None,
                sample_regex=None,
                forces=None,
                force_min=None, force_max=None,
                runs=None,
                max_series=None):

    if len(runs) != 1:
        colors = [
            (1.0, 0.25, 0.25),
            (0.7, 0.1, 0.1),
            (0.35, 0.0, 0.0),
            (0.25, 0.9, 0.25),
            (0.1, 0.6, 0.1),
            (0.0, 0.3, 0.0)
        ]
        linestyle="-"
        linewidth=1.5
    else:
        colors = [(0, 0, 0),
             (0, 0, 0),
             (0, 0, 0),
             (0, 0, 0),
             (0, 0, 0),
             (0, 0, 0)]
        linestyle=(0, (1, 1))
        linewidth=1.5

    data = long.copy()

    if samples is not None:
        if isinstance(samples, str):
            samples = [samples]
        data = data[data["sample_id"].isin(samples)]

    if sample_regex is not None:
        data = data[data["sample_id"].str.contains(sample_regex, regex=True, na=False)]

    if forces is not None:
        if not isinstance(forces, (list, tuple, set)):
            forces = [forces]
        data = data[data["force_nN"].isin(forces)]

    if force_min is not None:
        data = data[data["force_nN"] >= float(force_min)]
    if force_max is not None:
        data = data[data["force_nN"] <= float(force_max)]

    if runs is not None:
        if not isinstance(runs, (list, tuple, set)):
            runs = [runs]
        data = data[data["run"].isin(runs)]

    apply_plot_config(wide_plot_config)
    plt.figure()
    plotted = 0
    for series_name, g in data.groupby("series", sort=False):
        g = g.sort_values("time")
        distances = np.array(g["distance"])
        distances = [x/distances[0] for x in distances]
        window = 50
        run_average = []
        for i in range(0, len(distances)-window):
            run_average.append(sum(distances[i:i+window])/window)
        if series_name.find("4N") != -1:
            color_idx = 0
        else:
            color_idx = 3

        color_idx += (int(series_name.split("R")[1])-1)
        if series_name.find("4N") != -1:
            plt.plot(g["time"][0:len(run_average)],  run_average, linestyle=linestyle, color=colors[color_idx], linewidth=linewidth)
        else:
            plt.plot(g["time"][0:len(run_average)],  run_average,color=colors[color_idx], linewidth=linewidth)

        plotted += 1
        if max_series is not None and plotted >= max_series:
            break

    plt.xticks([0, 20000, 40000, 60000, 80000, 100000], ["0", "20", "40", "60", "80", "100"] )
    plt.xlabel("Time [ns]")
    plt.ylabel("Frac. end-to-end distance")
    if len(runs) == 1:
        handles = [Line2D([], [], color='black', label='First Step'),
                   Line2D([], [], color='black', linestyle=linestyle, label='Second Step')]
    else:
        handles = [Line2D([], [], color='green', label='First Step'),
                   Line2D([], [], color='red', linestyle=linestyle, label='Second Step')]

    plt.legend(handles=handles)

    plt.tight_layout()
    plt.savefig(f"{graphics_path}/{file_name}")

matplotlib.use("Agg")


CSV_PATH = Path("../ExtractedData/ete_distances_all.csv")
GRAPHICS_PATH = Path("../Graphics")


long_data = load_long(CSV_PATH)

# Create all EndToEnd Distance Plots
for sample in ["HP004", "HP044", "HP064", "HP0104"]:
    for force in [0.1, 2.0]:
        plot_subset(long_data, samples=[sample, f"{sample}N"], forces=[force], runs=[1,2,3],
                    file_name=f"{sample}_{str(float(force))}nN_EndToEndDistances.pdf", graphics_path=GRAPHICS_PATH)

# Create one plot for the main figure
plot_subset(long_data, samples=["HP064", f"HP064N"], forces=[0.1], runs=[3],
            file_name=f"HP064_0.1nN_EndToEndDistances_Main.pdf", graphics_path=GRAPHICS_PATH)
