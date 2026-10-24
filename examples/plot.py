"""Plot the quickstart outputs; published phase diagrams require ensembles."""
import argparse
import csv
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    with (args.input / "metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    length = metadata["L"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    for ax, model, color in zip(axes, ("cgpm", "spmm"), ("#008875", "#702e96")):
        params = metadata["models"][model]
        label = f"{model.upper()}: J = {params['strength']:g}, p = {params['p']:g}"
        with (args.input / f"{model}.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        t = np.array([float(r["layer"]) for r in rows])
        s = np.array([float(r["halfcut_entropy_bits"]) for r in rows])
        if (len(t) != metadata["layers"] or not np.isfinite(t).all()
                or not np.isfinite(s).all() or not np.all(np.diff(t) > 0)):
            raise ValueError("invalid trace")
        ax.plot(t, s, color=color, linewidth=1.8)
        ax.set(xlabel="Gate + measurement layer", ylabel="Half-cut entropy (bits)",
               title=label, ylim=(0, length / 2))
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Small-system examples: L = {length}, one trajectory per model", fontsize=12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("xb") as stream:
        fig.savefig(stream, format="png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
