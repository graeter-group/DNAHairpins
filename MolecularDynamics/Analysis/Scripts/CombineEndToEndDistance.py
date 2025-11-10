from pathlib import Path
import re
import pandas as pd

INPUT_DIR = Path("../ExtractedData/")
OUTPUT_CSV = INPUT_DIR / "ete_distances_all.csv"

FILENAME_RE = re.compile(r'^(?P<sample_id>HP\d+N?)_(?P<force>[\d.]+)nN_(?P<run>R\d+)', re.IGNORECASE)

def parse_xvg(xvg_path: Path):
    times, vals = [], []
    with xvg_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith(("@", "#")):
                continue
            parts = line.split()
            try:
                nums = [float(tok) for tok in parts]
            except ValueError:
                continue
            if len(nums) >= 2:
                times.append(nums[0])
                vals.append(nums[-1])
    return times, vals

def extract_meta(xvg_path: Path):
    m = FILENAME_RE.match(xvg_path.stem)
    if not m:
        raise ValueError(f"Filename does not match pattern: {xvg_path.name}")
    sample_id = m.group("sample_id")
    force = float(m.group("force"))
    run = int(m.group("run")[1:])
    label = f"{sample_id}_{m.group('force')}nN_R{run}"
    return sample_id, force, run, label

def pad_row(prefix_cells, values, target_len):
    row = list(prefix_cells) + list(values)
    row.extend([''] * (target_len - len(values)))
    return row

def main():
    xvg_files = sorted(INPUT_DIR.glob("*.xvg"))
    if not xvg_files:
        raise SystemExit(f"No .xvg files found in {INPUT_DIR}")

    series = []
    max_len = 0

    for xvg in xvg_files:
        sample_id, force, run, label = extract_meta(xvg)
        times, dists = parse_xvg(xvg)
        max_len = max(max_len, len(times), len(dists))
        series.append(((sample_id, force, run), label, times, dists))

    series.sort(key=lambda t: (t[0][0], t[0][1], t[0][2]))

    rows = []
    header = ["series", "kind"] + [f"v{i}" for i in range(max_len)]

    for (_keys, label, times, dists) in series:
        rows.append(pad_row([label, "time"], times, max_len))
        rows.append(pad_row([label, "distance"], dists, max_len))

    df = pd.DataFrame(rows, columns=header)
    df.to_csv(OUTPUT_CSV, index=False)

    print(f"Wrote {OUTPUT_CSV}")
    print(f"Rows: {len(rows)} (2 per file); Columns: {len(header)}")

if __name__ == "__main__":
    main()
