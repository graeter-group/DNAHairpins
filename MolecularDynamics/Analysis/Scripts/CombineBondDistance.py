import re
import csv
from pathlib import Path

input_dir = Path("../ExtractedData")
output_file = "../ExtractedData/backbone_bonds_all.csv"

rows = []
HDR_RE = re.compile(r"^(?P<bond>.+)\.Res(?P<res>[0-9.]+)\.(?P<strand>[A-Za-z]+):\s*$")
AVG_RE = re.compile(r"Average distance:\s+([0-9.eE+-]+)")

for txtfile in input_dir.glob("*.txt"):
    sample_name = txtfile.stem.strip("_BackboneBonds")
    bond_type = residue_index = strand = None

    with txtfile.open("r", encoding="utf-8") as f:
        for line in f:
            m = HDR_RE.match(line.strip())
            if m:
                bond_type = m.group("bond")
                residue_index = m.group("res")
                strand = m.group("strand")
                continue

            m = AVG_RE.search(line)
            if m and bond_type is not None:
                avg = m.group(1)
                rows.append([sample_name, bond_type, residue_index, strand, avg])

with open(output_file, "w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out)
    writer.writerow(["SampleName", "BondType", "ResidueIndex", "Strand", "AverageDistance"])
    writer.writerows(rows)

print(f"Wrote {len(rows)} rows to {output_file}")