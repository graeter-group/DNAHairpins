import re
import pandas as pd
import matplotlib
from PlottingPreferences import *

matplotlib.use("Agg")

def convert_to_forces(distances_in, bond_type_in):
    bond_names = ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"]
    equ_dists = [0.161, 0.141, 0.1526, 0.1562, 0.141, 0.161]
    force_constants = [192464, 267776, 259408, 259408, 267776, 192464]

    force_conversion = 10 ** (-2) / 6.022140
    force_constants = [a * force_conversion for a in force_constants]

    force_constants_dict = dict(zip(bond_names, force_constants))
    equ_dists_dict = dict(zip(bond_names, equ_dists))

    force_constant_selected = force_constants_dict[bond_type_in]
    equ_dists_selected = equ_dists_dict[bond_type_in]

    return [(d - equ_dists_selected) * force_constant_selected for d in distances_in]

def parse_name(name):
    pattern = re.compile(r"^(?P<Base>[A-Za-z0-9]+)_(?P<Force>[0-9.]+)nN_R(?P<Run>\d+)$")
    m = pattern.match(name)
    if m:
        return m.group("Base"), float(m.group("Force")), int(m.group("Run"))
    else:
        return None, None, None

def plot_forces(samples, strands, forces_in, bond_types, substract, graphics_dir_in, plot_name):
    for sample in samples:
        for strand in strands:
            for cur_force in forces_in:
                for bond_type in bond_types:
                    distances_collect_singles = []
                    for run in [1, 2, 3]:
                        subset = df[
                            (df["Base"] == sample)
                            & (df["Force"] == cur_force)
                            & (df["Run"] == run)
                            & (df["Strand"] == strand)
                            & (df["BondType"] == bond_type)
                            ].copy()

                        selection_distances = subset["AverageDistance"].to_numpy().tolist()
                        distances_collect_singles.append(selection_distances)
                        selection_residues = subset["ResidueIndex"].to_numpy().tolist()

                    distances_run_average = [(a + b + c) / 3 for a, b, c in
                                             zip(distances_collect_singles[0], distances_collect_singles[1],
                                                 distances_collect_singles[2])]

                    forces = convert_to_forces(distances_run_average, bond_type)

                    if substract:
                        distances_collect = []
                        for run in [1, 2, 3]:
                            subset = df[
                                (df["Base"] == sample)
                                & (df["Force"] == 0.1)
                                & (df["Run"] == run)
                                & (df["Strand"] == strand)
                                & (df["BondType"] == bond_type)
                                ].copy()
                            selection_residues = subset["ResidueIndex"].to_numpy().tolist()
                            selection_distances = subset["AverageDistance"].to_numpy().tolist()
                            distances_collect.append(selection_distances)
                        distances_run_average = [(a + b + c) / 3 for a, b, c in
                                                 zip(distances_collect[0], distances_collect[1], distances_collect[2])]
                        forces_small = convert_to_forces(distances_run_average, bond_type)
                        forces = [f_1 - f_2 for f_1, f_2 in zip(forces, forces_small)]

                    plt.plot(selection_residues, forces, label=bond_type)

                plt.xticks([10, 30, 50, 70, 90], ["-40", "-20", "0", "+20", "40"])

                plt.xlabel("Relative residue index")
                plt.ylabel("Force [nN]")
                plt.tight_layout()
                plt.legend(fontsize=6, ncol=3)
                plt.savefig(f"{graphics_dir_in}/{sample}_{cur_force}nN_{plot_name}.pdf")
                plt.clf()

def plot_strand_comparison(samples, strand_labels, zero_pos, break_length, graphics_dir_in):
    bond_type = "C3'-O3'"
    for sample in samples:
        for strand in ["NHP", "HP"]:
            distances_collect_singles = []
            for run in [1, 2, 3]:
                subset = df[
                    (df["Base"] == sample)
                    & (df["Force"] == 2.0)
                    & (df["Run"] == run)
                    & (df["Strand"] == strand)
                    & (df["BondType"] == bond_type)
                    ].copy()

                selection_distances = subset["AverageDistance"].to_numpy().tolist()
                distances_collect_singles.append(selection_distances)
            distances_run_average = [(a + b + c) / 3 for a, b, c in
                                     zip(distances_collect_singles[0], distances_collect_singles[1],
                                         distances_collect_singles[2])]
            forces = convert_to_forces(distances_run_average, bond_type)
            forces_1 = convert_to_forces(distances_collect_singles[0], bond_type)
            forces_2 = convert_to_forces(distances_collect_singles[1], bond_type)
            forces_3 = convert_to_forces(distances_collect_singles[2], bond_type)

            distances_collect = []
            for run in [1, 2, 3]:
                subset = df[
                    (df["Base"] == sample)
                    & (df["Force"] == 0.1)
                    & (df["Run"] == run)
                    & (df["Strand"] == strand)
                    & (df["BondType"] == bond_type)
                    ].copy()
                selection_residues = subset["ResidueIndex"].to_numpy().tolist()
                selection_distances = subset["AverageDistance"].to_numpy().tolist()
                distances_collect.append(selection_distances)
            distances_run_average = [(a + b + c) / 3 for a, b, c in
                                     zip(distances_collect[0], distances_collect[1], distances_collect[2])]
            forces_small = convert_to_forces(distances_run_average, bond_type)

            forces_corrected = [f_1 - f_2 for f_1, f_2 in zip(forces, forces_small)]
            forces_corrected_1 = [f_1 - f_2 for f_1, f_2 in zip(forces_1, forces_small)]
            forces_corrected_2 = [f_1 - f_2 for f_1, f_2 in zip(forces_2, forces_small)]
            forces_corrected_3 = [f_1 - f_2 for f_1, f_2 in zip(forces_3, forces_small)]

            if strand == "NHP":
                color = purple
                selection_residues = [-(x-zero_pos[1]) for x in selection_residues]
                label = strand_labels[1]
            else:
                label = strand_labels[0]
                selection_residues = [x - zero_pos[0] for x in selection_residues]
                color = teal

            if break_length and strand == "NHP":
                selection_residues_fixed_1 = []
                selection_residues_fixed_2 = []
                for x in selection_residues:
                    if x < 0:
                        selection_residues_fixed_1.append(x - break_length)
                    else:
                        selection_residues_fixed_2.append(x + break_length)
                plt.plot(selection_residues_fixed_2, forces_corrected[:len(selection_residues_fixed_1)], color=color, label=label)
                plt.plot(selection_residues_fixed_2, forces_corrected_1[:len(selection_residues_fixed_1)], color=color, alpha=0.3)
                plt.plot(selection_residues_fixed_2, forces_corrected_2[:len(selection_residues_fixed_1)], color=color, alpha=0.3)
                plt.plot(selection_residues_fixed_2, forces_corrected_3[:len(selection_residues_fixed_1)], color=color, alpha=0.3)
                plt.plot(selection_residues_fixed_1, forces_corrected[len(selection_residues_fixed_1):], color=color)
                plt.plot(selection_residues_fixed_1, forces_corrected_1[len(selection_residues_fixed_1):], color=color, alpha=0.3)
                plt.plot(selection_residues_fixed_1, forces_corrected_2[len(selection_residues_fixed_1):], color=color, alpha=0.3)
                plt.plot(selection_residues_fixed_1, forces_corrected_3[len(selection_residues_fixed_1):], color=color, alpha=0.3)
            else:
                plt.plot(selection_residues, forces_corrected, color=color, label=label)
                plt.plot(selection_residues, forces_corrected_1, color=color, alpha=0.3)
                plt.plot(selection_residues, forces_corrected_2, color=color, alpha=0.3)
                plt.plot(selection_residues, forces_corrected_3, color=color, alpha=0.3)

    if break_length:
        plt.axvspan(-break_length-0.5, break_length+0.5, facecolor='grey', alpha=0.10, zorder=0)

    plt.xticks([-40, -20, 0, 20, 40], ["-40", "-20", "0", "+20", "+40"])

    plt.xlabel("Relative residue index")
    plt.ylabel("Force [nN]")
    plt.tight_layout()
    plt.ylim(-0.3, 2.3)
    plt.legend()
    plt.savefig(f"{graphics_dir_in}/{sample}_Strand_Comparision.pdf")
    plt.clf()

def plot_sample_comparison(samples, strands, sample_translations, zero_pos, plot_name, graphics_dir_in):
    for sample in samples:
        for strand in strands:
                bond_type = "C3'-O3'"
                distances_collect_singles = []
                for run in [1, 2, 3]:
                    subset = df[
                        (df["Base"] == sample)
                        & (df["Force"] == 2.0)
                        & (df["Run"] == run)
                        & (df["Strand"] == strand)
                        & (df["BondType"] == bond_type)
                        ].copy()

                    selection_distances = subset["AverageDistance"].to_numpy().tolist()
                    distances_collect_singles.append(selection_distances)

                distances_run_average = [(a + b + c) / 3 for a, b, c in
                                         zip(distances_collect_singles[0], distances_collect_singles[1],
                                             distances_collect_singles[2])]
                forces = convert_to_forces(distances_run_average, bond_type)

                distances_collect = []
                for run in [1, 2, 3]:
                    subset = df[
                        (df["Base"] == sample)
                        & (df["Force"] == 0.1)
                        & (df["Run"] == run)
                        & (df["Strand"] == strand)
                        & (df["BondType"] == bond_type)
                        ].copy()
                    selection_residues = subset["ResidueIndex"].to_numpy().tolist()
                    selection_distances = subset["AverageDistance"].to_numpy().tolist()
                    distances_collect.append(selection_distances)
                distances_run_average = [(a + b + c) / 3 for a, b, c in
                                         zip(distances_collect[0], distances_collect[1], distances_collect[2])]
                forces_small = convert_to_forces(distances_run_average, bond_type)

                forces_corrected = [f_1 - f_2 for f_1, f_2 in zip(forces, forces_small)]
                selection_residues = [a-zero_pos[sample] for a in selection_residues]

                plt.plot(selection_residues, forces_corrected,
                         color=colors[sample], label=sample_translations[sample])

    plt.xticks([-40, -20, 0, 20, 40], ["-40", "-20", "0", "+20", "+40"])

    plt.xlabel("Relative residue index")
    plt.ylabel("Force [nN]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{graphics_dir_in}/{plot_name}.pdf")
    plt.clf()

def plot_force_variation(samples, strands, zero_pos, graphics_dir_in):
    for sample in samples:
        for strand in strands:
            zero_pos_cur = zero_pos[strand]
            for force in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
                bond_type = "C3'-O3'"
                distances_collect = []
                for run in [1, 2, 3]:
                    subset = df[
                        (df["Base"] == sample)
                        & (df["Force"] == force)
                        & (df["Run"] == run)
                        & (df["Strand"] == strand)
                        & (df["BondType"] == bond_type)
                        ].copy()

                    selection_distances = subset["AverageDistance"].to_numpy().tolist()
                    distances_collect.append(selection_distances)
                distances_run_average = [(a + b + c) / 3 for a, b, c in
                                         zip(distances_collect[0], distances_collect[1], distances_collect[2])]
                forces = convert_to_forces(distances_run_average, bond_type)

                distances_collect = []
                for run in [1, 2, 3]:
                    subset = df[
                        (df["Base"] == sample)
                        & (df["Force"] == 0.1)
                        & (df["Run"] == run)
                        & (df["Strand"] == strand)
                        & (df["BondType"] == bond_type)
                        ].copy()
                    selection_residues = subset["ResidueIndex"].to_numpy().tolist()
                    selection_distances = subset["AverageDistance"].to_numpy().tolist()
                    distances_collect.append(selection_distances)
                distances_run_average = [(a + b + c) / 3 for a, b, c in
                                         zip(distances_collect[0], distances_collect[1], distances_collect[2])]
                forces_small = convert_to_forces(distances_run_average, bond_type)

                forces_corrected = [f_1 - f_2 for f_1, f_2 in zip(forces, forces_small)]
                selection_residues = [a - zero_pos_cur for a in selection_residues]

                plt.plot(selection_residues, forces_corrected, label=f"{force} nN")

            plt.xticks([-40, -20, 0, 20, 40], ["-40", "-20", "0", "+20", "+40"])

            plt.legend()
            plt.xlabel("Relative residue index")
            plt.ylabel("Force [nN]")
            plt.tight_layout()
            plt.savefig(f"{graphics_dir_in}/{sample}_{strand}_ForceVariation.pdf")
            plt.clf()

apply_plot_config(wide_plot_config)
purple = (0.580, 0.404, 0.741)
teal = (0.149, 0.651, 0.604)

CSV_FILE = "../ExtractedData/backbone_bonds_all.csv"
GRAPHICS_DIR = "../Graphics/"
df = pd.read_csv(CSV_FILE)

df[["Base", "Force", "Run"]] = df["SampleName"].apply(lambda x: pd.Series(parse_name(x)))
df["AverageDistance"] = pd.to_numeric(df["AverageDistance"], errors="coerce")
df["ResidueIndex"] = pd.to_numeric(df["ResidueIndex"], errors="coerce")

colors = {
           "HP0104": (0.192, 0.627, 0.204, 1.0),
           "HP064": (0.580, 0.404, 0.741, 1.0),
           "HP044": (0.690196, 0.0, 0.125490, 1.0),
           "HP004": (1.000, 0.498, 0.0, 1.0),
           "HP0104N": (0.192, 0.627, 0.204, 1.0),
           "HP064N": (0.580, 0.404, 0.741, 1.0),
           "HP044N": (0.690196, 0.0, 0.125490, 1.0),
           "HP004N": (1.000, 0.498, 0.0, 1.0),
        }

all_samples = ["HP404", "HP004", "HP044", "HP064", "HP0104"]
samples_no_control_step1 = ["HP004", "HP044", "HP064", "HP0104"]
samples_no_control_step2 = ["HP004N", "HP044N", "HP064N", "HP0104N"]
zero_positions_hp = {"HP004N": 52.5, "HP044N": 56.5, "HP064N": 58.5, "HP0104N": 62.5, "HP004": 52.5, "HP044": 56.5, "HP064": 58.5, "HP0104": 62.5}
zero_positions_nhp = {"HP004": 50.5, "HP044": 50.5, "HP064": 50.5, "HP0104": 50.5, "HP004N": 50.5, "HP044N": 50.5, "HP064N": 50.5, "HP0104N": 50.5}


sample_trans = {"HP004": "Hp(0,0,4)", "HP044": "Hp(0,4,4)", "HP064": "Hp(0,6,4)", "HP0104": "Hp(0,10,4)",
                "HP004N": "Hp(0,0,4)", "HP044N": "Hp(0,4,4)", "HP064N": "Hp(0,6,4)", "HP0104N": "Hp(0,10,4)"}

apply_plot_config(wide_plot_config)

# Bond Forces (without substraction)
plot_forces(all_samples, ["NHP"], [0.1, 2.0], ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"], False, graphics_dir_in=GRAPHICS_DIR, plot_name="All_Bonds")

# With substraction
plot_forces(all_samples, ["NHP"], [2.0], ["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"], True, graphics_dir_in=GRAPHICS_DIR, plot_name="All_Bonds_Substracted")

# Bubble First Step
plot_strand_comparison(["HP404"], ["Upper Strand", "Lower Strand"], [52.5, 52.5], None,  GRAPHICS_DIR)

# HP064 First + Second Step
plot_strand_comparison(["HP064"], ["Hairpin Strand", "Non-Hairpin Strand"], [58.5, 50.5], 8, GRAPHICS_DIR)
plot_strand_comparison(["HP064N"], ["Hairpin Strand", "Non-Hairpin Strand"], [58.5, 50.5], 8, GRAPHICS_DIR)

# Variation of Forces
plot_force_variation(["HP064", "HP064N"], ["NHP", "HP"], {"HP": 58.5, "NHP": 50.5}, GRAPHICS_DIR)

# Stem Variation Comparision
plot_sample_comparison(samples_no_control_step1, ["HP"], sample_trans, zero_positions_hp,"Stem_Comparison_Step1_HP", GRAPHICS_DIR)
plot_sample_comparison(samples_no_control_step1, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Comparison_Step1_NHP", GRAPHICS_DIR)

plot_sample_comparison(samples_no_control_step2, ["HP"], sample_trans, zero_positions_hp, "Stem_Comparison_Step2_HP", GRAPHICS_DIR)
plot_sample_comparison(samples_no_control_step2, ["NHP"], sample_trans, zero_positions_nhp,"Stem_Comparison_Step2_NHP", GRAPHICS_DIR)
