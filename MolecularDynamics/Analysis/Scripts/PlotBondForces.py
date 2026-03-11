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
    pattern = re.compile(r"^(?P<Base>[^_]+)_(?P<Force>\d+(?:\.\d+)?)nN_R(?P<Run>\d+)$")
    m = pattern.match(str(name))
    if m:
        return m.group("Base"), float(m.group("Force")), int(m.group("Run"))
    else:
        return None, None, None


def get_available_runs(df, sample, force, strand, bond_type):
    subset = df[
        (df["Base"] == sample)
        & (df["Force"] == force)
        & (df["Strand"] == strand)
        & (df["BondType"] == bond_type)
    ]
    runs = subset["Run"].dropna().astype(int).unique().tolist()
    return runs


def get_run_distance_lists(df, sample, force, strand, bond_type):
    runs = get_available_runs(df, sample, force, strand, bond_type)

    run_data = []
    for run in runs:
        subset = df[
            (df["Base"] == sample)
            & (df["Force"] == force)
            & (df["Run"] == run)
            & (df["Strand"] == strand)
            & (df["BondType"] == bond_type)
        ].copy()

        residues = subset["ResidueIndex"].to_numpy().tolist()
        distances = subset["AverageDistance"].to_numpy().tolist()

        if len(distances) == 0:
            continue

        run_data.append((run, residues, distances))

    return run_data


def get_mean_distances(df, sample, force, strand, bond_type):
    run_data = get_run_distance_lists(df, sample, force, strand, bond_type)

    if not run_data:
        return [], []

    min_len = min(len(distances) for _, _, distances in run_data)

    residues = run_data[0][1][:min_len]
    distance_lists = [distances[:min_len] for _, _, distances in run_data]

    mean_distances = [sum(vals) / len(vals) for vals in zip(*distance_lists)]

    return residues, mean_distances


def plot_forces(samples, strands, forces_in, bond_types, substract, graphics_dir_in, plot_name):
    for sample in samples:
        for strand in strands:
            for cur_force in forces_in:
                plotted_any = False

                for bond_type in bond_types:
                    selection_residues, distances_run_average = get_mean_distances(
                        df, sample, cur_force, strand, bond_type
                    )

                    if not distances_run_average:
                        continue

                    forces = convert_to_forces(distances_run_average, bond_type)

                    if substract:
                        _, distances_run_average_small = get_mean_distances(
                            df, sample, 0.1, strand, bond_type
                        )

                        if not distances_run_average_small:
                            continue

                        forces_small = convert_to_forces(distances_run_average_small, bond_type)

                        min_len = min(len(selection_residues), len(forces), len(forces_small))
                        selection_residues = selection_residues[:min_len]
                        forces = [f_1 - f_2 for f_1, f_2 in zip(forces[:min_len], forces_small[:min_len])]

                    plt.plot(selection_residues, forces, label=bond_type)
                    plotted_any = True

                if plotted_any:
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
        plotted_any = False

        for strand in ["NHP", "HP"]:
            run_data = get_run_distance_lists(df, sample, 2.0, strand, bond_type)
            if not run_data:
                continue

            min_len = min(len(dists) for _, _, dists in run_data)
            selection_residues = run_data[0][1][:min_len]

            distance_lists = [dists[:min_len] for _, _, dists in run_data]
            distances_run_average = [sum(vals) / len(vals) for vals in zip(*distance_lists)]

            forces = convert_to_forces(distances_run_average, bond_type)
            forces_runs = [convert_to_forces(dists, bond_type) for dists in distance_lists]

            small_run_data = get_run_distance_lists(df, sample, 0.1, strand, bond_type)
            if not small_run_data:
                continue

            min_len_small = min(len(dists) for _, _, dists in small_run_data)
            distance_lists_small = [dists[:min_len_small] for _, _, dists in small_run_data]
            distances_run_average_small = [sum(vals) / len(vals) for vals in zip(*distance_lists_small)]
            forces_small = convert_to_forces(distances_run_average_small, bond_type)

            min_len_final = min(len(selection_residues), len(forces), len(forces_small))
            selection_residues = selection_residues[:min_len_final]
            forces = forces[:min_len_final]
            forces_small = forces_small[:min_len_final]

            forces_corrected = [f_1 - f_2 for f_1, f_2 in zip(forces, forces_small)]
            forces_runs_corrected = [
                [f_1 - f_2 for f_1, f_2 in zip(run_forces[:min_len_final], forces_small)]
                for run_forces in forces_runs
            ]

            if strand == "NHP":
                color = purple
                selection_residues = [-(x - zero_pos[1]) for x in selection_residues]
                label = strand_labels[1]
            else:
                color = teal
                selection_residues = [x - zero_pos[0] for x in selection_residues]
                label = strand_labels[0]

            if break_length and strand == "NHP":
                selection_residues_fixed_1 = []
                selection_residues_fixed_2 = []

                for x in selection_residues:
                    if x < 0:
                        selection_residues_fixed_1.append(x - break_length)
                    else:
                        selection_residues_fixed_2.append(x + break_length)

                split_idx = len(selection_residues_fixed_1)

                plt.plot(selection_residues_fixed_2, forces_corrected[:split_idx], color=color, label=label)
                for run_forces in forces_runs_corrected:
                    plt.plot(selection_residues_fixed_2, run_forces[:split_idx], color=color, alpha=0.3)

                plt.plot(selection_residues_fixed_1, forces_corrected[split_idx:], color=color)
                for run_forces in forces_runs_corrected:
                    plt.plot(selection_residues_fixed_1, run_forces[split_idx:], color=color, alpha=0.3)
            else:
                plt.plot(selection_residues, forces_corrected, color=color, label=label)
                for run_forces in forces_runs_corrected:
                    plt.plot(selection_residues, run_forces, color=color, alpha=0.3)

            plotted_any = True

        if plotted_any:
            if break_length:
                plt.axvspan(-break_length - 0.5, break_length + 0.5, facecolor="grey", alpha=0.10, zorder=0)

            plt.xticks([-40, -20, 0, 20, 40], ["-40", "-20", "0", "+20", "+40"])
            plt.xlabel("Relative residue index")
            plt.ylabel("Force [nN]")
            plt.tight_layout()
            plt.ylim(-0.3, 2.3)
            plt.legend()
            plt.savefig(f"{graphics_dir_in}/{sample}_Strand_Comparision.pdf")
            plt.clf()


def plot_sample_comparison(samples, strands, sample_translations, zero_pos, plot_name, graphics_dir_in):
    plotted_any = False

    for sample in samples:
        for strand in strands:
            bond_type = "C3'-O3'"

            selection_residues, distances_run_average = get_mean_distances(
                df, sample, 2.0, strand, bond_type
            )
            if not distances_run_average:
                continue

            forces = convert_to_forces(distances_run_average, bond_type)

            _, distances_run_average_small = get_mean_distances(
                df, sample, 0.1, strand, bond_type
            )
            if not distances_run_average_small:
                continue

            forces_small = convert_to_forces(distances_run_average_small, bond_type)

            min_len = min(len(selection_residues), len(forces), len(forces_small))
            selection_residues = selection_residues[:min_len]
            forces_corrected = [f_1 - f_2 for f_1, f_2 in zip(forces[:min_len], forces_small[:min_len])]
            selection_residues = [a - zero_pos[sample] for a in selection_residues]

            plt.plot(
                selection_residues,
                forces_corrected,
                color=colors[sample],
                label=sample_translations[sample]
            )
            plotted_any = True

    if plotted_any:
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
            plotted_any = False

            for force in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
                bond_type = "C3'-O3'"

                selection_residues, distances_run_average = get_mean_distances(
                    df, sample, force, strand, bond_type
                )
                if not distances_run_average:
                    continue

                forces = convert_to_forces(distances_run_average, bond_type)

                _, distances_run_average_small = get_mean_distances(
                    df, sample, 0.1, strand, bond_type
                )
                if not distances_run_average_small:
                    continue

                forces_small = convert_to_forces(distances_run_average_small, bond_type)

                min_len = min(len(selection_residues), len(forces), len(forces_small))
                selection_residues_cur = [a - zero_pos_cur for a in selection_residues[:min_len]]
                forces_corrected = [f_1 - f_2 for f_1, f_2 in zip(forces[:min_len], forces_small[:min_len])]

                plt.plot(selection_residues_cur, forces_corrected, label=f"{force} nN")
                plotted_any = True

            if plotted_any:
                plt.xticks([-40, -20, 0, 20, 40], ["-40", "-20", "0", "+20", "+40"])
                plt.legend()
                plt.xlabel("Relative residue index")
                plt.ylabel("Force [nN]")
                plt.tight_layout()
                plt.savefig(f"{graphics_dir_in}/{sample}_{strand}_ForceVariation.pdf")
                plt.clf()

def get_mean_and_std_distances(df, sample, force, strand, bond_type):
    run_data = get_run_distance_lists(df, sample, force, strand, bond_type)

    if not run_data:
        return [], [], []

    min_len = min(len(distances) for _, _, distances in run_data)

    residues = run_data[0][1][:min_len]
    distance_lists = [distances[:min_len] for _, _, distances in run_data]

    mean_distances = [sum(vals) / len(vals) for vals in zip(*distance_lists)]

    if len(distance_lists) > 1:
        std_distances = pd.DataFrame(distance_lists).std(axis=0, ddof=1).tolist()
    else:
        std_distances = [0.0] * min_len

    return residues, mean_distances, std_distances


def plot_sample_comparison_with_std(samples, strands, sample_translations, zero_pos, plot_name, graphics_dir_in):
    plotted_any = False

    for sample in samples:
        for strand in strands:
            bond_type = "C3'-O3'"

            selection_residues, mean_distances, std_distances = get_mean_and_std_distances(
                df, sample, 2.0, strand, bond_type
            )
            if not mean_distances:
                continue

            _, mean_distances_small, std_distances_small = get_mean_and_std_distances(
                df, sample, 0.1, strand, bond_type
            )
            if not mean_distances_small:
                continue

            forces_mean = convert_to_forces(mean_distances, bond_type)
            forces_small_mean = convert_to_forces(mean_distances_small, bond_type)

            forces_plus = convert_to_forces(
                [m + s for m, s in zip(mean_distances, std_distances)], bond_type
            )
            forces_minus = convert_to_forces(
                [m - s for m, s in zip(mean_distances, std_distances)], bond_type
            )

            forces_small_plus = convert_to_forces(
                [m + s for m, s in zip(mean_distances_small, std_distances_small)], bond_type
            )
            forces_small_minus = convert_to_forces(
                [m - s for m, s in zip(mean_distances_small, std_distances_small)], bond_type
            )

            min_len = min(
                len(selection_residues),
                len(forces_mean),
                len(forces_small_mean),
                len(forces_plus),
                len(forces_minus)
            )

            selection_residues = selection_residues[:min_len]
            selection_residues = [a - zero_pos[sample] for a in selection_residues]

            forces_corrected_mean = [
                f1 - f2 for f1, f2 in zip(forces_mean[:min_len], forces_small_mean[:min_len])
            ]

            forces_corrected_plus = [
                f1 - f2 for f1, f2 in zip(forces_plus[:min_len], forces_small_plus[:min_len])
            ]

            forces_corrected_minus = [
                f1 - f2 for f1, f2 in zip(forces_minus[:min_len], forces_small_minus[:min_len])
            ]

            color = colors[sample]

            # mean line
            plt.plot(
                selection_residues,
                forces_corrected_mean,
                color=color,
                label=sample_translations[sample]
            )

            # std "tube"
            plt.fill_between(
                selection_residues,
                forces_corrected_minus,
                forces_corrected_plus,
                color=color,
                alpha=0.25
            )

            plotted_any = True

    if plotted_any:
        plt.xticks([-40, -20, 0, 20, 40], ["-40", "-20", "0", "+20", "+40"])
        plt.xlabel("Relative residue index")
        plt.ylabel("Force [nN]")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{graphics_dir_in}/{plot_name}.pdf")
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
            "HP034AT": (0.498, 0.114, 0.114, 1.0),
            "HP034": (0.090, 0.090, 0.090, 1.0),
            "HP044AT": (0.086, 0.639, 0.290, 1.0),
            "HP054": (0.0, 0.898, 1.0, 1.0),
            "HP054AT": (0.737, 0.620, 0.404, 1.0),
            "HP064AT": (0.000, 0.569, 1.000, 1.0),
        }

all_samples = ["HP404", "HP004", "HP044", "HP064", "HP0104"]
all_samples_with_stem_type = ["HP404", "HP004", "HP044", "HP064", "HP0104", "HP034", "HP034AT", "HP044AT", "HP054", "HP054AT", "HP064AT"]
samples_no_control_step1 = ["HP004", "HP044", "HP064", "HP0104"]
samples_no_control_step2 = ["HP004N", "HP044N", "HP064N", "HP0104N"]

zero_positions_hp = {
    "HP004N": 52.5, "HP044N": 56.5, "HP064N": 58.5, "HP0104N": 62.5,
    "HP004": 52.5, "HP044": 56.5, "HP064": 58.5, "HP0104": 62.5
}
zero_positions_nhp = {
    "HP004": 50.5, "HP044": 50.5, "HP064": 50.5, "HP0104": 50.5,
    "HP004N": 50.5, "HP044N": 50.5, "HP064N": 50.5, "HP0104N": 50.5,
    "HP034": 50.5, "HP034AT": 50.5, "HP044AT": 50.5, "HP054": 50.5,
    "HP054AT": 50.5, "HP064AT": 50.5
}

stem_type_variation_samples_GC = ["HP004", "HP034", "HP044", "HP054", "HP064"]
stem_type_variation_samples_AT = ["HP004", "HP034AT", "HP044AT", "HP054AT", "HP064AT"]

sample_trans = {
    "HP004": "Hp(0,0,4)",
    "HP044": "Hp(0,4,4)",
    "HP064": "Hp(0,6,4)",
    "HP0104": "Hp(0,10,4)",
    "HP004N": "Hp(0,0,4)",
    "HP044N": "Hp(0,4,4)",
    "HP064N": "Hp(0,6,4)",
    "HP0104N": "Hp(0,10,4)",
    "HP034": "Hp(0,3,4)",
    "HP034AT": "Hp(0,3,4) AT",
    "HP044AT": "Hp(0,4,4) AT",
    "HP054": "Hp(0,5,4)",
    "HP054AT": "Hp(0,5,4) AT",
    "HP064AT": "Hp(0,6,4) AT"
}

apply_plot_config(wide_plot_config)

# Bond Forces (without subtraction)
plot_forces(all_samples_with_stem_type,["NHP"],[0.1, 2.0],["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"],False,graphics_dir_in=GRAPHICS_DIR,plot_name="All_Bonds")

# With subtraction
plot_forces(all_samples_with_stem_type,["NHP"],[2.0],["P-O5'", "O5'-C5'", "C5'-C4'", "C4'-C3'", "C3'-O3'", "O3'-P"],True,graphics_dir_in=GRAPHICS_DIR,plot_name="All_Bonds_Substracted")

# Bubble First Step
plot_strand_comparison(["HP404"], ["Upper Strand", "Lower Strand"], [52.5, 52.5], None, GRAPHICS_DIR)

# HP064 First + Second Step
plot_strand_comparison(["HP064"], ["Hairpin Strand", "Non-Hairpin Strand"], [58.5, 50.5], 8, GRAPHICS_DIR)
plot_strand_comparison(["HP064N"], ["Hairpin Strand", "Non-Hairpin Strand"], [58.5, 50.5], 8, GRAPHICS_DIR)

# Variation of Forces
plot_force_variation(["HP064", "HP064N"], ["NHP", "HP"], {"HP": 58.5, "NHP": 50.5}, GRAPHICS_DIR)

# Stem Variation Comparison
plot_sample_comparison(samples_no_control_step1, ["HP"], sample_trans, zero_positions_hp, "Stem_Comparison_Step1_HP", GRAPHICS_DIR)
plot_sample_comparison(samples_no_control_step1, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Comparison_Step1_NHP", GRAPHICS_DIR)
plot_sample_comparison_with_std(samples_no_control_step1,["NHP"],sample_trans,zero_positions_nhp,"Stem_Comparison_Step1_NHP_with_std",GRAPHICS_DIR)

plot_sample_comparison(samples_no_control_step2, ["HP"], sample_trans, zero_positions_hp, "Stem_Comparison_Step2_HP", GRAPHICS_DIR)
plot_sample_comparison(samples_no_control_step2, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Comparison_Step2_NHP", GRAPHICS_DIR)

# Stem Type Variation Comparison
plot_sample_comparison(stem_type_variation_samples_GC, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_GC", GRAPHICS_DIR)
plot_sample_comparison_with_std(stem_type_variation_samples_GC, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_GC_with_std", GRAPHICS_DIR)

plot_sample_comparison(stem_type_variation_samples_AT, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_AT", GRAPHICS_DIR)
plot_sample_comparison_with_std(stem_type_variation_samples_AT, ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_AT_with_std", GRAPHICS_DIR)

plot_sample_comparison(["HP034", "HP034AT"], ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_Pairwise3", GRAPHICS_DIR)
plot_sample_comparison_with_std(["HP034", "HP034AT"],["NHP"],sample_trans,zero_positions_nhp,"Stem_Type_Comparision_Pairwise3_with_std",GRAPHICS_DIR)
plot_sample_comparison(["HP044", "HP044AT"], ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_Pairwise4", GRAPHICS_DIR)
plot_sample_comparison_with_std(["HP044", "HP044AT"],["NHP"], sample_trans, zero_positions_nhp,"Stem_Type_Comparision_Pairwise4_with_std", GRAPHICS_DIR)
plot_sample_comparison(["HP054", "HP054AT"], ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_Pairwise5", GRAPHICS_DIR)
plot_sample_comparison_with_std(["HP054", "HP054AT"],["NHP"],sample_trans,zero_positions_nhp,"Stem_Type_Comparision_Pairwise5_with_std",GRAPHICS_DIR)
plot_sample_comparison(["HP064", "HP064AT"], ["NHP"], sample_trans, zero_positions_nhp, "Stem_Type_Comparision_Pairwise6", GRAPHICS_DIR)
plot_sample_comparison_with_std(["HP064", "HP064AT"],["NHP"],sample_trans,zero_positions_nhp,"Stem_Type_Comparision_Pairwise6_with_std",GRAPHICS_DIR)



