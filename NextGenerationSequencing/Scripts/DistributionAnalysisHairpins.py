from NGSPlotting import *
from NGSDataUtilities import *
from NGSFitting import *

import scipy.stats as stats


def remove_bars(data_in, bars_in):
    bar_1 = bars_in[0]
    bar_2 = bars_in[1]
    data_x = sorted(list(data_in.keys()))
    neighbour_hood_1 = [data_in[x] for x in data_x if bar_1 - 3 <= x <= bar_1 + 3 and x != bar_1]
    bar_1_replacement = sum(neighbour_hood_1) / len(neighbour_hood_1)
    neighbour_hood_2 = [data_in[x] for x in data_x if bar_2 - 3 <= x <= bar_2 + 3 and x != bar_2]
    bar_2_replacement = sum(neighbour_hood_2) / len(neighbour_hood_2)
    data_in[bar_1] = bar_1_replacement
    data_in[bar_2] = bar_2_replacement
    total_sum = sum(list(data_in.values()))
    for x_in in data_x:
        data_in[x_in] = data_in[x_in] / total_sum
    return data_in


INCLUDE_CONTROL = False
sample_groups = ["HP4.0.4", "HP0.0.4", "HP0.4.4", "HP0.6.4", "HP0.10.4", "HP0.6.0", "HP0.6.8", "HP4.6.4", "HP5.6.4", "HP0.6.4.GC"]

colors = {"HP0.10.4": (0.192, 0.627, 0.204, 1.0),
          "HP0.6.4": (0.580, 0.404, 0.741, 1.0),
          "HP0.4.4": (0.690196, 0.0, 0.125490, 1.0),
          "HP0.0.4": (1.000, 0.498, 0.0, 1.0),
          "HP0.6.0": (1.000, 0.498, 0.0, 1.0),
          "HP0.6.8": (1.000, 0.498, 0.0, 1.0),
          "HP4.6.4": (1.000, 0.498, 0.0, 1.0),
          "HP5.6.4": (1.000, 0.498, 0.0, 1.0),
          "HP0.6.4.GC": (1.000, 1.000, 1.0, 1.0)
          }

sample_name_dict = {
    "HP0.6.4": "Hp(0,6,4)",
    "HP064": "Hp(0,6,4)",

    "HP0.10.4": "Hp(0,10,4)",
    "HP0104": "Hp(0,10,4)",

    "HP0.4.4": "Hp(0,4,4)",
    "HP044": "Hp(0,4,4)",

    "HP0.0.4": "Hp(0,0,4)",
    "HP004": "Hp(0,0,4)",

    "HP4.0.4": "Hp(4,0,4)",
    "HP404": "Hp(4,0,4)",

    "HP0.6.0": "Hp(0,6,0)",
    "HP060": "Hp(0,6,0)",

    "HP0.6.8": "Hp(0,6,8)",
    "HP068": "Hp(0,6,8)",

    "HP4.6.4": "Hp(4,6,4)",
    "HP464": "Hp(4,6,4)",

    "HP5.6.4": "Hp(5,6,4)",
    "HP564": "Hp(5,6,4)",

    "HP0.6.4.GC": "Hp(0,6,4) mod",
    "HP064GC": "Hp(0,6,4) mod"
}

runs = [1, 2, 3]

proc_data_dir = "../ProcessedData/Hairpins/"
graphics_dir = "../Graphics/Hairpins/"

# Ligation positions
bars_top = {"HP0.0.4": [-28.5, 29.5], "HP0.4.4": [-32.5, 33.5], "HP0.10.4": [-38.5, 39.5],
            "HP4.0.4": [-28.5, 29.5], "HP0.6.4": [-34.5, 35.5], "HP0.6.0": [-32.5, 33.5],
            "HP0.6.8": [-36.5, 37.5], "HP4.6.4": [-34.5, 35.5], "HP5.6.4": [-34.5, 35.5], "HP0.6.4.GC": [-34.5, 35.5]}

bars_bot = {"HP0.0.4": [-30.5, 31.5], "HP0.4.4": [-30.5, 31.5], "HP0.10.4": [-30.5, 31.5],
            "HP4.0.4": [-32.5, 33.5], "HP0.6.4": [-30.5, 31.5], "HP0.6.0": [-30.5, 31.5],
            "HP0.6.8": [-30.5, 31.5], "HP4.6.4": [-32.5, 33.5], "HP5.6.4": [-33, 34], "HP0.6.4.GC": [-30.5, 31.5]}

widths_top = {}
centers_top = {}
widths_bot = {}
centers_bot = {}
plotter = PlottingNGS()

# Students t dist on all distributions
for sample_group in sample_groups:
    cur_group_width_top = []
    cur_group_width_bot = []
    cur_group_center_top = []
    cur_group_center_bot = []

    for run in runs:
        sample_name = f"{sample_group}_{run}"
        cur_breaking_data_norm_shifted = read_breaking_results([f"{sample_name}_Breaking_Dist_Norm_Shifted.csv"],
                                                               proc_data_dir)
        cur_breaking_data_norm_shifted_top = cur_breaking_data_norm_shifted[f"{sample_name}#TOP"]
        cur_breaking_data_norm_shifted_bot = cur_breaking_data_norm_shifted[f"{sample_name}#BOT"]
        cur_breaking_data_norm_shifted_top = remove_bars(cur_breaking_data_norm_shifted_top, bars_top[sample_group])
        cur_breaking_data_norm_shifted_bot = remove_bars(cur_breaking_data_norm_shifted_bot, bars_bot[sample_group])
        df_top, mu_top, sigma_top, IQR_top = fit_student(cur_breaking_data_norm_shifted_top)
        df_bot, mu_bot, sigma_bot, IQR_bot = fit_student(cur_breaking_data_norm_shifted_bot)
        rsquare_top = rsquares_student(cur_breaking_data_norm_shifted_top, df_top, mu_top, sigma_top)
        rsquare_bot = rsquares_student(cur_breaking_data_norm_shifted_bot, df_bot, mu_bot, sigma_bot)

        widths_bot[sample_name] = IQR_bot
        centers_bot[sample_name] = mu_bot
        widths_top[sample_name] = IQR_top
        centers_top[sample_name] = mu_top
        cur_group_width_top.append(IQR_top)
        cur_group_width_bot.append(IQR_bot)
        cur_group_center_top.append(mu_top)
        cur_group_center_bot.append(mu_bot)

    average_width_top = round(np.mean(cur_group_width_top), 2)
    average_width_bot = round(np.mean(cur_group_width_bot), 2)
    average_mu_top = round(np.mean(cur_group_center_top), 2)
    average_mu_bot = round(np.mean(cur_group_center_bot), 2)
    std_width_top = round(np.std(cur_group_width_top), 2)
    std_width_bot = round(np.std(cur_group_width_bot), 2)
    std_mu_top = round(np.std(cur_group_center_top), 2)
    std_mu_bot = round(np.std(cur_group_center_bot), 2)

    print(f"{sample_group}: Top: Mu: {average_mu_top} +/- {std_mu_top}, "
          f"IQR: {average_width_top} +/- {std_width_top}, Bot: Mu: {average_mu_bot} +/- {std_mu_bot}, "
          f"IQR: {average_width_bot} +/- {std_width_bot}")


sample_group_stem = ["HP0.0.4", "HP0.4.4", "HP0.6.4", "HP0.10.4"]
sample_group_loop = ["HP0.6.0", "HP0.6.4", "HP0.6.8"]
sample_group_unpaired = ["HP0.6.4", "HP4.6.4", "HP5.6.4"]
sample_group_sequence_change = ["HP0.0.4", "HP0.6.4", "HP0.6.4.GC"]
group_labels = ["Stem", "Loop", "Unpaired", "Sequence"]

colors_loop_top = [(0.1216, 0.3059, 0.8471, 1.0), (0.1490, 0.6510, 0.6039, 1.0), (1.0000, 0.8235, 0.2471, 1.0)]
colors_loop_bot = [(0.9569, 0.6510, 0.7569, 1.0), (0.580, 0.404, 0.741, 1.0), (0.5451, 0.3529, 0.1686, 1.0)]

colors_unpaired_top = [(0.1490, 0.6510, 0.6039, 1.0), (0.7608, 0.0941, 0.5451, 1.0), (1.0000, 0.4353, 0.3804, 1.0)]
colors_unpaired_bot = [(0.580, 0.404, 0.741, 1.0), (0.9020, 0.8275, 0.6392, 1.0), (0.7529, 0.7529, 0.7529, 1.0)]

colors_top = [{"HP0.0.4": (0.192, 0.627, 0.204, 1.0), "HP0.4.4": (0.690196, 0.0, 0.125490, 1.0), "HP0.6.4": (0.690196, 0.0, 0.125490, 1.0),"HP0.10.4": (0.580, 0.404, 0.741, 1.0)},
              {"HP0.6.0": (0.1216, 0.3059, 0.8471, 1.0), "HP0.6.4": (0.1490, 0.6510, 0.6039, 1.0), "HP0.6.8": (1.0000, 0.8235, 0.2471, 1.0)},
              {"HP0.6.4": (0.1490, 0.6510, 0.6039, 1.0), "HP4.6.4": (0.7608, 0.0941, 0.5451, 1.0), "HP5.6.4": (1.0000, 0.4353, 0.3804, 1.0)},
              {"HP0.0.4": (0.192, 0.627, 0.204, 1.0), "HP0.6.4": (0.1490, 0.6510, 0.6039, 1.0), "HP0.6.4.GC": (1.0000, 1.0000, 1.0000, 1.0)}]

colors_bot = [{"HP0.0.4": (0.192, 0.627, 0.204, 1.0), "HP0.4.4": (0.690196, 0.0, 0.125490, 1.0), "HP0.6.4": (0.690196, 0.0, 0.125490, 1.0),"HP0.10.4": (0.580, 0.404, 0.741, 1.0)},
              {"HP0.6.0": (0.9569, 0.6510, 0.7569, 1.0), "HP0.6.4": (0.580, 0.404, 0.741, 1.0), "HP0.6.8": (0.5451, 0.3529, 0.1686, 1.0)},
              {"HP0.6.4": (0.580, 0.404, 0.741, 1.0), "HP4.6.4": (0.9020, 0.8275, 0.6392, 1.0), "HP5.6.4": (0.7529, 0.7529, 0.7529, 1.0)},
              {"HP0.0.4": (0.192, 0.627, 0.204, 1.0), "HP0.6.4": (0.580, 0.404, 0.741, 1.0), "HP0.6.4.GC": (1.0000, 1.0000, 1.0000, 1.0)}]

for sample_group_variation, group_label, color_top, color_bot in zip([sample_group_stem, sample_group_loop, sample_group_unpaired, sample_group_sequence_change], group_labels, colors_top, colors_bot):
    # Statistics per variational group
    out_matrix_exact_mu_top = np.zeros((len(sample_group_variation), len(sample_group_variation)))
    out_matrix_exact_width_top = np.zeros((len(sample_group_variation), len(sample_group_variation)))
    out_matrix_exact_mu_bot = np.zeros((len(sample_group_variation), len(sample_group_variation)))
    out_matrix_exact_width_bot = np.zeros((len(sample_group_variation), len(sample_group_variation)))

    for i in range(0, len(sample_group_variation)):
        for j in range(0, len(sample_group_variation)):
            group_a = [centers_top[a] for a in centers_top.keys() if a.split("_")[0] == sample_group_variation[i]]
            group_b = [centers_top[a] for a in centers_top.keys() if a.split("_")[0] == sample_group_variation[j]]
            out_matrix_exact_mu_top[i][j] = stats.f_oneway(group_a, group_b).pvalue
            group_a = [widths_top[a] for a in widths_top.keys() if a.split("_")[0] == sample_group_variation[i]]
            group_b = [widths_top[a] for a in widths_top.keys() if a.split("_")[0] == sample_group_variation[j]]
            out_matrix_exact_width_top[i][j] = stats.f_oneway(group_a, group_b).pvalue
            group_a = [centers_bot[a] for a in centers_bot.keys() if a.split("_")[0] == sample_group_variation[i]]
            group_b = [centers_bot[a] for a in centers_bot.keys() if a.split("_")[0] == sample_group_variation[j]]
            out_matrix_exact_mu_bot[i][j] = stats.f_oneway(group_a, group_b).pvalue
            group_a = [widths_bot[a] for a in widths_bot.keys() if a.split("_")[0] == sample_group_variation[i]]
            group_b = [widths_bot[a] for a in widths_bot.keys() if a.split("_")[0] == sample_group_variation[j]]
            out_matrix_exact_width_bot[i][j] = stats.f_oneway(group_a, group_b).pvalue

        print_p_matrix(out_matrix_exact_width_top, sample_group_variation, f"{group_label}_Pvalues_Width_Top.csv", proc_data_dir)
        print_p_matrix(out_matrix_exact_width_bot, sample_group_variation, f"{group_label}_Pvalues_Width_Bot.csv", proc_data_dir)
        print_p_matrix(out_matrix_exact_mu_top, sample_group_variation, f"{group_label}_Pvalues_Mu_Top.csv", proc_data_dir)
        print_p_matrix(out_matrix_exact_mu_bot, sample_group_variation, f"{group_label}_Pvalues_Mu_Bot.csv", proc_data_dir)

        plotter.plot_p_value_matrix(f"{group_label}_Pvalues_Width_Top.csv", proc_data_dir,
                                    save_path=f"{graphics_dir}/{group_label}_Pvalues_Width_Top.pdf")
        plotter.plot_p_value_matrix(f"{group_label}_Pvalues_Width_Bot.csv", proc_data_dir,
                                    save_path=f"{graphics_dir}/{group_label}_Pvalues_Width_Bot.pdf")
        plotter.plot_p_value_matrix(f"{group_label}_Pvalues_Mu_Top.csv", proc_data_dir,
                                    save_path=f"{graphics_dir}/{group_label}_Pvalues_Mu_Top.pdf")
        plotter.plot_p_value_matrix(f"{group_label}_Pvalues_Mu_Bot.csv", proc_data_dir,
                                    save_path=f"{graphics_dir}/{group_label}_Pvalues_Mu_Bot.pdf")



    cur_group_width_top = {}
    cur_group_width_bot = {}
    cur_group_mu_top = {}
    cur_group_mu_bot = {}
    group_dict = {}


    for samples_group_var in sample_group_variation:
        for run in runs:
            group_dict[f"{samples_group_var}_{run}"] = samples_group_var
        group_dict[f"{samples_group_var}_Average"] = samples_group_var
    for sample_group_var in sample_group_variation:
        cur_dict_width_top = {key: float(widths_top[key]) for key in widths_top.keys() if
                              key.split("_")[0] == sample_group_var}
        cur_group_width_top = {**cur_group_width_top, **cur_dict_width_top,
                               f"{sample_group_var}_Average": float(np.mean(list(cur_dict_width_top.values())))}
        cur_dict_width_bot = {key: float(widths_bot[key]) for key in widths_bot.keys() if
                              key.split("_")[0] == sample_group_var}
        cur_group_width_bot = {**cur_group_width_bot, **cur_dict_width_bot,
                               f"{sample_group_var}_Average": float(np.mean(list(cur_dict_width_bot.values())))}
        cur_dict_mu_top = {key: float(centers_top[key]) for key in centers_top.keys() if
                           key.split("_")[0] == sample_group_var}
        cur_group_mu_top = {**cur_group_mu_top, **cur_dict_mu_top,
                            f"{sample_group_var}_Average": float(np.mean(list(cur_dict_mu_top.values())))}
        cur_dict_mu_bot = {key: float(centers_bot[key]) for key in centers_bot.keys() if
                           key.split("_")[0] == sample_group_var}
        cur_group_mu_bot = {**cur_group_mu_bot, **cur_dict_mu_bot,
                            f"{sample_group_var}_Average": float(np.mean(list(cur_dict_mu_bot.values())))}
    print_fit_parameters(cur_group_width_top, group_dict, f"{group_label}_Width_Top.csv", proc_data_dir)
    print_fit_parameters(cur_group_width_bot, group_dict, f"{group_label}_Width_Bot.csv", proc_data_dir)
    print_fit_parameters(cur_group_mu_top, group_dict, f"{group_label}_Mu_Top.csv", proc_data_dir)
    print_fit_parameters(cur_group_mu_bot, group_dict, f"{group_label}_Mu_Bot.csv", proc_data_dir)

    plotter.plot_fit_paramameter_scatter_single([f"{group_label}_Width_Top.csv"], proc_data_dir, "IQR",
                                                save_path=f"{graphics_dir}/{group_label}_Width_Top.pdf", colors=color_top,
                                                sample_names=sample_name_dict)
    plotter.plot_fit_paramameter_scatter_single([f"{group_label}_Width_Bot.csv"], proc_data_dir, "IQR",
                                                save_path=f"{graphics_dir}/{group_label}_Width_Bot.pdf", colors=color_bot,
                                                sample_names=sample_name_dict)
    plotter.plot_fit_paramameter_scatter_single([f"{group_label}_Mu_Top.csv"], proc_data_dir, "Mu",
                                                save_path=f"{graphics_dir}/{group_label}_Mu_Top.pdf", colors=color_top,
                                                sample_names=sample_name_dict)
    plotter.plot_fit_paramameter_scatter_single([f"{group_label}_Mu_Bot.csv"], proc_data_dir, "Mu",
                                                save_path=f"{graphics_dir}/{group_label}_Mu_Bot.pdf", colors=color_bot,
                                                sample_names=sample_name_dict)
