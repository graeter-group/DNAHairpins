from NGSPlotting import *
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


plotter = PlottingNGS()

graphics_output_directors = "../Graphics/Hairpins_Selfassembly/"
bars_bot = [-30.5, 31.5]

plotter.plot_single_breaking_data_average(
    "HPSA0.6.4_Breaking_Dist_Norm_Shifted_Average.csv", "../ProcessedData/Hairpins_Selfassembly/",
    "TOP", "STUDENT T", save_path=f"{graphics_output_directors}/HPSA0_Breaking_Dist_Norm_Shifted_Average_Top.pdf")
plotter.plot_single_breaking_data_average(
    "HPSA0.6.4_Breaking_Dist_Norm_Shifted_Average.csv", "../ProcessedData/Hairpins_Selfassembly/",
    "BOT", "STUDENT T", save_path=f"{graphics_output_directors}/HPSA0_Breaking_Dist_Norm_Shifted_Average_Bot.pdf")

IQRs = {"SA_BOT": [], "SA_TOP": [], "LI_BOT": [], "LI_TOP": []}
mus = {"SA_BOT": [], "SA_TOP": [], "LI_BOT": [], "LI_TOP": []}
self_assembly_IQRs_top = []
self_assembly_IQRs_bot = []
self_assembly_mus_top = []
self_assembly_mus_bot = []

for strand in ["TOP", "BOT"]:
    for sample_idx in [1, 2]:
        cur_breaking_data_norm_shifted = read_breaking_results([f"HPSA0.6.4_{sample_idx}_Breaking_Dist_Norm_Shifted.csv"],
                                  "../ProcessedData/Hairpins_Selfassembly/")[f"HPSA0.6.4_{sample_idx}#{strand}"]
        _, mu, _, IQR = fit_student(cur_breaking_data_norm_shifted)
        IQRs[f"SA_{strand}"].append(IQR)
        mus[f"SA_{strand}"].append(mu)

for strand in ["TOP", "BOT"]:
    for sample_idx in [1, 2, 3]:
        cur_breaking_data_norm_shifted = read_breaking_results([f"HP0.6.4_{sample_idx}_Breaking_Dist_Norm_Shifted.csv"],
                                                               "../ProcessedData/Hairpins/")[
            f"HP0.6.4_{sample_idx}#{strand}"]
        cur_breaking_data_norm_shifted = remove_bars(cur_breaking_data_norm_shifted, bars_bot)

        _, mu, _, IQR = fit_student(cur_breaking_data_norm_shifted)
        IQRs[f"LI_{strand}"].append(IQR)
        mus[f"LI_{strand}"].append(mu)

IQR_means = {"SA_BOT": [], "SA_TOP": [], "LI_BOT": [], "LI_TOP": []}
mu_means = {"SA_BOT": [], "SA_TOP": [], "LI_BOT": [], "LI_TOP": []}
IQR_stds = {"SA_BOT": [], "SA_TOP": [], "LI_BOT": [], "LI_TOP": []}
mu_stds = {"SA_BOT": [], "SA_TOP": [], "LI_BOT": [], "LI_TOP": []}

for key in IQRs.keys():
    IQR_means[key] = np.mean(IQRs[key])
    mu_means[key] = np.mean(mus[key])
    IQR_stds[key] = np.std(IQRs[key])
    mu_stds[key] = np.std(mus[key])
    print(
        f"{key}: mu = {round(np.mean(mus[key]), 2)} +/- {round(np.std(mus[key]), 2)}, IQR = {round(np.mean(IQRs[key]), 2)} +/- {round(np.std(IQRs[key]), 2)}")

print("p value on IQRs (ligation vs. self-assembly, top)",
      round(stats.f_oneway(IQRs["SA_TOP"], IQRs["LI_TOP"]).pvalue, 3))
print("p value on mus (ligation vs. self-assembly, top)", round(stats.f_oneway(mus["SA_TOP"], mus["LI_TOP"]).pvalue, 3))
print("p value on IQRs (ligation vs. self-assembly, bot)",
      round(stats.f_oneway(IQRs["SA_BOT"], IQRs["LI_BOT"]).pvalue, 3))
print("p value on mus (ligation vs. self-assembly, bot)", round(stats.f_oneway(mus["SA_BOT"], mus["LI_BOT"]).pvalue, 3))
