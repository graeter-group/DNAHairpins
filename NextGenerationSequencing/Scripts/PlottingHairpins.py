from NGSPlotting import *

sample_groups = ["HP4.0.4", "HP0.0.4", "HP0.4.4", "HP0.6.4", "HP0.10.4", "HP0.6.0", "HP0.6.8", "HP4.6.4", "HP5.6.4", "HP0.6.4.GC"]

graphics_dir = "../Graphics/Hairpins/"
proc_data_dir = "../ProcessedData/Hairpins/"

runs = [1, 2, 3]
plotter = PlottingNGS()

# Ligation positions
bars_top = {"HP0.0.4": [-28.5, 29.5], "HP0.4.4": [-32.5, 33.5], "HP0.10.4": [-38.5, 39.5],
            "HP4.0.4": [-28.5, 29.5], "HP0.6.4": [-34.5, 35.5], "HP0.6.0": [-32.5, 33.5],
            "HP0.6.8": [-36.5, 37.5], "HP4.6.4": [-34.5, 35.5], "HP5.6.4": [-34.5, 35.5], "HP0.6.4.GC": [-34.5, 35.5]}

bars_bot = {"HP0.0.4": [-30.5, 31.5], "HP0.4.4": [-30.5, 31.5], "HP0.10.4": [-30.5, 31.5],
            "HP4.0.4": [-32.5, 33.5], "HP0.6.4": [-30.5, 31.5], "HP0.6.0": [-30.5, 31.5],
            "HP0.6.8": [-30.5, 31.5], "HP4.6.4": [-32.5, 33.5], "HP5.6.4": [-33, 34], "HP0.6.4.GC": [-30.5, 31.5]}

stem_samples = ["HP0.10.4", "HP0.6.4", "HP0.4.4", "HP0.0.4"]
loop_samples = ["HP0.6.0", "HP0.6.4", "HP0.6.8"]
unpaired_samples = ["HP0.6.4", "HP4.6.4", "HP5.6.4"]

colors = {"HP0.10.4": (0.192, 0.627, 0.204, 1.0),
          "HP0.6.4": (0.580, 0.404, 0.741, 1.0),
          "HP0.4.4": (0.690196, 0.0, 0.125490, 1.0),
          "HP0.0.4": (1.000, 0.498, 0.0, 1.0),
          "HP0.6.0": (0.1216, 0.3059, 0.8471, 1.0),
          "HP0.6.8": (1.0000, 0.8235, 0.2471, 1.0),
          "HP4.6.4": (1.0, 1.0, 1.0, 1.0),
          "HP5.6.4": (1.0, 1.0, 1.0, 1.0),
          "HP0.6.4.GC": (1.0, 1.0, 1.0, 1.0)}

sample_name_dict = {
    "HP0.6.4": "Hp(0,6,4)",
    "HP064": "Hp(0,6,4)",
    "HP0.10.4": "Hp(0,10,4)",
    "HP0104": "Hp(0,10,4)",
    "HP0.4.4": "Hp(0,4,4)",
    "HP044": "Hp(0,4,4)",
    "HP0.0.4": "Hp(0,0,4)",
    "HP004": "Hp(0,0,4)",
    "HP060": "Hp(0,6,0)",
    "HP0.6.0": "Hp(0,6,0)",
    "HP068": "Hp(0,6,8)",
    "HP0.6.8": "Hp(0,6,8)",
    "HP464": "Hp(4,6,4)",
    "HP4.6.4": "Hp(4,6,4)",
    "HP564": "Hp(5,6,4)",
    "HP5.6.4": "Hp(5,6,4)",
    "HP4.0.4": "Hp(4,0,4)",
    "HP040": "Hp(4,0,4)",
    "HP0.6.4.GC": "Hp(0,6,4) mod"
}



# Calculate "Selectivities" for bot strand HP(0,6,4)
sample_names = [f"HP0.6.4_{run}_Breaking_Dist_Full_Shifted.csv" for run in [1, 2, 3]]
breaking_results = read_breaking_results(sample_names, proc_data_dir)
sample_names_top = [f"HP0.6.4_{run}#BOT" for run in [1, 2, 3]]
selectivities = []
for break_dist in [breaking_results[key] for key in sample_names_top]:
    total_reads = sum([break_dist[key] for key in break_dist.keys() if -300<=key<=300])
    reads_in_window = sum([break_dist[key] for key in break_dist.keys() if -2<=key<=2])
    selectivity_cur = 100*reads_in_window/total_reads
    selectivities.append(selectivity_cur)
selectivity_std = np.std(selectivities)
selectivity_mean = np.mean(selectivities)
print("Yields HP(0,6,4) bot std: ", selectivity_std, " mean: ", selectivity_mean)

# Calculate "Selectivities" for dsDNA
sample_names = [f"dsDNA_{run}_Breaking_Dist_Full_Shifted.csv" for run in [1, 2]]
breaking_results = read_breaking_results(sample_names, "../ProcessedData/Double_Stranded_DNA/")
sample_names_bot = [f"dsDNA_{run}#BOT" for run in [1, 2]]
selectivities = []
for break_dist in [breaking_results[key] for key in sample_names_bot]:
    total_reads = sum([break_dist[key] for key in break_dist.keys() if -300<=key<=300])
    reads_in_window = sum([break_dist[key] for key in break_dist.keys() if -2<=key<=2])
    selectivity_cur = 100*reads_in_window/total_reads
    selectivities.append(selectivity_cur)
selectivity_std = np.std(selectivities)
selectivity_mean = np.mean(selectivities)
print("Yields dsDNA bot std: ", selectivity_std, " mean: ", selectivity_mean)

# Calculate "Selectivities" for nickedDNA
sample_names = [f"704_{run}_Breaking_Dist_Full_Shifted.csv" for run in [1, 2, 3]]
breaking_results = read_breaking_results(sample_names, "../ProcessedData/704_AT_GC_400_1500_ATMi/")
sample_names_bot = [f"704_{run}#TOP" for run in [1, 2, 3]]
selectivities = []
for break_dist in [breaking_results[key] for key in sample_names_bot]:
    total_reads = sum([break_dist[key] for key in break_dist.keys() if -300<=key<=300])
    reads_in_window = sum([break_dist[key] for key in break_dist.keys() if -2<=key<=2])
    selectivity_cur = 100*reads_in_window/total_reads
    selectivities.append(selectivity_cur)
selectivity_std = np.std(selectivities)
selectivity_mean = np.mean(selectivities)
print("Yields 704* top, std: ", selectivity_std, " mean: ", selectivity_mean)


# Overlay histogram stem
files = [f"{base_name}_Breaking_Dist_Norm_Shifted_Average.csv" for base_name in stem_samples]
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="TOP", colors=list(colors.values()),
                                bars=[bars_top[base_name] for base_name in stem_samples],
                                save_path=f"{graphics_dir}/Stem_Variation_Overlay_Histogram_Top.pdf",
                                sample_names=sample_name_dict, legend_size=3, fits=["t-student"] * len(stem_samples))
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="BOT", colors=list(colors.values()),
                                bars=[bars_bot[base_name] for base_name in stem_samples],
                                save_path=f"{graphics_dir}/Stem_Variation_Overlay_Histogram_Bot.pdf",
                                sample_names=sample_name_dict, legend_size=3, fits=["t-student"] * len(stem_samples))

# Overlay histogram loop
colors_loop_top = [(0.1216, 0.3059, 0.8471, 1.0), (0.1490, 0.6510, 0.6039, 1.0), (1.0000, 0.8235, 0.2471, 1.0)]
colors_loop_bot = [(0.9569, 0.6510, 0.7569, 1.0), (0.580, 0.404, 0.741, 1.0), (0.5451, 0.3529, 0.1686, 1.0)]
files = [f"{base_name}_Breaking_Dist_Norm_Shifted_Average.csv" for base_name in loop_samples]
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="TOP", colors=colors_loop_top,
                                bars=[bars_top[base_name] for base_name in loop_samples],
                                save_path=f"{graphics_dir}/Loop_Variation_Overlay_Histogram_Top.pdf",
                                sample_names=sample_name_dict, legend_size=3, fits=["t-student", "t-student", "t-bimodal"])
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="BOT", colors=colors_loop_bot,
                                bars=[bars_bot[base_name] for base_name in loop_samples],
                                save_path=f"{graphics_dir}/Loop_Variation_Overlay_Histogram_Bot.pdf",
                                sample_names=sample_name_dict, legend_size=3, fits=["t-student"] * len(loop_samples))

# Overlay histogram upaired
colors_unpaired_top = [(0.1490, 0.6510, 0.6039, 1.0), (0.7608, 0.0941, 0.5451, 1.0), (1.0000, 0.4353, 0.3804, 1.0)]
colors_unpaired_bot = [(0.580, 0.404, 0.741, 1.0), (0.9020, 0.8275, 0.6392, 1.0), (0.7529, 0.7529, 0.7529, 1.0)]
files = [f"{base_name}_Breaking_Dist_Norm_Shifted_Average.csv" for base_name in unpaired_samples]
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="TOP", colors=colors_unpaired_top,
                                bars=[bars_top[base_name] for base_name in unpaired_samples],
                                save_path=f"{graphics_dir}/Unpaired_Variation_Overlay_Histogram_Top.pdf",
                                sample_names=sample_name_dict, legend_size=3, fits=["t-student", "t-student", "t-student"])
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="BOT", colors=colors_unpaired_bot,
                                bars=[bars_bot[base_name] for base_name in unpaired_samples],
                                save_path=f"{graphics_dir}/Unpaired_Variation_Overlay_Histogram_Bot.pdf",
                                sample_names=sample_name_dict, legend_size=3, fits=["t-student"] * len(unpaired_samples))

# Breaking distributions
for sample_group in sample_groups:
    if sample_group == "HP0.6.8":
        fit_type_top = "t-bimodal"
    else:
        fit_type_top = "t-student"

    if sample_group == "HP0.6.4.GC":
        lengend = "t-student"
    else:
        lengend = None

    plotter.plot_three_breaking([f"{sample_group}_{i}_Breaking_Dist_Norm_Shifted.csv" for i in range(1, 4)],
                                proc_data_dir, "TOP", bars=bars_top[sample_group],
                                save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Three_Top.pdf",
                                fit=fit_type_top)
    plotter.plot_three_breaking([f"{sample_group}_{i}_Breaking_Dist_Norm_Shifted.csv" for i in range(1, 4)],
                                proc_data_dir, "BOT", bars=bars_bot[sample_group],
                                save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Three_Bot.pdf",
                                fit="t-student")

    plotter.plot_single_breaking_data_average(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                              "TOP", bars=bars_top[sample_group],
                                              save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Top.pdf",
                                              fit=fit_type_top, legend=lengend)
    plotter.plot_single_breaking_data_average(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                              "BOT", bars=bars_bot[sample_group],
                                              save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Bot.pdf",
                                              fit="t-student", legend=lengend)

# Quality plotting
for sample_group in sample_groups:
    sample_names = [f"{sample_group}_{run}" for run in [1, 2, 3]]
    plotter.plot_quality_values([f"{sample_name}_Quality.csv" for sample_name in sample_names], proc_data_dir, sample_name_dict,
                                save_path=f"{graphics_dir}/{sample_group}_Quality.pdf")

# Plot nicked DNA data
graphics_directory = "../Graphics/704_AT_GC_400_1500_ATMi"
plotter.plot_single_breaking_data_average(
    file_name_in=f"704_Breaking_Dist_Norm_Shifted_Average.csv",
    proc_directory_in="../ProcessedData/704_AT_GC_400_1500_ATMi/", strand="TOP", fit="t-student",
    save_path=f"{graphics_directory}/Reference_Breaking_Dist_Norm_Shifted_Average_Top.pdf", legend="t-student")

# Scheme plots
for sample_group in ["HP4.0.4", "HP0.6.4"]:
    plotter.plot_single_breaking_scheme(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                 "TOP", bars=bars_top[sample_group], save_path=f"{graphics_dir}/{sample_group}_Scheme_Top.pdf", color=(178/255, 14/255, 14/255))
    plotter.plot_single_breaking_scheme(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                 "BOT", bars=bars_bot[sample_group], save_path=f"{graphics_dir}/{sample_group}_Scheme_Bot.pdf", color=(0/255, 84/255, 133/255))


