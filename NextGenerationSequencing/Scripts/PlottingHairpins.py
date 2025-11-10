from NGSPlotting import *

sample_groups = ["HP4.0.4", "HP0.0.4", "HP0.4.4", "HP0.6.4", "HP0.10.4"]

graphics_dir = "../Graphics/Hairpins/"
proc_data_dir = "../ProcessedData/Hairpins/"

runs = [1, 2, 3]
plotter = PlottingNGS()

# Ligation positions
bars_top = {"HP0.0.4": [-28.5, 29.5], "HP0.4.4": [-32.5, 33.5], "HP0.10.4": [-38.5, 39.5],
            "HP4.0.4": [-28.5, 29.5], "HP0.6.4": [-34.5, 35.5]}

bars_bot = {"HP0.0.4": [-30.5, 31.5], "HP0.4.4": [-30.5, 31.5], "HP0.10.4": [-30.5, 31.5],
            "HP4.0.4": [-32.5, 33.5], "HP0.6.4": [-30.5, 31.5]}

stem_samples = ["HP0.10.4", "HP0.6.4", "HP0.4.4", "HP0.0.4"]

colors = {"HP0.10.4": (0.192, 0.627, 0.204, 1.0),
          "HP0.6.4": (0.580, 0.404, 0.741, 1.0),
          "HP0.4.4": (0.690196, 0.0, 0.125490, 1.0),
          "HP0.0.4": (1.000, 0.498, 0.0, 1.0)}

sample_name_dict = {
    "HP0.6.4": "Hp(0,6,4)",
    "HP064": "Hp(0,6,4)",
    "HP0.10.4": "Hp(0,10,4)",
    "HP0104": "Hp(0,10,4)",
    "HP0.4.4": "Hp(0,4,4)",
    "HP044": "Hp(0,4,4)",
    "HP0.0.4": "Hp(0,0,4)",
    "HP004": "Hp(0,0,4)",
}

# Overlay histogram

files = [f"{base_name}_Breaking_Dist_Norm_Shifted_Average.csv" for base_name in stem_samples]
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="TOP", colors=list(colors.values()),
                                bars=[bars_top[base_name] for base_name in stem_samples],
                                save_path=f"{graphics_dir}/Stem_Variation_Overlay_Histogram_Top.pdf",
                                sample_names=sample_name_dict, legend_size=3)
plotter.plot_overlay_histogramm(files, proc_data_dir, strand="BOT", colors=list(colors.values()),
                                bars=[bars_bot[base_name] for base_name in stem_samples],
                                save_path=f"{graphics_dir}/Stem_Variation_Overlay_Histogram_Bot.pdf",
                                sample_names=sample_name_dict)

# Breaking distributions
for sample_group in sample_groups:
    plotter.plot_three_breaking([f"{sample_group}_{i}_Breaking_Dist_Norm_Shifted.csv" for i in range(1, 4)],
                                proc_data_dir, "TOP", bars=bars_top[sample_group],
                                save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Three_Top.pdf",
                                fit="ON")
    plotter.plot_three_breaking([f"{sample_group}_{i}_Breaking_Dist_Norm_Shifted.csv" for i in range(1, 4)],
                                proc_data_dir, "BOT", bars=bars_bot[sample_group],
                                save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Three_Bot.pdf",
                                fit="ON")

    plotter.plot_single_breaking_data_average(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                              "TOP", bars=bars_top[sample_group],
                                              save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Top.pdf",
                                              fit="ON")
    plotter.plot_single_breaking_data_average(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                              "BOT", bars=bars_bot[sample_group],
                                              save_path=f"{graphics_dir}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Bot.pdf",
                                              fit="ON")

# Quality plotting
for sample_group in sample_groups:
    sample_names = [f"{sample_group}_{run}" for run in [1, 2, 3]]
    plotter.plot_quality_values([f"{sample_name}_Quality.csv" for sample_name in sample_names], proc_data_dir,
                                save_path=f"{graphics_dir}/{sample_group}_Quality.pdf")

# Plot nicked DNA data
graphics_directory = "../Graphics/704_AT_GC_400_1500_ATMi"
plotter.plot_single_breaking_data_average(
    file_name_in=f"704_Breaking_Dist_Norm_Shifted_Average.csv",
    proc_directory_in="../ProcessedData/704_AT_GC_400_1500_ATMi/", strand="TOP", fit="ON",
    save_path=f"{graphics_directory}/Reference_Breaking_Dist_Norm_Shifted_Average_Top.pdf", legend="ON")

# Scheme plots
for sample_group in ["HP4.0.4", "HP0.6.4"]:
    plotter.plot_single_breaking_scheme(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                 "TOP", bars=bars_top[sample_group], save_path=f"{graphics_dir}/{sample_group}_Scheme_Top.pdf", color=(178/255, 14/255, 14/255))
    plotter.plot_single_breaking_scheme(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_data_dir,
                                 "BOT", bars=bars_bot[sample_group], save_path=f"{graphics_dir}/{sample_group}_Scheme_Bot.pdf", color=(0/255, 84/255, 133/255))