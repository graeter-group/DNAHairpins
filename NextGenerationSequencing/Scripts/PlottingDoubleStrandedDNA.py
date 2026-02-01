from NGSPlotting import *

runs = [1, 2]
plotter = PlottingNGS()

proc_directory_in = "../ProcessedData/Double_Stranded_DNA/"
graphics_directory = "../Graphics/Double_Stranded_DNA/"

sample_group = "dsDNA"
sample_names = [f"dsDNA_{run}" for run in runs]
plotter.plot_quality_values([f"{sample_name}_Quality.csv" for sample_name in sample_names], proc_directory_in,
                    save_path=f"{graphics_directory}/dsDNA_Quality.pdf")


plotter.plot_single_breaking_data(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_directory_in,
                                          "TOP",
                                          save_path=f"{graphics_directory}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Top.pdf")
plotter.plot_single_breaking_data(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_directory_in,
                                          "BOT",
                                          save_path=f"{graphics_directory}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Bot.pdf")

plotter.plot_single_breaking_data_window(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_directory_in,
                                          "TOP",
                                          save_path=f"{graphics_directory}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Top_Window.pdf")
plotter.plot_single_breaking_data_window(f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", proc_directory_in,
                                          "BOT",
                                          save_path=f"{graphics_directory}/{sample_group}_Breaking_Dist_Norm_Shifted_Average_Bot_Window.pdf")