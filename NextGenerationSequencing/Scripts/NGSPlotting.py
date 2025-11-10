# Utilities for plotting results from the NGS analysis section
# Author: Boris N. Schüpp

import matplotlib
from PlottingPreferences import *
from NGSDataUtilities import *
from NGSFitting import *

matplotlib.use("Agg")


class PlottingNGS:
    def __init__(self):
        self.colors = {"BlackOpaque": (0, 0, 0, 0.8),
                       "BlackFaded": (0, 0, 0, 0.5),
                       "DarkGrey": (0.43921569, 0.43921569, 0.44313725, 1),
                       "LightGrey": (0.85490196, 0.85490196, 0.85882352, 1)}

    @staticmethod
    def get_complementary_strand(strand_in):
        complement_35 = ""
        for base in strand_in:
            if base == "A":
                complement_35 += "T"
            if base == "T":
                complement_35 += "A"
            if base == "G":
                complement_35 += "C"
            if base == "C":
                complement_35 += "G"
        complement_53 = complement_35[::-1]
        return complement_53

    def plot_three_breaking(self, file_names_in, directory_in, strand, fit=None, bars=None, save_path=None):

        breaking_data = read_breaking_results(file_names_in, directory_in)
        data = [breaking_data[i] for i in breaking_data.keys() if i.find(strand) != -1]

        if len(data) == 0:
            data = [i for i in breaking_data.values()]

        y_limit = max([max(list(data[i].values())) for i in range(0, 3)]) + 0.001
        apply_plot_config(wide_plot_config)

        fig, axs = plt.subplots(1, 3)
        axis = None

        for i in range(0, 3):
            cur_data = data[i]
            cur_data_x = sorted(list(cur_data.keys()))
            cur_data_y = [cur_data[x_i] for x_i in cur_data_x]
            axis = axs[i]

            x_all = np.array(sorted({float(v) for v in cur_data_x}), dtype=float)
            if len(x_all) >= 2:
                step = float(np.min(np.diff(x_all)))
            else:
                step = 1.0
            bar_width = step * 0.9

            def to_left_edges(xc):
                xc = np.asarray(xc, dtype=float)
                return xc - bar_width / 2.0

            if bars:
                old_bar_values = []
                bar_1 = bars[0]
                bar_2 = bars[1]
                neighbour_hood_1 = [cur_data_y[cur_data_x.index(x)] for x in cur_data_x if
                                    bar_1 - 3 <= x <= bar_1 + 3 and x != bar_1]
                bar_1_replacement = sum(neighbour_hood_1) / len(neighbour_hood_1)
                old_bar_values.append(cur_data_y[cur_data_x.index(bar_1)])
                neighbour_hood_2 = [cur_data_y[cur_data_x.index(x)] for x in cur_data_x if
                                    bar_2 - 3 <= x <= bar_2 + 3 and x != bar_2]
                bar_2_replacement = sum(neighbour_hood_2) / len(neighbour_hood_2)
                old_bar_values.append(cur_data_y[cur_data_x.index(bar_2)])
                cur_data_y[cur_data_x.index(bar_1)] = bar_1_replacement
                cur_data_y[cur_data_x.index(bar_2)] = bar_2_replacement
                total_sum = sum(cur_data_y)
                cur_data_y = [y_i / total_sum for y_i in cur_data_y]

                axis.bar(to_left_edges([bar_1, bar_2]), old_bar_values,
                         color=self.colors["LightGrey"],
                         edgecolor=self.colors["LightGrey"], width=bar_width, linewidth=0.2)

            plt.sca(axis)
            axis.bar(to_left_edges(cur_data_x), cur_data_y, linewidth=0.2, width=bar_width, align='edge')

            axis.set_ylabel("Fraction of reads")
            axis.set_xlabel("Relative fragment break index")

            axis.set_xlim(min(cur_data_x), max(cur_data_x))

            x_ticks = [i for i in range(int(min(cur_data_x)), int(max(cur_data_x))) if i % 20 == 0]
            x_tick_labels = []
            for x_tick in x_ticks:
                if x_tick > 0:
                    x_tick_labels.append(f"+{str(int(x_tick))}")
                else:
                    x_tick_labels.append(f"{str(int(x_tick))}")

            axis.set_xticks(x_ticks, x_tick_labels)

            if fit:
                df, mu, sigma, iqr, x_fit, y_fit = fit_student({x_i: y_i for x_i, y_i in zip(cur_data_x, cur_data_y)},
                                                               True)
                plt.plot(x_fit, y_fit, color="red", linewidth=1.0)
                y_limit = max(max(y_fit) + 0.001, y_limit)

        plt.tight_layout()

        axis.set_ylim(0, y_limit + 0.001)

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_single_breaking_data_average(self, file_name_in, proc_directory_in, strand, fit=None, bars=None,
                                          color=(0, 0, 0, 1), save_path=None, legend=None):

        # Read data
        breaking_data = read_breaking_results([file_name_in], proc_directory_in)
        if strand == "TOP":
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1 and i.find("AVERAGE") != -1][
                0]
            data_er = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1 and i.find("STD") != -1][0]
        else:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1 and i.find("AVERAGE") != -1][
                0]
            data_er = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1 and i.find("STD") != -1][0]

        if len(data) == 0:
            print("Data is not correct, please try again")

        data = {i: data[i] for i in data.keys() if data[i] != 0}
        x_data = list(data.keys())
        y_data = [data[i] for i in x_data]
        data_er = [data_er[i] for i in x_data]

        apply_plot_config(default_plot_config)

        fig, ax = plt.subplots()

        x_all = np.array(sorted({float(v) for v in x_data}), dtype=float)
        if len(x_all) >= 2:
            step = float(np.min(np.diff(x_all)))
        else:
            step = 1.0
        bar_width = step * 0.75

        def to_left_edges(xc):
            xc = np.asarray(xc, dtype=float)
            return xc - bar_width / 2.0

        # Grey out artifacts
        old_data_y = [0]
        if bars:
            old_bar_values = []
            old_bar_values_er = []
            bar_1 = bars[0]
            bar_2 = bars[1]
            neighbour_hood_1 = [y_data[x_data.index(x)] for x in x_data if bar_1 - 3 <= x <= bar_1 + 3 and x != bar_1]
            bar_1_replacement = (sum(neighbour_hood_1) / len(neighbour_hood_1)) if neighbour_hood_1 else y_data[
                x_data.index(bar_1)]
            old_bar_values.append(y_data[x_data.index(bar_1)])
            neighbour_hood_2 = [y_data[x_data.index(x)] for x in x_data if bar_2 - 3 <= x <= bar_2 + 3 and x != bar_2]
            bar_2_replacement = (sum(neighbour_hood_2) / len(neighbour_hood_2)) if neighbour_hood_2 else y_data[
                x_data.index(bar_2)]
            old_bar_values.append(y_data[x_data.index(bar_2)])
            y_data[x_data.index(bar_1)] = bar_1_replacement
            y_data[x_data.index(bar_2)] = bar_2_replacement
            old_bar_values_er.append(data_er[x_data.index(bar_1)])
            old_bar_values_er.append(data_er[x_data.index(bar_2)])
            data_er[x_data.index(bar_1)] = 0
            data_er[x_data.index(bar_2)] = 0
            total_sum = sum(y_data)
            if total_sum > 0:
                y_data = [y_i / total_sum for y_i in y_data]
                data_er = [y_er_i / total_sum for y_er_i in data_er]

            ax.bar(
                to_left_edges([bar_1, bar_2]), old_bar_values,
                yerr=old_bar_values_er,
                color=self.colors["LightGrey"], edgecolor=self.colors["LightGrey"], linewidth=0.1,
                width=bar_width, align='edge', antialiased=False,
                capsize=1,
                error_kw={
                    'ecolor': self.colors["LightGrey"],
                    'elinewidth': 0.2,  # thin stem (use 0 to hide completely)
                    'capthick': 0.2
                },
            )
            old_data_y = [y + y_er + 0.001 for y, y_er in zip(old_bar_values, old_bar_values_er)]

        y_plus_err = [y + x for y, x in zip(y_data, data_er)]
        y_limit = max(y_plus_err) + 0.001
        y_limit = max(y_limit, max(old_data_y))

        # main bars with the same width + caps-only error bars
        ax.bar(
            to_left_edges(x_data), y_data,
            yerr=data_er, color=color,
            edgecolor='none', linewidth=1,  # match your overlay look
            width=bar_width, align='edge', antialiased=False,
            capsize=1,
            error_kw={
                'ecolor': self.colors["DarkGrey"],
                'elinewidth': 0.2,  # thin stem (use 0 to hide)
                'capthick': 0.2
            },
        )

        # Plot fit
        mu = None
        iqr = None
        if fit:
            df, mu, sigma, iqr, x_fit, y_fit = fit_student(dict(zip(x_data, y_data)), True)
            label = f"Student's $t$ fit\nIQR = {iqr:.2f}"
            ax.plot(x_fit, y_fit, color="red", linewidth=0.6, label=label)
            y_limit = max(max(y_fit) + 0.001, y_limit)

        # Format figure
        ax.set_xlabel("Relative fragment break index")
        ax.set_ylabel("Fraction of reads")
        ax.set_xlim(min(x_data), max(x_data))
        ax.set_ylim(0, y_limit)
        x_ticks = [i for i in range(int(min(x_data)), int(max(x_data))) if i % 20 == 0]
        x_tick_labels = [str(round(int(i))) for i in x_ticks]
        x_tick_labels_fixed = [("+" + t_i) if int(t_i) > 0 else t_i for t_i in x_tick_labels]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels_fixed)

        if legend:
            ax.legend(loc="upper left", handlelength=1.0, handletextpad=0.6)  # top-right corner

        plt.tight_layout()

        # Save image
        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

        if fit:
            return mu, iqr
        else:
            return None, None

    @staticmethod
    def plot_fit_paramameter_scatter_single(file_names_in, directory_in, parameter_name, sample_names,
                                            colors, save_path=None):

        fit_parameters = read_fit_parameters(file_names_in, directory_in)
        num_of_samples = len(list(fit_parameters.keys()))
        positions = range(1, num_of_samples + 1)

        apply_plot_config(default_plot_config)
        plt.figure(figsize=(single_column, single_column))

        for sample, position in zip(fit_parameters.keys(), positions):
            sample_id = sample.strip("*")
            current_color = colors[sample_id]
            for sample_num in fit_parameters[sample].keys():
                if fit_parameters[sample][sample_num] == -100 or fit_parameters[sample][sample_num] == 0.0:
                    continue
                if sample_num != "Average":
                    plt.scatter(position, fit_parameters[sample][sample_num], alpha=0.5, color=current_color,
                                marker="o", s=25, edgecolor='black', linewidth=0.4)
                else:
                    plt.scatter(position, fit_parameters[sample][sample_num], alpha=1, color=current_color,
                                marker="*", s=30, edgecolor='black', linewidth=0.4)
        if len(list(fit_parameters.keys())) <= 5:
            plt.xticks(positions, [sample_names[i] for i in fit_parameters.keys()])
        else:
            plt.xticks(positions, [sample_names[i] for i in fit_parameters.keys()])

        if parameter_name == "Mu":
            plt.ylabel("$\\mu$ difference from center position")
        elif parameter_name == "IQR":
            plt.ylabel("IQR")

        total_min = min([min(fit_parameters[i].values()) for i in fit_parameters.keys()])
        total_max = max([max(fit_parameters[i].values()) for i in fit_parameters.keys()])
        y_ticks = [i for i in range(int(round(total_min, 0)), int(round(total_max, 0)) + 1)]
        plt.yticks(y_ticks, [str(i) for i in y_ticks])

        plt.xlabel("Sample")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    def plot_overlay_histogramm(self, file_names_in, directory_in, colors, strand, sample_names, save_path=None,
                                bars=None, legend_size=None):

        data_breaking_all = read_breaking_results(file_names_in, directory_in)
        plotting_dict = {}

        for key in data_breaking_all.keys():
            sample = key.split("#")[0]
            position = key.split("#")[1]
            data_type = key.split('#')[2]
            if position != strand:
                continue
            if sample not in plotting_dict:
                plotting_dict[sample] = {}
            plotting_dict[sample][data_type] = data_breaking_all[key]

        apply_plot_config(default_plot_config)
        fig, ax = plt.subplots()

        ax.axhline(y=0.03, color='black', linewidth=1.0, zorder=1)
        ax.axhline(y=0.06, color='black', linewidth=1.0, zorder=1)

        all_x = []
        for s in plotting_dict:
            all_x.extend(plotting_dict[s]["AVERAGE"].keys())
        all_x = np.array(sorted({float(v) for v in all_x}), dtype=float)
        if len(all_x) >= 2:
            step = float(np.min(np.diff(all_x)))
        else:
            step = 1.0
        bar_width = step

        def to_left_edges(xc):
            xc = np.asarray(xc, dtype=float)
            return xc - bar_width / 2.0

        shift = len(file_names_in) * 0.03 - 0.03
        x_min, x_max = 1e9, -1e9

        for sample, color in zip(plotting_dict.keys(), colors):
            x_centers = list(plotting_dict[sample]["AVERAGE"].keys())
            y = list(plotting_dict[sample]["AVERAGE"].values())
            y_er = list(plotting_dict[sample]["STD"].values())
            sample_name = sample_names[sample]

            # Grey out artifacts
            if bars:
                cur_bars = bars[colors.index(color)]
                old_bar_values = []
                old_bar_values_er = []
                bar_1 = cur_bars[0]
                bar_2 = cur_bars[1]
                neighbour_hood_1 = [y[x_centers.index(x_i)] for x_i in x_centers
                                    if bar_1 - 3 <= x_i <= bar_1 + 3 and x_i != bar_1]
                bar_1_replacement = sum(neighbour_hood_1) / len(neighbour_hood_1) if neighbour_hood_1 else y[
                    x_centers.index(bar_1)]
                old_bar_values.append(y[x_centers.index(bar_1)])
                neighbour_hood_2 = [y[x_centers.index(x_i)] for x_i in x_centers
                                    if bar_2 - 3 <= x_i <= bar_2 + 3 and x_i != bar_2]
                bar_2_replacement = sum(neighbour_hood_2) / len(neighbour_hood_2) if neighbour_hood_2 else y[
                    x_centers.index(bar_2)]
                old_bar_values.append(y[x_centers.index(bar_2)])

                y[x_centers.index(bar_1)] = bar_1_replacement
                y[x_centers.index(bar_2)] = bar_2_replacement

                old_bar_values_er.append(y_er[x_centers.index(bar_1)])
                old_bar_values_er.append(y_er[x_centers.index(bar_2)])
                y_er[x_centers.index(bar_1)] = 0
                y_er[x_centers.index(bar_2)] = 0

                total_sum = sum(y)
                if total_sum > 0:
                    y = [yy / total_sum for yy in y]
                    y_er = [ee / total_sum for ee in y_er]

                ax.bar(
                    to_left_edges([bar_1, bar_2]), old_bar_values,
                    bottom=shift, yerr=old_bar_values_er,
                    color=self.colors["LightGrey"], edgecolor='black', linewidth=0.1,
                    width=bar_width * 0.75, align='edge', antialiased=False,
                    capsize=1,
                    error_kw={
                        'ecolor': 'black',
                        'elinewidth': 0.2,
                        'capthick': 0.2
                    },
                )

            ax.bar(
                to_left_edges(x_centers), y,
                bottom=shift, label=sample_name, zorder=2, yerr=y_er,
                color=color, width=bar_width * 0.75, align='edge',
                edgecolor='none', linewidth=1, antialiased=False,
                capsize=1,
                error_kw={
                    'ecolor': self.colors["BlackFaded"],
                    'elinewidth': 0.2,
                    'capthick': 0.2
                },
            )

            x_min = min(x_min, float(min(x_centers)))
            x_max = max(x_max, float(max(x_centers)))

            df, mu, sigma, iqr, x_fit, y_fit = fit_student({x_i: y_i for x_i, y_i in zip(x_centers, y)}, True)
            y_fit = [v + shift for v in y_fit]
            ax.plot(x_fit, y_fit, color="black", linewidth=0.6)
            shift -= 0.03
            ax.set_ylim(0)

        ax.set_xlim(left=-49, right=49)
        ax.yaxis.set_ticks([])
        for s in ('top', 'right', 'bottom', 'left'):
            ax.spines[s].set_visible(True)
        ax.set_ylabel("Fraction of reads")

        x_ticks = [i for i in range(int(x_min), int(x_max)) if i % 20 == 0]
        x_tick_labels = [f"+{int(t_i)}" if t_i > 0 else f"{int(t_i)}" for t_i in x_ticks]
        plt.xticks(x_ticks, x_tick_labels)

        if not legend_size:
            plt.legend()
        else:
            plt.legend(fontsize=legend_size)

        ax.set_xlabel("Relative break index")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    @staticmethod
    def plot_quality_values(file_names_in, directory_in, save_path=None):
        quality_results = read_quality_values(file_names_in, directory_in)
        apply_plot_config(default_plot_config)
        plt.figure()
        for sample_name in quality_results.keys():
            cur_sample = quality_results[sample_name]
            current_x = sorted(list(cur_sample.keys()))
            current_y = [cur_sample[i] for i in current_x]
            plt.plot(current_x, current_y, label=sample_name)

        plt.xlabel("Distance from 5'-end [bases]")
        plt.ylabel("Phred quality value")
        plt.legend()

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    @staticmethod
    def plot_p_value_matrix(file_name_in, directory_in, save_path=None):
        matrix, names_fixed = read_p_matrix(file_name_in, directory_in)
        apply_plot_config(pmatrix_plot_config)
        fig, ax = plt.subplots()

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                addition = ""
                if matrix[i][j] < 0.05:
                    addition = "*"
                if matrix[i][j] < 0.01:
                    addition = "**"
                if matrix[i][j] < 0.001:
                    addition = "***"
                if matrix[i][j] < 0.0001:
                    addition = "****"
                if addition != "":
                    if matrix[i][j] < 0.001:
                        ax.text(j, i, f"{matrix[i][j]:.1e}{addition}", va='center', ha='center', color='black',
                                fontweight="bold", fontsize=2)
                    else:
                        ax.text(j, i, f"{matrix[i][j]:.3f}{addition}", va='center', ha='center', color='black',
                                fontweight="bold", fontsize=2)

                else:
                    ax.text(j, i, f"{matrix[i][j]:.3f}", va='center', ha='center', color='black', fontsize=2)
        ax.xaxis.set_ticks_position('top')
        ax.grid(which='minor', color='black', linewidth=0.3)
        ax.set_xticks([a - 0.5 for a in range(len(names_fixed))], minor=True)
        ax.set_yticks([a - 0.5 for a in range(len(names_fixed))], minor=True)
        ax.set_xticks(range(len(names_fixed)))
        ax.set_yticks(range(len(names_fixed)))
        ax.set_xlim(-0.5, len(names_fixed) - 0.5)
        ax.set_ylim(len(names_fixed) - 0.5, -0.5)
        ax.set_xticklabels(names_fixed, rotation=90)
        ax.set_yticklabels(names_fixed)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()

    @staticmethod
    def plot_single_breaking_scheme(file_name_in, proc_directory_in, strand, bars=None,
                                    color=(0, 0, 0, 1), save_path=None):
        breaking_data = read_breaking_results([file_name_in], proc_directory_in)
        if strand == "TOP":
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("TOP") != -1 and i.find("AVERAGE") != -1][
                0]
        else:
            data = [breaking_data[i] for i in breaking_data.keys() if i.find("BOT") != -1 and i.find("AVERAGE") != -1][
                0]

        data = {i: data[i] for i in data.keys() if data[i] != 0}
        x_data = list(data.keys())
        y_data = [data[i] for i in x_data]

        apply_plot_config(scheme_plot_config)

        fig, ax = plt.subplots(figsize=(single_column_scheme, single_column_scheme))

        x_all = np.array(sorted({float(v) for v in x_data}), dtype=float)
        if len(x_all) >= 2:
            step = float(np.min(np.diff(x_all)))
        else:
            step = 1.0
        bar_width = step * 0.75

        def to_left_edges(xc):
            xc = np.asarray(xc, dtype=float)
            return xc - bar_width * 5 / 2.0

        if bars:
            bar_1 = bars[0]
            bar_2 = bars[1]
            neighbour_hood_1 = [y_data[x_data.index(x)] for x in x_data if bar_1 - 3 <= x <= bar_1 + 3 and x != bar_1]
            bar_1_replacement = (sum(neighbour_hood_1) / len(neighbour_hood_1)) if neighbour_hood_1 else y_data[
                x_data.index(bar_1)]
            neighbour_hood_2 = [y_data[x_data.index(x)] for x in x_data if bar_2 - 3 <= x <= bar_2 + 3 and x != bar_2]
            bar_2_replacement = (sum(neighbour_hood_2) / len(neighbour_hood_2)) if neighbour_hood_2 else y_data[
                x_data.index(bar_2)]
            y_data[x_data.index(bar_1)] = bar_1_replacement
            y_data[x_data.index(bar_2)] = bar_2_replacement
            total_sum = sum(y_data)
            if total_sum > 0:
                y_data = [y_i / total_sum for y_i in y_data]

        y_compound = []
        x_compound = []
        for i in range(int(len(y_data) / 5)):
            y_compound.append(np.sum(y_data[i * 5:i * 5 + 5]) / 5)
            x_compound.append(x_data[i * 5 + 2])

        ax.bar(
            to_left_edges(x_compound), y_compound, color=color,
            edgecolor='none', linewidth=1,
            width=bar_width * 5, align='edge', antialiased=False
        )

        y_limit = 0.07

        # Format figure
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xlim(-50, 50)
        ax.set_ylim(0, y_limit)
        ax.tick_params(bottom=False, left=False, top=False, right=False,
                       labelbottom=False, labelleft=False)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        else:
            plt.show()
            plt.clf()
        plt.close()


if __name__ == "__main__":
    pass
