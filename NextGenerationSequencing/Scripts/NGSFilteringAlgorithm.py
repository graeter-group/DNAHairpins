from NGSConfigurationHandeler import *
from NGSLogger import *
from NGSDataUtilities import *
from collections import defaultdict


def read_input_file(input_file_name_in):
    input_file = None
    try:
        input_file = open(input_file_name_in, "r")
    except FileNotFoundError:
        print(f"ERROR: File {input_file_name_in} not found! Please check the file and the configuration!. Aborting.")
        exit(3)

    # Fastq file format: Sequence id, sequence, not relevant field, quality data
    raw_data_out = {}
    running_count = 0
    current_read = []
    for line in input_file:
        line = line.strip("\n")
        if running_count != 2:
            current_read.append(line)
        if running_count == 3:
            index = current_read[0].split()[0]
            if index in raw_data_out.keys():
                print(f"WARNING: Multiple entries for same index {index}!")
            raw_data_out[index] = [current_read[1], current_read[2]]
            current_read = []
        running_count = (running_count + 1) % 4

    input_file.close()
    return raw_data_out


def translate_quality_values(quality_string_in):
    quality_dict = {"!": 0, "\"": 1, "#": 2, "$": 3, "%": 4, "&": 5, "'": 6, "(": 7, ")": 8, "*": 9, "+": 10, ",": 11,
                    "-": 12, ".": 13, "/": 14, "0": 15, "1": 16, "2": 17, "3": 18, "4": 19, "5": 20, "6": 21, "7": 22,
                    "8": 23, "9": 24, ":": 25, ";": 26, "<": 27, "=": 28, ">": 29, "?": 30, "@": 31, "A": 32, "B": 33,
                    "C": 34, "D": 35, "E": 36, "F": 37, "G": 38, "H": 39, "I": 40}

    quality_numbers_list_out = []
    for quality_char in quality_string_in:
        quality_numbers_list_out.append(quality_dict[quality_char])

    return quality_numbers_list_out


def calculate_quality_average(quality_values_in):
    longest_read = max([len(i) for i in quality_values_in])
    quality_dict_out = {}
    for i in range(0, longest_read):
        current_count = 0
        current_quality_sum = 0
        for quality_value in quality_values_in:
            if len(quality_value) > i:
                current_quality_sum += quality_value[i]
                current_count += 1
        quality_dict_out[i] = current_quality_sum / current_count

    return quality_dict_out


def match_to_reference_single_read_single_sequence(read_in, reference_strand):
    # Check if the initial bases are found on the reference
    initial_bases = read_in[0]
    pos_strand = reference_strand.find(initial_bases)
    if pos_strand == -1:
        return -1
    else:
        result_pos = pos_strand + 1
        return result_pos


def match_reads_single(reads_in, sequences_in, sample_name_in, output_dir_in, centroids_in, gaussian_fit_window_in):
    top_sequence = sequences_in[0]
    bot_sequence = sequences_in[1]
    top_indices = defaultdict(int)
    bot_indices = defaultdict(int)
    reads_matched = defaultdict(list)
    count_invalid = 0

    for data_set in reads_in:
        for read_idx in data_set.keys():
            read_seq = data_set[read_idx]
            res_1 = match_to_reference_single_read_single_sequence(read_seq, top_sequence)
            res_2 = match_to_reference_single_read_single_sequence(read_seq, bot_sequence)

            if res_1 != -1 and res_2 != -1:
                print("WARNING: Found a single read on both strands!")

            if res_1 == -1 and res_2 == -1:
                count_invalid += 1

            if res_1 != -1:
                top_indices[res_1] += 1
                reads_matched[read_idx].append((0, res_1))

            if res_2 != -1:
                bot_indices[res_2] += 1
                reads_matched[read_idx].append((1, res_2))

    filtering_accuarcy = 100 * count_invalid / sum([len(data) for data in reads_in])
    ngs_logger.write_log(f"{sample_name_in}: Invalid reads "
                         f"{round(filtering_accuarcy, 2)} %")
    top_indices = {key: top_indices[key] for key in sorted(list(top_indices.keys()))}
    bot_indices = {key: bot_indices[key] for key in sorted(list(bot_indices.keys()))}

    # Full (unmodified)
    breaking_dists_full = {f"{sample_name}#TOP": top_indices,
                           f"{sample_name}#BOT": bot_indices}
    print_breaking(breaking_dists_full, False, f"{sample_name_in}_Breaking_Dist_Full.csv", output_dir_in)

    # Full (shifted)
    centroid_top = centroids_in[0]
    centroid_bot = centroids_in[1]
    top_indices_shifted = {key - centroid_top: top_indices[key] for key in sorted(list(top_indices.keys()))}
    bot_indices_shifted = {key - centroid_bot: bot_indices[key] for key in sorted(list(bot_indices.keys()))}
    breaking_dists_full_shifted = {f"{sample_name}#TOP": top_indices_shifted,
                                   f"{sample_name}#BOT": bot_indices_shifted}
    print_breaking(breaking_dists_full_shifted, True, f"{sample_name_in}_Breaking_Dist_Full_Shifted.csv", output_dir_in)

    # Norm (within the window)
    top_indices_restricted = {key: top_indices[key] for key in sorted(list(top_indices.keys())) if
                              centroid_top - gaussian_fit_window_in <= key <= centroid_top + gaussian_fit_window_in}
    bot_indices_restricted = {key: bot_indices[key] for key in sorted(list(bot_indices.keys())) if
                              centroid_bot - gaussian_fit_window_in <= key <= centroid_bot + gaussian_fit_window_in}
    top_indices_restricted_normalized = {key: top_indices[key] / sum(list(top_indices_restricted.values())) for key in
                                         top_indices_restricted.keys()}
    bot_indices_restricted_normalized = {key: bot_indices[key] / sum(list(bot_indices_restricted.values())) for key in
                                         bot_indices_restricted.keys()}
    breaking_dists_norm = {f"{sample_name}#TOP": top_indices_restricted_normalized,
                           f"{sample_name}#BOT": bot_indices_restricted_normalized}
    print_breaking(breaking_dists_norm, False, f"{sample_name_in}_Breaking_Dist_Norm.csv", output_dir_in)

    # Norm shifted (within the window)
    top_indices_shifted_normalized = {key - centroid_top: top_indices_restricted_normalized[key] for key in
                                      sorted(list(top_indices_restricted_normalized.keys()))}
    bot_indices_shifted_normalized = {key - centroid_bot: bot_indices_restricted_normalized[key] for key in
                                      sorted(list(bot_indices_restricted_normalized.keys()))}
    breaking_dists_norm_shifted = {f"{sample_name}#TOP": top_indices_shifted_normalized,
                                   f"{sample_name}#BOT": bot_indices_shifted_normalized}
    print_breaking(breaking_dists_norm_shifted, False, f"{sample_name_in}_Breaking_Dist_Norm_Shifted.csv",
                   output_dir_in)

    return breaking_dists_norm_shifted, reads_matched, filtering_accuarcy


def average_break_counts(break_counts_group_in):
    x_values_all = [list(i.keys()) for i in break_counts_group_in]
    combined_x_values = []
    for x_values in x_values_all:
        combined_x_values += x_values
    combined_x_values = list(set(combined_x_values))
    averaged_break_counts = {"Average": {}, "Standard Deviation": {}}
    for x_value in combined_x_values:
        values = []
        for break_counts in break_counts_group_in:
            if x_value in break_counts.keys():
                values.append(break_counts[x_value])
        averaged_break_counts["Average"][x_value] = statistics.mean(values)
        if len(values) > 1:
            averaged_break_counts["Standard Deviation"][x_value] = statistics.stdev(values)
        else:
            averaged_break_counts["Standard Deviation"][x_value] = 0

    return averaged_break_counts


# Global properties
config_file_path = "../RunConfigurations/Config_704_AT_GC_400_1500_ATMi.txt"
normalization_window = 50

# Read in configuration
configuration_good, configuration_data = read_config(config_file_path)

if not configuration_good:
    print("Configuration faulty. Aborting.")
    exit(2)

(sample_file_names, sample_names, sample_sequences, sample_sequences_second, sample_nicks,
 data_dir, graphics_dir, data_export_dir) = configuration_data

# Start logging
ngs_logger = NGSLogger(data_export_dir)
ngs_logger.write_log(
    f"Starting analysis always using the entire read and normalization_window = {normalization_window}.")

# Main analysis

sample_groups = {sample_name: sample_name.split("_")[0] for sample_name in sample_names}
paired_read_sample_paths = [[f"{data_dir}{a[0]}", f"{data_dir}{a[1]}"] for a in sample_file_names]

for (paired_read_sample_path, sample_name, sample_sequence, secondary_sequence,
     sample_nick) in zip(paired_read_sample_paths, sample_names, sample_sequences,
                         sample_sequences_second, sample_nicks):
    data_read_1 = read_input_file(paired_read_sample_path[0])
    data_read_2 = read_input_file(paired_read_sample_path[1])
    centroids = [sample_nick[1], sample_nick[3]]
    ngs_logger.write_log(f"{sample_name}: Number of raw single reads are {len(data_read_1)} (read 1) "
                         f"and {len(data_read_2)} (read 2).")

    # Quality analysis
    quality_values = [data_read_1[i][1] for i in data_read_1.keys()]
    quality_values += [data_read_2[i][1] for i in data_read_2.keys()]
    quality_values_translated = [translate_quality_values(i) for i in quality_values]
    quality_values_averages = calculate_quality_average(quality_values_translated)
    quality_values_averages_dict = {sample_name: quality_values_averages}
    print_quality_results(quality_values_averages_dict, f"{sample_name}_Quality.csv", data_export_dir)

    # Single read analysis
    both_sequences = [sample_sequence, secondary_sequence]
    cur_breaking_data_norm_shifted, cur_reads_matched, cur_filtering_acc = match_reads_single(
        [data_read_1, data_read_2],
        both_sequences, sample_name, data_export_dir,
        centroids, 50)

# Grouped analysis
group_sample_names_norm_shifted = []

for sample_group in set(list(sample_groups.values())):
    group_sample_names_norm_shifted = [f"{sample}_Breaking_Dist_Norm_Shifted.csv" for sample in sample_groups.keys() if
                                       sample_groups[sample] == sample_group]
    data_breaking_group_norm_shifted_dict = read_breaking_results(group_sample_names_norm_shifted, data_export_dir)
    data_breaking_group_norm_shifted_averaged = {}
    for strand in ["TOP", "BOT"]:
        data_breaking_group_norm_shifted = [data_breaking_group_norm_shifted_dict[i] for i in
                                            data_breaking_group_norm_shifted_dict.keys() if i.find(f"{strand}") != -1]
        data_breaking_group_norm_shifted = average_break_counts(data_breaking_group_norm_shifted)
        data_breaking_group_norm_shifted_averaged[f"{sample_group}#{strand}#AVERAGE"] = data_breaking_group_norm_shifted["Average"]
        data_breaking_group_norm_shifted_averaged[f"{sample_group}#{strand}#STD"] = data_breaking_group_norm_shifted[
            "Standard Deviation"]

    print_breaking(data_breaking_group_norm_shifted_averaged, True,
                   f"{sample_group}_Breaking_Dist_Norm_Shifted_Average.csv", data_export_dir)

ngs_logger.total_time()
