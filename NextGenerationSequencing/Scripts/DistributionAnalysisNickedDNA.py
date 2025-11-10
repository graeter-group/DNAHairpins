from NGSDataUtilities import *
from NGSFitting import *
import numpy as np

sample_groups = ["AT", "GC", "704", "ATMi0", "ATMi1", "ATMi2", "ATMi3", "ATMi4", "400", "1500"]
runs = [1, 2, 3]
processed_data_dir = "../ProcessedData/704_AT_GC_400_1500_ATMi/"

mus = {}
iqrs = {}
for sample_group in sample_groups:
    for run in runs:
        sample_name = f"{sample_group}_{run}"
        cur_breaking_data_norm_shifted = \
            read_breaking_results([f"{sample_name}_Breaking_Dist_Norm_Shifted.csv"], processed_data_dir)[
                f"{sample_name}#TOP"]
        df, mu, sigma, IQR = fit_student(cur_breaking_data_norm_shifted)
        mus[sample_name] = float(mu)
        iqrs[sample_name] = float(IQR)

    mus[f"{sample_group}_Average"] = round(
        float(np.mean([mus[key] for key in mus.keys() if key.startswith(sample_group)])), 2)
    mus[f"{sample_group}_Std"] = round(
        float(np.std([mus[key] for key in mus.keys() if key.startswith(sample_group) and key.find("Average") == -1])),
        2)
    iqrs[f"{sample_group}_Average"] = round(
        float(np.mean([iqrs[key] for key in iqrs.keys() if key.startswith(sample_group)])), 2)
    iqrs[f"{sample_group}_Std"] = round(
        float(np.std([iqrs[key] for key in iqrs.keys() if key.startswith(sample_group) and key.find("Average") == -1])),
        2)
    print(f"{sample_group}: Mu: {mus[f"{sample_group}_Average"]} +/- {mus[f"{sample_group}_Std"]}, "
          f"IQR: {iqrs[f"{sample_group}_Average"]} +/- {iqrs[f"{sample_group}_Std"]}")
