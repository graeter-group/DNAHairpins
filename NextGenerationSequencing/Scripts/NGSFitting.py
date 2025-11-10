# Script for handling Fitting
# Author: Boris N. Schüpp

from scipy.optimize import curve_fit
import numpy as np
from scipy.stats import t


def student_t_dist(x, df, mu, sigma):
    return t.pdf(x, df, loc=mu, scale=sigma)


def t_widths(df, loc=0, scale=1):
    iqr = t.ppf(0.75, df, loc, scale) - t.ppf(0.25, df, loc, scale)
    return iqr


def fit_student(data_in, provide_data=False):
    x = [x_i for x_i in sorted(list(data_in.keys())) if data_in[x_i] != 0]
    y = [data_in[x_i] for x_i in x]
    popt_t = curve_fit(student_t_dist, x, y, p0=[5, 0, 10])

    df, mu_t, sigma_t = popt_t[0]

    x_fit = np.linspace(min(x), max(x), 1000)
    y_fit = [student_t_dist(i, df, mu_t, sigma_t) for i in x_fit]

    if provide_data:
        return df, mu_t, sigma_t, t_widths(df, mu_t, sigma_t), x_fit, y_fit
    else:
        return df, mu_t, sigma_t, t_widths(df, mu_t, sigma_t)


def rsquares_student(data_in, df_in, mu_in, sigma_in):
    x = [x_i for x_i in sorted(list(data_in.keys())) if data_in[x_i] != 0]
    y = [data_in[x_i] for x_i in x]
    mean = sum(y) / len(y)
    ssres = sum([(y_i - student_t_dist(x_i, df_in, mu_in, sigma_in)) ** 2 for x_i, y_i in zip(x, y)])
    sstot = sum([(y_i - mean) ** 2 for y_i in y])
    return 1 - ssres / sstot
