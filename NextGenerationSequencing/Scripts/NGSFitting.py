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

def student_t_mix2(x, a1, df1, mu1, sigma1, a2, df2, mu2, sigma2):
    x = np.asarray(x)
    return (a1 * student_t_dist(x, df1, mu1, sigma1) +
            a2 * student_t_dist(x, df2, mu2, sigma2))

def fit_student_bimodal(data_in, provide_data=False, p0=None, bounds=None):
    x = np.array([xi for xi in sorted(data_in.keys()) if data_in[xi] != 0], dtype=float)
    y = np.array([data_in[xi] for xi in x], dtype=float)

    if x.size < 6:
        raise ValueError("Need more data points to fit 2 components reliably.")

    if p0 is None:
        idx_sorted = np.argsort(y)[::-1]
        mu1 = x[idx_sorted[0]]
        mu2 = None
        for j in idx_sorted[1:]:
            if abs(x[j] - mu1) > 0.1 * (x.max() - x.min() + 1e-12):
                mu2 = x[j]
                break
        if mu2 is None:
            mu2 = x[idx_sorted[1]]

        a1 = float(y.max())
        a2 = float(np.partition(y, -2)[-2]) if y.size >= 2 else float(y.max()) * 0.5

        spread = float(np.std(x)) if np.std(x) > 0 else (x.max() - x.min()) / 10.0
        sigma1 = spread if spread > 0 else 1.0
        sigma2 = spread if spread > 0 else 1.0

        df1 = 5.0
        df2 = 5.0

        p0 = [a1, df1, mu1, sigma1, a2, df2, mu2, sigma2]

    if bounds is None:
        pad = 0.1 * (x.max() - x.min() + 1e-12)
        lower = [0.0, 1e-3, x.min() - pad, 1e-6,
                 0.0, 1e-3, x.min() - pad, 1e-6]
        upper = [np.inf, 200.0, x.max() + pad, np.inf,
                 np.inf, 200.0, x.max() + pad, np.inf]
        bounds = (lower, upper)

    popt, pcov = curve_fit(
        student_t_mix2, x, y,
        p0=p0,
        bounds=bounds,
        maxfev=20000
    )

    a1, df1, mu1, sigma1, a2, df2, mu2, sigma2 = popt

    w1 = t_widths(df1, mu1, sigma1)
    w2 = t_widths(df2, mu2, sigma2)

    if mu2 < mu1:
        a1, a2 = a2, a1
        df1, df2 = df2, df1
        mu1, mu2 = mu2, mu1
        sigma1, sigma2 = sigma2, sigma1
        w1, w2 = w2, w1

    if provide_data:
        x_fit = np.linspace(x.min(), x.max(), 1000)
        y_fit = student_t_mix2(x_fit, a1, df1, mu1, sigma1, a2, df2, mu2, sigma2)
        y1 = a1 * student_t_dist(x_fit, df1, mu1, sigma1)
        y2 = a2 * student_t_dist(x_fit, df2, mu2, sigma2)
        return (a1, df1, mu1, sigma1, w1), (a2, df2, mu2, sigma2, w2), x_fit, y_fit, y1, y2
    else:
        return (a1, df1, mu1, sigma1, w1), (a2, df2, mu2, sigma2, w2)


def rsquares_student(data_in, df_in, mu_in, sigma_in):
    x = [x_i for x_i in sorted(list(data_in.keys())) if data_in[x_i] != 0]
    y = [data_in[x_i] for x_i in x]
    mean = sum(y) / len(y)
    ssres = sum([(y_i - student_t_dist(x_i, df_in, mu_in, sigma_in)) ** 2 for x_i, y_i in zip(x, y)])
    sstot = sum([(y_i - mean) ** 2 for y_i in y])
    return 1 - ssres / sstot
