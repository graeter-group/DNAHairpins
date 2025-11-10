import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image

single_column_mm = 50
single_column_scheme_mm = 70
double_column_mm = 180
single_column_large_mm = 88

def mm_to_inch(x_mm):
    return x_mm * 0.0393700787


single_column, single_column_scheme, double_column, single_column_large = mm_to_inch(single_column_mm), mm_to_inch(
    single_column_scheme_mm), mm_to_inch(double_column_mm), mm_to_inch(single_column_large_mm)

sns.set_context("paper")
sns.set_theme(style="ticks")

default_plot_config = {
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'axes.titlesize': 5,
    'axes.labelsize': 5,
    'legend.fontsize': 5,
    'figure.titlesize': 5,

    'lines.linewidth': 1.2,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2,
    'ytick.major.size': 2,

    'xtick.top': True,
    'ytick.right': True,
    'xtick.bottom': True,
    'ytick.left': True,

    'axes.spines.top': True,
    'axes.spines.right': True,
    'axes.spines.bottom': True,
    'axes.spines.left': True,

    'xtick.direction': 'inout',
    'ytick.direction': 'inout',

    'lines.markersize': 2,
    'font.sans-serif': ['Roboto'],
    'axes.labelcolor': 'black',

    'axes.edgecolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',

    'legend.frameon': False,
    'legend.framealpha': 0,

    'savefig.dpi': 300,
    'savefig.bbox': 'tight',

    'grid.color': '#EAEAEA',
    'grid.linestyle': '-',
    'grid.linewidth': 0.7,

    'pdf.fonttype': 42,

    "figure.figsize": [single_column, single_column]
}

wide_plot_config = {
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'legend.fontsize': 7,
    'figure.titlesize': 7,

    'lines.linewidth': 1.5,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,

    'xtick.top': True,
    'ytick.right': True,
    'xtick.bottom': True,
    'ytick.left': True,

    'axes.spines.top': True,
    'axes.spines.right': True,
    'axes.spines.bottom': True,
    'axes.spines.left': True,

    'xtick.direction': 'inout',
    'ytick.direction': 'inout',

    'lines.markersize': 2,
    'font.sans-serif': ['Roboto'],
    'axes.labelcolor': 'black',

    'axes.edgecolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',

    'legend.frameon': False,
    'legend.framealpha': 0,

    'savefig.dpi': 300,
    'savefig.bbox': 'tight',

    'grid.color': '#EAEAEA',
    'grid.linestyle':  '-',
    'grid.linewidth': 0.7,

    'pdf.fonttype': 42,

    "figure.figsize": [single_column_large, single_column_large/1.333333]

}

pmatrix_plot_config = {
    'xtick.labelsize': 3,
    'ytick.labelsize': 3,
    'axes.titlesize': 3,
    'axes.labelsize': 3,
    'legend.fontsize': 3,
    'figure.titlesize': 3,

    'lines.linewidth': 0.6,
    'axes.linewidth': 0.2,
    'xtick.major.width': 0.2,
    'ytick.major.width': 0.2,
    'xtick.major.size': 0.2,
    'ytick.major.size': 0.2,

    'xtick.top': True,
    'ytick.right': True,
    'xtick.bottom': True,
    'ytick.left': True,

    'axes.spines.top': True,
    'axes.spines.right': True,
    'axes.spines.bottom': True,
    'axes.spines.left': True,

    'xtick.direction': 'inout',
    'ytick.direction': 'inout',

    'lines.markersize': 0.5,
    'font.sans-serif': ['Roboto'],
    'axes.labelcolor': 'black',

    'axes.edgecolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',

    'legend.frameon': False,
    'legend.framealpha': 0,

    'savefig.dpi': 300,
    'savefig.bbox': 'tight',

    'grid.color': '#EAEAEA',
    'grid.linestyle': '-',
    'grid.linewidth': 0.7,

    'pdf.fonttype': 42,

    "figure.figsize": [single_column, single_column]
}

scheme_plot_config = {
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'axes.titlesize': 5,
    'axes.labelsize': 5,
    'legend.fontsize': 5,
    'figure.titlesize': 5,

    'lines.linewidth': 1.2,
    'axes.linewidth': 2.4,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2,
    'ytick.major.size': 2,

    'xtick.top': True,
    'ytick.right': True,
    'xtick.bottom': True,
    'ytick.left': True,

    'axes.spines.top': True,
    'axes.spines.right': True,
    'axes.spines.bottom': True,
    'axes.spines.left': True,

    'xtick.direction': 'inout',
    'ytick.direction': 'inout',

    'lines.markersize': 2,
    'font.sans-serif': ['Roboto'],
    'axes.labelcolor': 'black',

    'axes.edgecolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',

    'legend.frameon': False,
    'legend.framealpha': 0,

    'savefig.dpi': 300,
    'savefig.bbox': 'tight',

    'grid.color': '#EAEAEA',
    'grid.linestyle': '-',
    'grid.linewidth': 0.7,

    'pdf.fonttype': 42,

    "figure.figsize": [single_column_scheme, single_column_scheme]

}


def apply_plot_config(plot_config):
    plt.rcParams.update(plot_config)


def init_roboto_font(path_to_roboto_folder=None):
    import matplotlib.font_manager as fm
    from pathlib import Path

    cwd = Path(__file__).parent
    if not path_to_roboto_folder:
        path_to_roboto_folder = cwd / "assets" / "roboto"

    font_dirs = [path_to_roboto_folder]
    font_files = fm.findSystemFonts(fontpaths=font_dirs)

    for font_file in font_files:
        fm.fontManager.addfont(font_file)


def rgb2hex(r, g, b):
    return f"#{r:02x}{g:02x}{b:02x}"


def auto_init():
    apply_plot_config(plot_config=default_plot_config)
    init_roboto_font()
    return 0


def convert_to_rgb(outputfile):
    img = Image.open(outputfile)
    print(f"File {outputfile} mode BEFORE: {img.mode}")

    print("Convert to RGB")
    img = img.convert("RGB")
    img.save(outputfile)

    img = Image.open(outputfile)
    print(f"File {outputfile} mode AFTER: {img.mode}")

    mode_after = img.mode
    return mode_after


if __name__ == "__main__":
    pass
