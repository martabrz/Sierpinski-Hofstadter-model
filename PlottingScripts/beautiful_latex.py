import matplotlib as mpl
import matplotlib.pyplot as plt

def latex_plot(scale = 1, fontsize = 12):
    """ Changes the size of a figure and fonts for the publication-quality plots. """
    fig_width_pt = 246.0
    inches_per_pt = 1.0/72.27
    golden_mean = (np.sqrt(5.0)-1.0)/2.0
    fig_width = fig_width_pt*inches_per_pt*scale
    fig_height = fig_width*golden_mean
    fig_size = [fig_width, fig_height]
    eps_with_latex = {
        "pgf.texsystem": "pdflatex", "text.usetex": True, "font.family": "serif", "font.serif": [], "font.sans-serif": [], "font.monospace": [], "axes.labelsize": fontsize, "font.size": fontsize, "legend.fontsize": fontsize, "xtick.labelsize": 12, "ytick.labelsize": 12, "figure.figsize": fig_size
        }
    mpl.rcParams.update(eps_with_latex)
