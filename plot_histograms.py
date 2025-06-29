from param_calibration import *
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import os
import importlib
import math
from util import get_fig_ncols


def calc_self_distance(kde, n_samples, x_min, x_max, num_x_pts):
    x_vals = np.linspace(x_min, x_max, num_x_pts)
    dx = np.array([x_vals[i] - x_vals[i - 1] for i in range(1, len(x_vals))])
    E_Dself = 0.5 * np.sqrt(2 / n_samples / np.pi) * np.sum(np.sqrt(kde(x_vals)[1:] * dx))
    return E_Dself


def calc_hist_distance(kde, kde_ref, x_min, x_max, num_x_pts):
    x_vals = np.linspace(x_min, x_max, num_x_pts)
    dx = np.array([x_vals[i] - x_vals[i - 1] for i in range(1, len(x_vals))])
    D = 0.5 * np.sum(np.abs(kde(x_vals)[1:] - kde_ref(x_vals)[1:]) * dx)
    return D


def write_multicolor_word(fig, x, y, word, colors, fontsize=10, fontweight='normal', additional_text=None):

    if len(colors) != len(word):
        raise Exception("Length of 'colors' (%d) does not match length of 'word' (%d)" % (len(colors), len(word)))

    # Font properties
    fontprops = FontProperties(size=fontsize, weight=fontweight)

    # Measure full word width (includes kerning)
    tp_word = TextPath((0, 0), word, prop=fontprops)
    word_width = tp_word.get_extents().width / 72 / fig.get_size_inches()[0]

    # Measure individual letter widths
    letter_widths = []
    for letter in word:
        tp_letter = TextPath((0, 0), letter, prop=fontprops)
        letter_widths.append(tp_letter.get_extents().width / 72 / fig.get_size_inches()[0])  # convert pts to figure coords

    # Estimate average inter-letter spacing from kerning
    letter_spacing = 0 if len(word) <= 1 else (word_width - sum(letter_widths)) / (len(word) - 1)

    # Draw each character
    for i, (letter, color) in enumerate(zip(word, colors)):
        fig.text(x, y, letter, color=color, ha='left', fontsize=fontsize, fontweight=fontweight)
        # Move x-position forward by width of character + additional spacing
        x += letter_widths[i] + letter_spacing

    # Add additional text, if exists
    if isinstance(additional_text, str):
        fig.text(x, y, additional_text, ha='left', fontsize=fontsize, fontweight=fontweight)
    elif isinstance(additional_text, dict):
        text = additional_text.pop('text')
        ha = additional_text.pop('ha', 'left')
        fs = additional_text.pop('fontsize', fontsize)
        fw = additional_text.pop('fontweight', fontweight)
        fig.text(x, y, text, ha=ha, fontsize=fs, fontweight=fw, **additional_text)


if __name__ == '__main__':

    samples_ALL = []
    cell_line_dirs = ['HT1080', 'U2OS']
    bw_adjust = [3.0, 2.0]
    for cell_line_dir, bw_adj in zip(cell_line_dirs, bw_adjust):

        # Set the 'path' variable to the directory where the SIM_DATA.csv, run_ferroptosis_pydream.py, and expt data file are
        path = os.path.join('RESULTS', cell_line_dir)  # os.getcwd()

        # import everything from run_ferroptosis_pydream.py file that's in the path
        run_pydream_file = os.path.join(path, 'run_ferroptosis_pydream.py')
        import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
        module = importlib.import_module(import_string)  # Import the module

        # get the path to the experimental data file referenced in the run_ferroptosis_pydream.py file that's in the path
        if not os.path.isabs(module.exp_data_file):
            exp_data_file = os.path.normpath(os.path.join(path, module.exp_data_file))
        else:
            exp_data_file = os.path.normpath(module.exp_data_file)

        logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))
        samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))

        calibrator = ParameterCalibration(module.model,
                                          exp_data_file,
                                          module.multi_exp_injection,
                                          priors=module.custom_priors,
                                          no_sample=module.no_sample)

        _, samples, _ = calibrator.create_figures(
            logps_files, samples_files, obs_labels=module.obs_labels, show_plots=True,
            plot_ll_args={'cutoff': 2, 'save_plot': 'fig_PyDREAM_log_ps_%s' % cell_line_dir},
            plot_pd_args={'sharex': 'all', 'bw_adjust': bw_adj},
            which_plots=2)

        samples_ALL.append(samples)

    # Create figure with parameter histograms overlaid
    n_params = len(calibrator.parameter_idxs)
    ncols = get_fig_ncols(n_params)
    nrows = math.ceil(n_params / ncols)
    fig = plt.figure(constrained_layout=True, figsize=(0.65 * ncols * 6.4, 0.5 * nrows * 4.8))
    #
    labels = [calibrator.model.parameters[i].name for i in calibrator.parameter_idxs]
    labelsize = 10 * max(1, (2 / 5 * np.ceil(nrows / 2)))
    fontsize = 10 * max(1, (3 / 5 * np.ceil(nrows / 2)))
    colors = sns.color_palette(n_colors=n_params)
    axes = []
    reference_ax = None
    #
    D = []  # histogram distances
    E_Dself = []  # self distances
    for n in range(n_params):
        print(n, end=' ')
        # share x-axis with first subplot
        share_x_with = None if n == 0 else reference_ax
        ax = fig.add_subplot(nrows, ncols, n + 1, sharex=share_x_with)
        if n == 0:
            reference_ax = ax
        axes.append(ax)
        for i, (samples, color) in enumerate(zip(samples_ALL, [colors[n], 'k'])):
            sns.kdeplot(samples[:, n], color=color, fill=True, common_norm=False, ax=ax, bw_adjust=bw_adjust[i])
            x_vals = sorted(ax.collections[i].get_paths()[0].vertices[:, 0])  # get x-axis values from seaborn plot
            # get kernel density estimate (KDE) for calculating histogram distance and self distance
            kde = stats.gaussian_kde(samples[:, n])
            kde.set_bandwidth(kde.factor * bw_adjust[i])
            if i == 0:
                # calculate self distance (expected value)
                kde_ref = kde
                x_min_ref = x_vals[0]
                x_max_ref = x_vals[-1]
                E_Dself.append(calc_self_distance(kde_ref, len(samples[:, n]), x_min_ref, x_max_ref, 1000))
            else:
                # calculate histogram distance relative to the reference (use 2x the points, just to be safe)
                D.append(calc_hist_distance(kde, kde_ref, min(x_vals[0], x_min_ref), max(x_vals[-1], x_max_ref), 2000))
        empty_handle = Line2D([], [], linestyle="none")
        legend = ax.legend([empty_handle, empty_handle],
                           ['D: %.3f' % D[-1], r'E[D$^\mathrm{self}$]: %.3f' % E_Dself[-1]], fontsize=0.9 * labelsize,
                           loc='best', handlelength=0, labelspacing=0.4)
        legend.set_frame_on(False)
        ax.set_yticklabels([])
        ax.set_ylabel(None)
        ax.tick_params(axis='x', labelsize=labelsize)
        ax.label_outer()
        ax.set_title(labels[n], fontsize=labelsize)
    fig.supxlabel(r'log$_{10}$ value', fontsize=fontsize)
    fig.supylabel('Density', fontsize=fontsize)

    # Create a common figure legend
    additional_text = {'text': ": HT-1080", 'fontweight': 'normal'}
    write_multicolor_word(fig, 0.82, 0.1, "Multicolor", sns.color_palette(n_colors=len("Multicolor")),
                          fontsize=fontsize, fontweight='bold', additional_text=additional_text)
    additional_text = {'text': ": U2-OS", 'fontweight': 'normal'}
    write_multicolor_word(fig, 0.82, 0.078, "Black", ['k'] * len("Black"),
                          fontsize=fontsize, fontweight='bold', additional_text=additional_text)

    # Plot histogram distances with self distances
    NCOLS = 2
    sorted_idxs = np.argsort(D)[::-1]  # sort from largest to smallest
    sorted_labels = [labels[i] for i in sorted_idxs]
    table_data = []
    for col in range(NCOLS):
        start = col * n_params // NCOLS
        end = (col + 1) * n_params // NCOLS
        print('start: %d, end: %d, len(sorted_labels): %d' % (start, end, len(sorted_labels)))
        for row in range(end - start):
            print('row: %d, len(table_data): %d' % (row, len(table_data)))
            if row == len(table_data):
                table_data.append([])
            if len(table_data[row]) < col:
                table_data[row] += [""]
            table_data[row] += ["%d. %s" % (start + row, sorted_labels[start + row])]
    for row in table_data:
        print(row)

    # save overlaid histograms figure
    plt.savefig('fig_PyDREAM_hist_overlay_%s' % str.join('_', cell_line_dirs))

    plt.figure(constrained_layout=True)
    plt.bar(np.arange(len(sorted_idxs)), [D[i] for i in sorted_idxs])
    for i, e in enumerate([E_Dself[j] for j in sorted_idxs]):
        plt.plot([i - 0.4, i + 0.4], [e, e], color='r', lw=2)
    plt.xlabel('Index')
    plt.ylabel('Histogram Distance')

    # Add table in the corner
    table = plt.table(
        cellText=table_data,
        colLabels=None,
        loc='upper right',
        cellLoc='left'
    )
    table.set_fontsize(8)  # Set font size
    table.scale(1, 0.8)  # Stretch horizontally and vertically
    table.auto_set_column_width([i for i in range(NCOLS)])  # Adjust column widths automatically
    # Remove all borders
    for cell in table.get_celld().values():
        cell.set_linewidth(0)

    # save histogram distances bar plot
    plt.savefig('fig_PyDREAM_hist_distances_%s' % str.join('_', cell_line_dirs))

    plt.show()
