import numpy as np
import matplotlib.pyplot as plt
import os
import importlib

# Helper function for getting a good x-axis upper limit for dose-response plots
def round_up_nice(x):
    if x == 0:
        return 1
    exponent = np.floor(np.log10(x))
    fraction = x / 10**exponent
    nice_fraction = np.select(
        [fraction <= 1, fraction <= 2, fraction <= 5],
        [1, 2, 5],
        default=10
    )
    return nice_fraction * 10 ** exponent


samples_ALL = []
for results_dir in ['HT1080', 'U2OS']:

    # Set the 'path' variable to the directory where the SIM_DATA.csv, run_ferroptosis_pydream.py, and expt data file are
    path = os.path.join('RESULTS', results_dir)  # os.getcwd()

    # import everything from run_ferroptosis_pydream.py file that's in the path
    run_pydream_file = os.path.join(path,'run_ferroptosis_pydream.py')
    import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
    module = importlib.import_module(import_string)  # Import the module
    # the following loop emulates `from run_ferroptosis_pydream import *`
    # for name in getattr(module, '__all__', dir(module)):
    #     if not name.startswith('_'):
    #         globals()[name] = getattr(module, name)

    # get the path to the experimental data file referenced in the run_ferroptosis_pydream.py file that's in the path
    if not os.path.isabs(module.exp_data_file):
        exp_data_file = os.path.normpath(os.path.join(path, module.exp_data_file))
    else:
        exp_data_file = os.path.normpath(module.exp_data_file)

    conc_erastin = [list(d.values())[0][0][1] for d in module.time_perturb_value]
    print(conc_erastin)

    expt_data = np.genfromtxt(exp_data_file, dtype=None, delimiter=',', names=True,
                              encoding="utf_8_sig")
    print(expt_data.dtype.names)

    # get strings from the expt datafile name to use for the figure name and legend labels
    cell_type_marker, _ = os.path.splitext(os.path.basename(exp_data_file))
    label = '-'.join(cell_type_marker.upper().split('_')[:2])
    # 'ht_1080_gsh' -> ['HT', '1080', 'GSH'] -> HT-1080

    sim_data = np.genfromtxt(os.path.join(path, 'SIM_DATA.csv'), dtype=None, delimiter=',', names=True,
                             encoding="utf_8_sig")
    print(sim_data.dtype.names)

    plt.figure(constrained_layout=True, figsize=(6.4*0.7, 4.8*0.9))
    # sim data
    yval_min = sim_data['yval_min'] * 100
    yval_max = sim_data['yval_max'] * 100
    p = plt.plot(conc_erastin, (yval_min + yval_max)/2, ls='--', label='%s (sim)' % label)
    plt.fill_between(conc_erastin, yval_min, yval_max, alpha=0.25, color=p[0].get_color(), label='x')
    # expt data
    avg = expt_data['average'] * 100
    stderr = expt_data['stderr'] * 100
    plt.errorbar(conc_erastin, avg, yerr=stderr, fmt='o', ms=8, capsize=6, color=p[0].get_color(),
                 label='%s (expt)' % label)

    plt.xlabel(r'Erastin ($\mu$M)', fontsize=16)
    plt.ylabel('GSH level (%)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim(right=round_up_nice(max(conc_erastin)))
    plt.ylim(bottom=0)
    plt.legend(loc='best')

    # merge line and fill_between legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    n_sims = len([label for label in labels if 'sim' in label])
    new_handles = [(handles[n], handles[n + 1]) for n in range(0, n_sims * 2 - 1, n_sims)] + list(handles[n_sims * 2:])
    new_labels = [labels[n] for n in range(0, n_sims * 2 - 1, n_sims)] + list(labels[n_sims * 2:])
    plt.legend(new_handles, new_labels, loc='best')

    plt.savefig('fig_PyDREAM_DRC_%s' % cell_type_marker)

    ##########
    # Plot parameter histograms
    from param_calibration import *
    from scipy import stats
    from matplotlib.lines import Line2D
    from plotting import *

    # Helper function for getting the optimal number of columns for the histogram figure
    def get_ncols(ndims):
        if not isinstance(ndims, int) or ndims <= 0:
            raise ValueError("'ndims' must be a positive integer")
        r1 = round(math.sqrt(ndims))  # round() returns an int
        r2 = math.ceil(math.sqrt(ndims))  # math.ceil() also returns an int
        while r1 * r2 >= ndims:
            r1 -= 1
            r2 += 1
        return min(r1 + 1, r2 - 1)  # the smaller of the two integers is the # of columns


    logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))
    samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))

    calibrator = ParameterCalibration(module.model,
                                      exp_data_file,
                                      module.multi_exp_injection,
                                      priors=module.custom_priors,
                                      no_sample=module.no_sample)

    _, samples, _ = calibrator.create_figures(logps_files, samples_files, obs_labels=module.obs_labels, show_plots=True,
                                              plot_ll_args={'cutoff': 2}, plot_pd_args={'sharex': 'all'},
                                              which_plots=2)

    samples_ALL.append(samples)

# Create figure with parameter histograms overlaid
n_params = len(calibrator.parameter_idxs)
labels = [calibrator.model.parameters[i].name for i in calibrator.parameter_idxs]
ncols = get_ncols(n_params)
nrows = math.ceil(n_params / ncols)
labelsize = 10 * max(1, (2/5 * np.ceil(nrows / 2)))
fontsize = 10 * max(1, (3/5 * np.ceil(nrows / 2)))
colors = sns.color_palette(n_colors=n_params)
fig = plt.figure(constrained_layout=True, figsize=(0.65 * ncols * 6.4, 0.5 * nrows * 4.8))
axes = []
reference_ax = None
for n in range(n_params):
    print(n, end=' ')
    # share x-axis with first subplot
    share_x_with = None if n == 0 else reference_ax
    ax = fig.add_subplot(nrows, ncols, n + 1, sharex=share_x_with)
    if n == 0:
        reference_ax = ax
    axes.append(ax)
    for i, (samples, color) in enumerate(zip(samples_ALL, [colors[n], 'k'])):
        bw_adjust = 3.0 if i == 0 else 2.0
        sns.kdeplot(samples[:, n], color=color, fill=True, common_norm=False, ax=ax, bw_adjust=bw_adjust)
        x_vals = sorted(ax.collections[i].get_paths()[0].vertices[:, 0])  # get x-axis values from seaborn plot
        # get kernel density estimate (KDE) for calculating histogram distance and self distance
        kde = stats.gaussian_kde(samples[:, n])
        kde.set_bandwidth(kde.factor * bw_adjust)
        # TODO: save values of E_Dself and D, so we can plot them in descending order
        if i == 0:
            # calculate self distance (expected value)
            kde_ref = kde
            x_min_ref = x_vals[0]
            x_max_ref = x_vals[-1]
            E_Dself = calc_self_distance(kde_ref, len(samples[:, n]), x_min_ref, x_max_ref, 1000)
        else:
            # calculate histogram distance relative to the reference (use 2x the points, just to be safe)
            D = calc_hist_distance(kde, kde_ref, min(x_vals[0], x_min_ref), max(x_vals[-1], x_max_ref), 2000)
    empty_handle = Line2D([], [], linestyle="none")
    legend = ax.legend([empty_handle, empty_handle], ['D: %.3f' % D, r'E[D$^\mathrm{self}$]: %.3f' % E_Dself],
                       fontsize=0.9*labelsize, loc='best', handlelength=0, labelspacing=0.4)
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

plt.show()
