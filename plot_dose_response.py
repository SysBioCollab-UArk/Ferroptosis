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

plt.show()
