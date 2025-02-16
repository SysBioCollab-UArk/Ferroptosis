import numpy as np
import matplotlib.pyplot as plt
import os

# read in experimental data
datafiles = ['ht_1080_gsh.csv', 'u2_os_gsh.csv']
path = 'DATA'

erastin_conc = {'ht_1080_gsh.csv': [0, 0.056755444, 0.113703658, 0.228802973, 0.455863303],
                'u2_os_gsh.csv': [0, 0.216165686, 0.441596578, 0.896487612, 1.80157234]}

for datafile in datafiles:
    if datafile == 'ht_1080_gsh.csv':
        label = 'HT-1080'
        right = 0.5
    else:
        label = 'U-2 OS'
        right = 2.0
    data = np.genfromtxt(os.path.join(path, datafile), dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    print(datafile)
    print(data.dtype.names)
    observables = np.unique([d['observable'] for d in data])
    for obs in observables:
        plt.figure(constrained_layout=True)
        yvals = [d['average'] for d in data if d['observable'] == obs]
        stderr = [d['stderr'] for d in data if d['observable'] == obs]
        plt.errorbar(erastin_conc[datafile], yvals, yerr=stderr, fmt='o', ms=8, capsize=6, label=label)
        plt.xlim(right=right)
        plt.ylim(bottom=0)
        plt.xlabel('Erastin (uM)')
        plt.ylabel('GSH (fraction)')
        plt.legend(loc='best')
    print()

plt.show()
