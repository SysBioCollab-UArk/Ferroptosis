from ferroptosis import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import os

# read in experimental data
datafiles = ['ht_1080_gsh.csv', 'u2_os_gsh.csv']
path = 'DATA'

# erastin concentrations (uM)
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
        plt.figure('%s_%s' % (datafile, obs), constrained_layout=True)
        yvals = [d['average'] for d in data if d['observable'] == obs]
        stderr = [d['stderr'] for d in data if d['observable'] == obs]
        plt.errorbar(erastin_conc[datafile], yvals, yerr=stderr, fmt='o', ms=8, capsize=6, label=label)
        plt.xlim(right=right)
        plt.ylim(bottom=0)
        plt.xlabel('Erastin (uM)')
        plt.ylabel('GSH (fraction)')
        plt.legend(loc='best')
    print()

# run simulations
sim = ScipyOdeSimulator(model, verbose=True)

# equilibration simulation
tspan = np.linspace(0, 1e3, 1001)
output = sim.run(tspan=tspan)

plt.figure(constrained_layout=True, figsize=(6.4*1.5, 4.8))
for i in range(len(model.species)):
    plt.plot(tspan, output.all['__s%d' % i], lw=2, label='s%d' % i)
plt.yscale('log')
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc="upper left", bbox_to_anchor=(1.02, 1), ncols=2)

# perturbation simulations
initials = output.species[-1]  # species concentrations to start next simulation
erastin_idx = [str(sp) for sp in model.species].index("Erastin(sys_xc=None)")  # index for erastin

gsh_vals = [output.observables['GSH_Obs'][-1]]  # initialize w/ conc of GSH for erastin=0 (to scale by)

plt.figure(constrained_layout=True)
tspan = np.linspace(0, 24*60, 1001)
for conc in erastin_conc['ht_1080_gsh.csv'][1:]:
    conc *= 1e6  # uM -> M concentration
    initials[erastin_idx] = conc
    output = sim.run(tspan=tspan, initials=initials)
    gsh_vals.append(output.observables['GSH_Obs'][-1])
    #####
    plt.plot(tspan, output.observables['GSH_Obs'], lw=2, label='conc = %g' % conc)
plt.xlabel('time (min)')
plt.ylabel('conc (M)')
plt.legend(loc='best')

# Plot simulated GSH levels to compare to data
gsh_vals = np.array(gsh_vals)
plt.figure('%s_%s' % ('ht_1080_gsh.csv', 'GSH_Obs'))
plt.plot(erastin_conc['ht_1080_gsh.csv'], gsh_vals / gsh_vals[0], '^', ms=10, label='simulation')
plt.legend(loc='best')

plt.show()
