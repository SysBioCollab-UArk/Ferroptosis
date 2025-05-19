from util import plot_drc

basepath = 'RESULTS'
directories = ['HT1080', 'U2OS']
run_pydream_filename = 'run_ferroptosis_pydream.py'
expts_doses = [
    # experiment 0
    {'ht1080_erastin': [0, 0.056755444, 0.113703658, 0.228802973, 0.455863303],
      'xlabel': r'Erastin ($\mu$M)',
      'title': 'HT-1080'},
    # experiment 1
    {'u2os_erastin': [0, 0.216165686, 0.441596578, 0.896487612, 1.80157234],
     'xlabel': r'Erastin ($\mu$M)',
     'title': 'U2-OS'}
]
label_dict = {'GSH_Obs': 'GSH level', 'ht1080_erastin': None, 'u2os_erastin': None}

plot_drc(basepath, directories, run_pydream_filename, expts_doses, label_dict)
