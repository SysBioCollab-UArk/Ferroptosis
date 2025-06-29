from ferroptosis import model
from param_calibration import *
from pysb.simulator import ScipyOdeSimulator
from SIM_PROTOCOLS.sim_protocols import *

solver = ScipyOdeSimulator(model)

# concentrations in uM
# HT-1080 cells
# time_perturb_value = [{0: ("Erastin(sys_xc=None)", 0)},  # control expt
#                       {0: ("Erastin(sys_xc=None)", 0.056755444)},  # expt 1
#                       {0: ("Erastin(sys_xc=None)", 0.113703658)},  # expt 2
#                       {0: ("Erastin(sys_xc=None)", 0.228802973)},  # expt 3
#                       {0: ("Erastin(sys_xc=None)", 0.455863303)}]  # expt 4

# U2-OS cells
time_perturb_value = [{0: ("Erastin(sys_xc=None)", 0)},  # control expt
                      {0: ("Erastin(sys_xc=None)", 0.216165686)},  # expt 1
                      {0: ("Erastin(sys_xc=None)", 0.441596578)},  # expt 2
                      {0: ("Erastin(sys_xc=None)", 0.896487612)},  # expt 3
                      {0: ("Erastin(sys_xc=None)", 1.80157234)}]  # expt 4

scale_by_eidx_time = {"GSH_Obs": {"eidx": 0, "time": 300}}  # scale by output at t=300 min in expt 0
multi_exp_injection = ParallelExperiments(solver, t_equil=1e3, time_perturb_value=time_perturb_value,
                                          scale_by_eidx_time=scale_by_eidx_time)

custom_priors = {}  # 'N': ('uniform', 0.3)}
no_sample = ["Glu_intra_0", "Cystine_extra_0", "Cys_0", "Glu_Cys_GCL_Product_0", "Gly_0", "GSH_0", "GSSG_0",
             "NADPplus_0", "LOOH_0", "LOH_0", "Lipid_metab_0", "LO_0", "Erastin_0", "Trx_TrxR_0",  "RSL3_0",
             "Iron_chelators_0", "Iron_storage_0", "Ferrostatin_0", 'kf_RSL3_GPX4', 'kr_RSL3_GPX4', 'kf_Ferrostatin_LO',
             'kr_Ferrostatin_LO', 'kf_Storage', 'kr_Storage', 'kf_Chelators', 'kr_Chelators', "kcat_Cys_TrxR",
             "km_Cys_TrxR"]

obs_labels = {'GSH_obs': 'GSH level'}

exp_data_file = os.path.join('DATA', 'u2_os_gsh.csv')  # 'ht_1080_gsh.csv')

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      multi_exp_injection,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, obs_labels=obs_labels, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True})
