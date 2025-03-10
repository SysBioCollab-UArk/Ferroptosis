from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A, pi

Model()

# Volume
Diameter = 15  # micrometers https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079008/
Volume = 4/3 * pi * (Diameter / 1e6 / 2)**3  # m^3

# Metabolites
Monomer('Glu', ['loc'], {'loc': ['intra', 'extra']})
Monomer('Cystine', ['loc'], {'loc': ['intra', 'extra']})
Monomer("Cys")
Monomer("Glu_Cys_GCL_Product")
Monomer("Gly")
Monomer("GSH")
Monomer("GSSG")
Monomer("NADPplus")
Monomer("NADPH")
Monomer("LOOH")
Monomer("LOH")
Monomer("Lipid_metab")
Monomer("LO", ["ferrostatin"])

Parameter("Glu_intra_0", 0)  # 0.1109/1e3 * Volume * N_A) #0.1109 is in mM
Parameter("Cystine_extra_0", 0)  # 0.0009/1e3 * Volume * N_A * 1000) #0.0009 is in mM
Parameter("Cys_0", 0)
Parameter("Glu_Cys_GCL_Product_0", 0)
Parameter("Gly_0", 0)  # originally 401805906
Parameter("GSH_0", 0)  # 1e5)
Parameter("GSSG_0", 0)  # 1e3)  # originally 1e5, changed to fix scale of Gly,GSH,GSSG graph
Parameter("NADPplus_0", 0)
Parameter("NADPH_0", 1e4)  # 500 3197682) # TODO: this must be non-zero (may want to add a source term)
Parameter("LOOH_0", 0)  # 1000 532947)
Parameter("LOH_0", 0)  # 50
Parameter("Lipid_metab_0", 0)  # 500 1)
Parameter("LO_0", 0)

Initial(Glu(loc='intra'), Glu_intra_0)
Initial(Cystine(loc='extra'), Cystine_extra_0)
Initial(Cys(), Cys_0)
Initial(Glu_Cys_GCL_Product(), Glu_Cys_GCL_Product_0)
Initial(Gly(), Gly_0)
Initial(GSH(), GSH_0)
Initial(GSSG(), GSSG_0)
Initial(NADPplus(), NADPplus_0)
Initial(NADPH(), NADPH_0)
Initial(LOOH(), LOOH_0)
Initial(LOH(), LOH_0)
Initial(Lipid_metab(), Lipid_metab_0)
Initial(LO(ferrostatin=None), LO_0)

# Enzymes
Monomer("System_Xc", ["erastin"])
Monomer("GCL")
Monomer("GSS")
Monomer("GPX4", ["rsl3"])
Monomer("GR")
Monomer("G6PD")
Monomer("PPARg")
Monomer("Iron", ["b"])

Parameter("System_Xc_0", 100)  # 50*42 - estimating 50 ppm times 42 million proteins/cell
Parameter("GCL_0", 100)  # origianlly 1918609200
Parameter("GSS_0", 1e3)  # 250 15 is in mM, originally 15/1e3 * Volume * N_A
Parameter("GPX4_0", 2e6)  # CHANGE THIS BACK5/1e3 * Volume * N_A) #5 is in mM
Parameter("GR_0", 2e6)  # 9593046)
Parameter("G6PD_0", 2e6)  # 19719039)
Parameter("PPARg_0", 100)  # 1*42)  # estimating 1 ppm times 42 million proteins/cell
Parameter("Iron_0", 100)  # 0.61/1e6 * Volume * N_A)  # 0.61 is in micromolar

Initial(System_Xc(erastin=None), System_Xc_0)
Initial(GCL(), GCL_0)
Initial(GSS(), GSS_0)
Initial(GPX4(rsl3=None), GPX4_0)
Initial(GR(), GR_0)
Initial(G6PD(), G6PD_0)
Initial(PPARg(), PPARg_0)
Initial(Iron(b=None), Iron_0)

# Inhibitors
Monomer("Erastin", ["sys_xc"])
Monomer("Trx_TrxR")
Monomer("RSL3", ["gpx4"])
Monomer("Iron_chelators", ["iron"])
Monomer("Iron_storage", ["iron"])
Monomer("Ferrostatin", ["lo"])

Parameter("Erastin_0", 0)  # 1e9
Parameter("Trx_TrxR_0", 0)  # 100 originally 0
Parameter("RSL3_0", 0)  # 6e6)  # 2e6 real value
Parameter("Iron_chelators_0", 0)
Parameter("Iron_storage_0", 0)  # 10
Parameter("Ferrostatin_0", 0)

Initial(Erastin(sys_xc=None), Erastin_0)
Initial(Trx_TrxR(), Trx_TrxR_0)
Initial(RSL3(gpx4=None), RSL3_0)
Initial(Iron_chelators(iron=None), Iron_chelators_0)
Initial(Iron_storage(iron=None), Iron_storage_0)
Initial(Ferrostatin(lo=None), Ferrostatin_0)

# Observables
Observable("Cys_Obs", Cys())
Observable("Glu_Cys_GCL_Product_Obs", Glu_Cys_GCL_Product())
Observable("Gly_Obs", Gly())
Observable("GSH_Obs", GSH())
Observable("LOOH_Obs", LOOH())
Observable("NADPH_Obs", NADPH())
Observable("GSSG_Obs", GSSG())
Observable("NADPplus_Obs", NADPplus())
Observable("Lipid_metab_Obs", Lipid_metab())
Observable("LO_Obs", LO(ferrostatin=None))

'''Observable("Cystine_extra_Obs", Cystine(loc='extra')) # Cystine_extra())
Observable("Cystine_intra_Obs", Cystine(loc='intra')) # Cystine_intra())
Observable("Glu_intra_Obs", Glu(loc='intra')) # Glu_intra())
Observable("Glu_extra_Obs", Glu(loc='extra')) # Glu_extra())
Observable("LOH_Obs", LOH())
Observable("GPX_4_free", GPX4(rsl3=None))
Observable("NADPtotal", NADPH()+NADPplus())
Observable("Irontotal", Iron())'''

# Rules
# Glu (intracellular) + Cystine (extracellular) + System Xc-->
# Glu (extracellular) + Cystine (intracellular) + System Xc-

# Parameter("k_Glu_intra", 2e3)
# Rule("Glu_Intra_Synthesis", None >> Glu_intra(), k_Glu_intra)
# Rule("Glu_Intra_Synthesis", None >> Glu(loc='intra'), k_Glu_intra)
Parameter("k_Glu_intra_syn", 2e3)
Parameter("k_Glu_intra_deg", 1e-3)
Rule("Glu_Intra_Syn_Deg", None | Glu(loc='intra'), k_Glu_intra_syn, k_Glu_intra_deg)

Parameter("k_Glu_extra", 0.1)  # 5 originally
# Rule("Glu_Extra_Degradation", Glu_extra() >> None, k_Glu_extra)
Rule("Glu_Extra_Degradation", Glu(loc='extra') >> None, k_Glu_extra)

# Parameter("k_Cystine_extra", 1e3)
# Rule("Cystine_Extra_Synthesis", None >> Cystine_extra(), k_Cystine_extra)
# Rule("Cystine_Extra_Synthesis", None >> Cystine(loc='extra'), k_Cystine_extra)
Parameter("k_Cystine_extra_syn", 1e3)
Parameter("k_Cystine_extra_deg", 1e-3)
Rule("Cystine_Extra_Syn_Deg", None | Cystine(loc='extra'), k_Cystine_extra_syn, k_Cystine_extra_deg)

Parameter("kcat_Glu_Cystine", 1e-7)
# Parameter("km_Glu_Cystine", 1e7)
# Expression("keff_Glu_Cystine", kcat_Glu_Cystine / (km_Glu_Cystine+Glu_intra_Obs+Cystine_extra_Obs))
# Rule("Glu_Cystine_Transport",
#      Glu_intra() + Cystine_extra() + System_Xc(erastin=None) >>
#      Glu_extra() + Cystine_intra() + System_Xc(erastin=None), kcat_Glu_Cystine)
Rule("Glu_Cystine_Transport",
     Glu(loc='intra') + Cystine(loc='extra') + System_Xc(erastin=None) >>
     Glu(loc='extra') + Cystine(loc='intra') + System_Xc(erastin=None), kcat_Glu_Cystine)

# Cystine (intracellular)  -> Cys
Parameter("k_Cystine_intra", 0.1)
# Rule("Cystine_intra_to_cys", Cystine_intra() >> Cys(), k_Cystine_intra)
Rule("Cystine_intra_to_cys", Cystine(loc='intra') >> Cys(), k_Cystine_intra)

# Cys  + Trx/TrxR -> ? + Trx/TrxR  # Note: Giving value here is important to stabilize
Parameter("kcat_Cys_TrxR", 1)
Parameter("km_Cys_TrxR", 100)
Expression("keff_Cys_TrxR", kcat_Cys_TrxR / (km_Cys_TrxR + Cys_Obs))
Rule("Cys_deg_TrxR", Cys() + Trx_TrxR() >> Trx_TrxR(), keff_Cys_TrxR)  # kcat_Cys_TrxR)

# Glu + Cys + GCL -> Glu_Cys_GCL_Product + GCL
Parameter("kcat_Glu_Cys_GCL", 2e-7)
# Parameter("km_Glu_Cys_GCL", 1e7) # 100
# Expression("keff_Glu_Cys_GCL", kcat_Glu_Cys_GCL / (km_Glu_Cys_GCL+Cys_Obs+Glu_intra_Obs))
# Rule("Glu_Cys_GCL", Glu_intra() + Cys() + GCL() >> Glu_Cys_GCL_Product() + GCL(), kcat_Glu_Cys_GCL)
Rule("Glu_Cys_GCL", Glu(loc='intra') + Cys() + GCL() >> Glu_Cys_GCL_Product() + GCL(), kcat_Glu_Cys_GCL)

# Glu_Cys_GCL_Product + Gly + GSS -> GSH + GSS
Parameter("kcat_Gly_Glu_Cys_GCL_Product", 7e-5)  # 2e-7, 1e-4
Parameter("km_Gly_Glu_Cys_GCL_Product", 100)  # 100
Expression("keff_Gly_Glu_Cys_GCL_Product",
           kcat_Gly_Glu_Cys_GCL_Product / (km_Gly_Glu_Cys_GCL_Product + Gly_Obs + Glu_Cys_GCL_Product_Obs))
Rule("Gly_Glu_Cys_GCL_Product", Glu_Cys_GCL_Product() + Gly() + GSS() >> GSH() + GSS(),
     keff_Gly_Glu_Cys_GCL_Product)

# Glu_Cys_GCL_Product -> None
Parameter('kdeg_Glu_Cys_GCL_Prod', 100)
Rule('Glu_Cys_GCL_Product_degradation', Glu_Cys_GCL_Product() >> None, kdeg_Glu_Cys_GCL_Prod)

Parameter("k_Gly_synth", 1e3)  # 20 100
Parameter('k_Gly_deg', 100)
Rule("Gly_Synth_Deg", None | Gly(), k_Gly_synth, k_Gly_deg)

# LOOH + GSH + GPX4 -> LOH + GSSG + GPX4
Parameter("kcat_LOOH_GSH", 2.3e-7)
Parameter("km_LOOH_GSH", 100)  # 100
Expression("keff_LOOH_GSH", kcat_LOOH_GSH / (km_LOOH_GSH + LOOH_Obs + GSH_Obs))
Rule("LOOH_LOH", LOOH() + GSH() + GPX4(rsl3=None) >> LOH() + GSSG() + GPX4(rsl3=None), keff_LOOH_GSH)

# GSH degradation
Parameter('k_deg_GSH', 0.05)
Rule('deg_GSH', GSH() >> None, k_deg_GSH)

# LOOH synthesis
Parameter('k_synth_LOOH', 1000)
Rule('synth_LOOH', None >> LOOH(), k_synth_LOOH)

# LOOH degradation
Parameter('k_deg_LOOH', 0.05)
Rule('deg_LOOH', LOOH() >> None, k_deg_LOOH)

# LOH degradation
Parameter('k_deg_LOH', 0.08)
Rule('deg_LOH', LOH() >> None, k_deg_LOH)

# PPARG -> PPARG + Lipid Metabolism
Parameter("k_PPARG_Lipid_Metab", 3)  # 1
Rule("PPARG_Lipid_metab", PPARg() >> PPARg() + Lipid_metab(), k_PPARG_Lipid_Metab)

# GSSG + NADPH + GR -> GSH + NADP+ + GR
Parameter("kcat_GSSG_NADPH", 5e-6)  # 5e-7
Parameter("km_GSSG_NADPH", 100)  # 1e9
Expression("keff_GSSG_NADPH", kcat_GSSG_NADPH / (km_GSSG_NADPH + GSSG_Obs + NADPH_Obs))
Rule("GSSG_GSH_NADPH", GSSG() + NADPH() + GR() >> GSH() + NADPplus() + GR(), keff_GSSG_NADPH)

# NADP+ + G6PD -> NADPH + G6PD (only if one compound is required for reaction to occur)
Parameter("kcat_NADPplus_NADPH", 1e-4)  # 1e-3)
Parameter("km_NADPplus_NADPH", 100)  # 100
Expression("keff_NADPplus_NADPH", kcat_NADPplus_NADPH / (km_NADPplus_NADPH + NADPplus_Obs))
Rule("NADPplus_NADPH", NADPplus() + G6PD() >> NADPH() + G6PD(), keff_NADPplus_NADPH)

# LOOH + Lipid_metab + Iron -> LO. + Iron
Parameter("kcat_LOOH_Lipid_metab", 1e-2)  # 1e3
Parameter("km_LOOH_Lipid_metab", 100)  # 100
Expression("keff_LOOH_Lipid_metab", kcat_LOOH_Lipid_metab / (km_LOOH_Lipid_metab + LOOH_Obs + Lipid_metab_Obs))
Rule("LOOH_Lipid_metab", LOOH() + Lipid_metab() + Iron(b=None) >> LO(ferrostatin=None) + Iron(b=None),
     keff_LOOH_Lipid_metab)

# Lipid_metab degradation
Parameter('k_deg_Lipid_metab', 0.04) #0.04
Rule('deg_Lipid_metab', Lipid_metab() >> None, k_deg_Lipid_metab)

# LO degradation
Parameter('k_deg_LO', 0.5)
Rule('deg_LO', LO() >> None, k_deg_LO)

# System Xc- + Erastin   <--> System Xc- : Erastin
Parameter("kf_Xc_Erastin", 1)  # 1
Parameter("kr_Xc_Erastin", 100)  # 1
Rule("Xc_Erastin", System_Xc(erastin=None) + Erastin(sys_xc=None) | System_Xc(erastin=1) % Erastin(sys_xc=1),
     kf_Xc_Erastin, kr_Xc_Erastin)

# RSL3 + GPX4 <--> RSL3:GPX4
Parameter("kf_RSL3_GPX4", 1)  # 1
Parameter("kr_RSL3_GPX4", 100)  # 100
Rule("RSL3_binds_GPX4", RSL3(gpx4=None) + GPX4(rsl3=None) | RSL3(gpx4=1) % GPX4(rsl3=1), kf_RSL3_GPX4, kr_RSL3_GPX4)

# Ferrostatin-1 + LO. <--> Ferrostatin-1:LO.
Parameter("kf_Ferrostatin_LO", 1)  # 1
Parameter("kr_Ferrostatin_LO", 100)  # 1
Rule("Ferrostatin_LO", Ferrostatin(lo=None) + LO(ferrostatin=None) | Ferrostatin(lo=1) % LO(ferrostatin=1),
     kf_Ferrostatin_LO, kr_Ferrostatin_LO)

# Iron storage & trafficking heme metabolism + Iron <--> Iron storage & trafficking heme metabolism : Iron
Parameter("kf_Storage", 10)  # 1
Parameter("kr_Storage", 100)  # 1
Rule("Storage_hememetabolism", Iron_storage(iron=None) + Iron(b=None) | Iron_storage(iron=1) % Iron(b=1),
     kf_Storage, kr_Storage)

# Iron Chelators + Iron  <--> Iron Chelators : Iron
Parameter("kf_Chelators", 1)  # 1
Parameter("kr_Chelators", 100)  # 1
Rule("Chelators_Iron", Iron_chelators(iron=None) + Iron(b=None) | Iron(b=1) % Iron_chelators(iron=1),
     kf_Chelators, kr_Chelators)

if __name__ == '__main__':

    # run simulation
    Tspan = np.linspace(0, 10000, 100001)
    sim = ScipyOdeSimulator(model, Tspan, verbose=False)

    Erastin_equil = []
    Erastin_conc = [0, 3e8, 5e8, 7e8, 1e9]
    for km, e0 in zip([100] * 5, Erastin_conc): # , 1e3, 1e4]:
        print('km_LOOH_GSH:', km)
        print('Erastin_0:', e0)
        result = sim.run(param_values={"km_LOOH_GSH": km, 'Erastin_0': e0})

        print(result.observables['GSH_Obs'][-1])
        Erastin_equil.append(result.observables['GSH_Obs'][-1])

        # plot results
        obs2plot = [#["Cystine_extra_Obs", "Cystine_intra_Obs", "Cys_Obs"],
                    ["Cystine_intra_Obs", "Cys_Obs"],
                    ["Glu_intra_Obs", "Glu_extra_Obs"],
                    # ["Gly_Obs","GSSG_Obs" ,"GSH_Obs"],
                    ["GSH_Obs"],
                    ["LOOH_Obs", "LOH_Obs"],
                    ["NADPH_Obs","NADPplus_Obs","NADPtotal"],
                    ["Irontotal","Lipid_metab_Obs","LO_Obs"]]

        fig, axs = plt.subplots(nrows=len(obs2plot), ncols=1, sharex='all', constrained_layout=True,
                                figsize = (9.6, 2.4*len(obs2plot)))
        # default figsize = (6.4,4.8)
        fig.suptitle("km_LOOH_GSH = %g, Erastin_0 = %g" % (km, e0), fontsize=16)
        fig.supxlabel("Time", fontsize=16)
        fig.supylabel("Concentration", fontsize=16)

        row = 0
        col = 0
        for obs_group in obs2plot:
            for obs_name in obs_group:
                axs[row].plot(Tspan, result.observables[obs_name], lw=2, label=obs_name)
                axs[row].legend(loc="best", fontsize=16)
                axs[row].xaxis.set_tick_params(labelsize=16)
                axs[row].yaxis.set_tick_params(labelsize=16)
            if (col+1) % 1 == 0:
                row += 1
                col = 0
            else:
                col += 1

        # plt.savefig("km_Glu_Cys_%g.pdf" % km,format="pdf")

    # plot equilibrium Erastin values
    plt.figure(constrained_layout=True)
    plt.plot(Erastin_conc, np.array(Erastin_equil) / Erastin_equil[0] * 100, 'ko', ms=8)
    plt.ylim(bottom=0)
    plt.xlabel("[Erastin]", fontsize=16)
    plt.ylabel("GSH level (%)", fontsize=16)

    plt.show()

    # Tspan=np.linspace(0,1000,10001)
    # print(Tspan)
    #
    # for km in [1e5,1e7,1e8,1e9]:
    #     sim=ScipyOdeSimulator(model,Tspan,verbose=True)
    #     result=sim.run(param_values={"km_LOOH_GSH": km})
    #
    #     #obs2plot=["LOOH_Obs","LOH_Obs","LO_Obs"]
    #     obs2plot=[["LO_Obs"],
    #               ["LOOH_Obs","LOH_Obs"],
    #               ["GSH_Obs","GSSG_Obs"],
    #               ["NADPplus_Obs","NADPH_Obs"],
    #               ["Cystine_extra_Obs","Cystine_intra_Obs","Cys_Obs"],
    #               ["Glu_intra_Obs"]]#,"Glu_extra_Obs"]]
    #
    #     fig,axs=plt.subplots(nrows=3,ncols=2,sharex=True,figsize=(9.6,9.6))
    #     fig.suptitle("km_LOOH_GSH = %g" % km, fontsize=15)
    #
    #     #use this code for plotting multiple km values
    #     #fig.suptitle("km_Glu_Cystine = %g,\n km_Glu_Cystine = %g,\n km_Glu_Cystine = %g" %
    #                  #(km_Glu_Cystine.value,km_Glu_Cystine.value,km_Glu_Cystine.value), fontsize=15)
    #
    #     row=0
    #     col=0
    #     for obs_group in obs2plot:
    #         #plt.figure()
    #         for obs_name in obs_group:
    #             #model.observables:
    #             #plt.figure()
    #             maxconc=LOOH_0.value
    #             #plt.plot(Tspan,result.observables[obs_name]/maxconc,lw=2,label=obs_name)
    #             axs[row,col].plot(Tspan, result.observables[obs_name], lw=2, label=obs_name)
    #             if row == 2:
    #                 axs[row,col].set_xlabel("Time")
    #             axs[row,col].set_ylabel("Concentration")
    #             axs[row,col].legend(loc="best")
    #
    #         if (col+1) % 2 == 0:
    #             row += 1
    #             col = 0
    #         else:
    #             col += 1
    #
    #             plt.tight_layout()
    #     plt.savefig("km_Glu_Cys_%g.pdf" % km,format="pdf")
    #
    #
    #
    # #Use this code if planning to save files with multiple km changes
    # #plt.savefig("km_1_%g,km_2_%g,km_3_%g.pdf" % (km_1,km_2,km_3),format="pdf")
    # plt.show()
