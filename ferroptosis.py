from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A, pi
#test commit
Model()

#Volume
Diameter = 15 #micrometers https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079008/
Volume = 4/3 * pi * (Diameter / 1e6 / 2)**3 #m^3



#Metabolites
Monomer("Glu_intra") #done
Monomer("Cystine_extra")
Monomer("Glu_extra") #done
Monomer("Cystine_intra")
Monomer("Cys") #done
Monomer("Gly")
Monomer("GSH") #done
Monomer("GSSG") #done
Monomer("NADPplus")
Monomer("NADPH")
Monomer("Iron_chelators",["iron"])
Monomer("LOOH")
Monomer("LOH")
Monomer("Ferrostatin",["lo"])
Monomer("LO",["ferrostatin"])
Monomer("Glu_Cys_GCL_Product")
Monomer("Lipid_metab")

#enzymes
Monomer("System_Xc",["erastin"])
Monomer("Trx_TrxR")
Monomer("GCL")
Monomer("GSS")
Monomer("GPX4",["rsl3"])
Monomer("GR")
Monomer("G6PD")
Monomer("PPARg") #Add into model!!!
Monomer("Iron_storage",["iron"])

#drugs
Monomer("Erastin",["sys_xc"])
Monomer("RSL3",["gpx4"])
Monomer("Iron",["b"])

#Paramaters
Parameter("Glu_intra_0",0.1109/1e3 * Volume * N_A) #0.1109 is in mM
Parameter("Cystine_extra_0",0.0009/1e3 * Volume * N_A * 1000) #0.0009 is in mM
Parameter("System_Xc_0",50*42) #estimating 50 ppm times 42 million proteins/cell
Parameter("Erastin_0",0)
Parameter("RSL3_0",6e6) #2e6 real value
Parameter("GPX4_0",0) # CHANGE THIS BACK5/1e3 * Volume * N_A) #5 is in mM
Parameter("Cys_0",0)
Parameter("Trx_TrxR_0",0)
Parameter("Iron_storage_0",0)
Parameter("Iron_0",0.61/1e6 * Volume * N_A) #0.61 is in micromolar
Parameter("Iron_chelators_0",0)
Parameter("Ferrostatin_0",0)
Parameter("LO_0",0)
Parameter("GCL_0",1918609200)
Parameter("Glu_Cys_GCL_Product_0",0)
Parameter("LOOH_0",532947)
Parameter("Lipid_metab_0",1)
Parameter("Gly_0",401805906)
Parameter("GSS_0",15/1e3 * Volume * N_A) #15 is in mM
Parameter("GSH_0",0)
Parameter("GSSG_0",0)
Parameter("NADPH_0",3197682)
Parameter("GR_0",9593046)
Parameter("NADPplus_0",0)
Parameter("G6PD_0",19719039)
Parameter("LOH_0",0)
Parameter("PPARg_0",1*42) #estimating 1 ppm times 42 million proteins/cell



#Initials
Initial(Glu_intra(),Glu_intra_0)
Initial(Cystine_extra(),Cystine_extra_0)
Initial(System_Xc(erastin=None),System_Xc_0)
Initial(Erastin(sys_xc=None),Erastin_0)
Initial(RSL3(gpx4=None),RSL3_0)
Initial(GPX4(rsl3=None),GPX4_0)
Initial(Cys(),Cys_0)
Initial(Trx_TrxR(),Trx_TrxR_0)
Initial(Iron_storage(iron=None),Iron_storage_0)
Initial(Iron(b=None),Iron_0)
Initial(Iron_chelators(iron=None),Iron_chelators_0)
Initial(Ferrostatin(lo=None),Ferrostatin_0)
Initial(PPARg(),PPARg_0)
Initial(LO(ferrostatin=None),LO_0)
Initial(GCL(),GCL_0)
Initial(Glu_Cys_GCL_Product(),Glu_Cys_GCL_Product_0)
Initial(LOOH(),LOOH_0)
Initial(Lipid_metab(),Lipid_metab_0)
Initial(Gly(),Gly_0)
Initial(GSS(),GSS_0)
Initial(GSH(),GSH_0)
Initial(GSSG(),GSSG_0)
Initial(NADPH(),NADPH_0)
Initial(GR(),GR_0)
Initial(NADPplus(),NADPplus_0)
Initial(G6PD(),G6PD_0)
Initial(LOH(),LOH_0)



#Observables
Observable("LO_Obs",LO(ferrostatin=None))
Observable("GSH_Obs",GSH())
Observable("Cys_Obs",Cys())
Observable("Cystine_extra_Obs",Cystine_extra())
Observable("Cystine_intra_Obs",Cystine_intra())
Observable("NADPplus_Obs",NADPplus())
Observable("LOOH_Obs",LOOH())
Observable("LOH_Obs",LOH())
Observable("Glu_intra_Obs",Glu_intra())
Observable("Glu_extra_Obs",Glu_extra())
Observable("GPX_4_free",GPX4(rsl3=None))
Observable("Gly_Obs",Gly())
Observable("Glu_Cys_GCL_Product_Obs",Glu_Cys_GCL_Product())
Observable("Lipid_metab_Obs",Lipid_metab())
Observable("GSSG_Obs",GSSG())
Observable("NADPH_Obs",NADPH())




#Reactions
# Glu (intracellular) + Cystine (extracellular) + System Xc--> Glu (extracellular) + Cystine (intracellular) + System Xc-

Parameter("kcat_Glu_Cystine",1)
Parameter("km_Glu_Cystine",1e7)
Expression("keff_Glu_Cystine",kcat_Glu_Cystine / (km_Glu_Cystine+Glu_intra_Obs+Cystine_extra_Obs))
Rule("Glu_Cystine_Transport",Glu_intra() + Cystine_extra() + System_Xc(erastin=None) >> Glu_extra() + Cystine_intra() + System_Xc(erastin=None), keff_Glu_Cystine)

# System Xc- + Erastin   <--> System Xc- : Erastin
Parameter("kf_Xc_Erastin",1)
Parameter("kr_Xc_Erastin",1)
Rule("Xc_Erastin",System_Xc(erastin=None) + Erastin(sys_xc=None) | System_Xc(erastin=1) % Erastin(sys_xc=1),kf_Xc_Erastin,kr_Xc_Erastin)

# Cystine (intracellular)  -> Cys
Parameter("k_Cystine_intra",1)
Rule("Cystine_intra_to_cys",Cystine_intra() >> Cys(),k_Cystine_intra)

# Cys  + Trx/TrxR -> ? + Trx/TrxR

Parameter("kcat_Cys_TrxR",1)
Parameter("km_Cys_TrxR",100)
Expression("keff_Cys_TrxR",kcat_Cys_TrxR / (km_Cys_TrxR+Cys_Obs))
Rule("Cys_TrxR",Cys() + Trx_TrxR() >> Trx_TrxR(),keff_Cys_TrxR)

# Iron storage & trafficking heme metabolism + Iron <--> Iron storage & trafficking heme metabolism : Iron
Parameter("kf_Storage",1)
Parameter("kr_Storage",1)
Rule("Storage_hememetabolism",Iron_storage(iron=None) + Iron(b=None) | Iron_storage(iron=1) % Iron(b=1),kf_Storage,kr_Storage)

# Iron Chelators + Iron  <--> Iron Chelators : Iron
Parameter("kf_Chelators",1)
Parameter("kr_Chelators",1)
Rule("Chelators_Iron", Iron_chelators(iron=None) + Iron(b=None) | Iron(b=1) % Iron_chelators(iron=1),kf_Chelators,kr_Chelators)

# Ferrostatin-1 + LO. <--> Ferrostatin-1:LO.
Parameter("kf_Ferrostatin_LO",1)
Parameter("kr_Ferrostatin_LO",1)
Rule("Ferrostatin_LO",Ferrostatin(lo=None) + LO(ferrostatin=None) | Ferrostatin(lo=1) % LO(ferrostatin=1),kf_Ferrostatin_LO,kr_Ferrostatin_LO)

# Glu + Cys + GCL -> Glu_Cys_GCL_Product + GCL DONE
Parameter("kcat_Glu_Cys_GCL",1)
Parameter("km_Glu_Cys_GCL",100)
Expression("keff_Glu_Cys_GCL",kcat_Glu_Cys_GCL / (km_Glu_Cys_GCL+Cys_Obs+Glu_intra_Obs))
Rule("Glu_Cys_GCL",Glu_intra() + Cys()+ GCL() >> Glu_Cys_GCL_Product() + GCL(),keff_Glu_Cys_GCL)

# LOOH + Lipid_metab + Iron -> LO. + Iron
Parameter("kcat_LOOH_Lipid_metab",1)
Parameter("km_LOOH_Lipid_metab",100)
Expression("keff_LOOH_Lipid_metab",kcat_LOOH_Lipid_metab / (km_LOOH_Lipid_metab + LOOH_Obs + Lipid_metab_Obs))
Rule("LOOH_Lipid_metab",LOOH() + Lipid_metab() + Iron(b=None) >> LO(ferrostatin=None) + Iron(b=None),keff_LOOH_Lipid_metab)

# Glu_Cys_GCL_Product + Gly + GSS -> GSH + GSS DONE
Parameter("kcat_Gly_Glu_Cys_GCL_Product",1)
Parameter("km_Gly_Glu_Cys_GCL_Product",100)
Expression("keff_Gly_Glu_Cys_GCL_Product",kcat_Gly_Glu_Cys_GCL_Product / (km_Gly_Glu_Cys_GCL_Product+Gly_Obs+Glu_Cys_GCL_Product_Obs))
Rule("Gly_Glu_Cys_GCL_Product",Glu_Cys_GCL_Product() + Gly() + GSS() >> GSH() + GSS(),keff_Gly_Glu_Cys_GCL_Product)

# Example: LOOH + GSH + GPX4 -> LOH + GSSG + GPX4
Parameter("kcat_LOOH_GSH",1)
Parameter("km_LOOH_GSH",100)
Expression("keff_LOOH_GSH", kcat_LOOH_GSH / (km_LOOH_GSH + LOOH_Obs + GSH_Obs))
Rule("LOOH_LOH",LOOH() + GSH() + GPX4(rsl3=None) >> LOH() + GSSG() + GPX4(rsl3=None),keff_LOOH_GSH)

# RSL3 + GPX4 <--> RSL3:GPX4
Parameter("kf_RSL3_GPX4",1)
Parameter("kr_RSL3_GPX4",100)
Rule("RSL3_binds_GPX4",RSL3(gpx4=None) + GPX4(rsl3=None) | RSL3(gpx4=1) % GPX4(rsl3=1),kf_RSL3_GPX4,kr_RSL3_GPX4)

# PPARG -> PPARG + Lipid Metabolism
Parameter("k_PPARG_Lipid_Metab",1)
Rule("PPARG_Lipid_metab",PPARg() >> PPARg() + Lipid_metab(),k_PPARG_Lipid_Metab)

# Example 2: GSSG + NADPH + GR -> GSH + NADP + GR
Parameter("kcat_GSSG_NADPH",1)
Parameter("km_GSSG_NADPH",100)
Expression("keff_GSSG_NADPH", kcat_GSSG_NADPH / (km_GSSG_NADPH + GSSG_Obs + NADPH_Obs))
Rule("GSSG_GSH_NADPH",GSSG() + NADPH() + GR() >> GSH() + NADPplus() + GR(),keff_GSSG_NADPH)

# Example 3: NADP+ + G6PD -> NADPH + G6PD (only if one compound is required for reaction to occur)
Parameter("kcat_NADPplus_NADPH",1)
Parameter("km_NADPplus_NADPH",100)
Expression("keff_NADPplus_NADPH",kcat_NADPplus_NADPH / (km_NADPplus_NADPH + NADPplus_Obs))
Rule("NADPplus_NADPH",NADPplus() + G6PD() >> NADPH() + G6PD(),keff_NADPplus_NADPH)








Tspan=np.linspace(0,10,1001)
print(Tspan)

for km in [1e7,1e8]:
    sim=ScipyOdeSimulator(model,Tspan,verbose=True)
    result=sim.run(param_values={"km_Glu_Cystine": km})

    #obs2plot=["LOOH_Obs","LOH_Obs","LO_Obs"]
    obs2plot=[["LO_Obs"],
              ["LOOH_Obs","LOH_Obs"],
              ["GSH_Obs","GSSG_Obs"],
              ["NADPplus_Obs","NADPH_Obs"],
              ["Cystine_extra_Obs","Cystine_intra_Obs","Cys_Obs"],
              ["Glu_intra_Obs","Glu_extra_Obs"]]

    fig,axs=plt.subplots(nrows=3,ncols=2,sharex=True,figsize=(9.6,9.6))
    fig.suptitle("km_Glu_Cystine = %g" % km, fontsize=15)

    #use this code for plotting multiple km values
    #fig.suptitle("km_Glu_Cystine = %g,\n km_Glu_Cystine = %g,\n km_Glu_Cystine = %g" %
                 #(km_Glu_Cystine.value,km_Glu_Cystine.value,km_Glu_Cystine.value), fontsize=15)

    row=0
    col=0
    for obs_group in obs2plot:
        #plt.figure()
        for obs_name in obs_group:
            #model.observables:
            #plt.figure()
            maxconc=LOOH_0.value
            #plt.plot(Tspan,result.observables[obs_name]/maxconc,lw=2,label=obs_name)
            axs[row,col].plot(Tspan, result.observables[obs_name], lw=2, label=obs_name)
            if row == 2:
                axs[row,col].set_xlabel("Time")
            axs[row,col].set_ylabel("Concentration")
            axs[row,col].legend(loc="best")

        if (col+1) % 2 == 0:
            row += 1
            col = 0
        else:
            col += 1

            plt.tight_layout()
    plt.savefig("km_Glu_Cys_%g.pdf" % km,format="pdf")

#Use this code if planning to save files with multiple km changes
#plt.savefig("km_1_%g,km_2_%g,km_3_%g.pdf" % (km_1,km_2,km_3),format="pdf")
plt.show()


