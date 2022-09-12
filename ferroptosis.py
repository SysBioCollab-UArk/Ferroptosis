from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
#test commit
Model()

#Metabolites
Monomer("Glu_intra") #done
Monomer("Cystine_extra")
Monomer("Glu_extra") #done
Monomer("Cystine_intra")
Monomer("Cys") #done
Monomer("Glu") #done
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

Parameter("Glu_intra_0",100)
Parameter("Cystine_extra_0",100)
Parameter("System_Xc_0",10)
Parameter("Erastin_0",0)
Parameter("RSL3_0",0)
Parameter("GPX4_0",100)
Parameter("Cys_0",100)
Parameter("Trx_TrxR_0",100)
Parameter("Iron_storage_0",100)
Parameter("Iron_0",100)
Parameter("Iron_chelators_0",100)
Parameter("Ferrostatin_0",100)
Parameter("LO_0",100)
Parameter("Glu_0",100)
Parameter("GCL_0",100)
Parameter("Glu_Cys_GCL_Product_0",100)
Parameter("LOOH_0",100)
Parameter("Lipid_metab_0",100)
Parameter("Gly_0",100)
Parameter("GSS_0",100)
Parameter("GSH_0",100)
Parameter("GSSG_0",100)
Parameter("NADPH_0",100)
Parameter("GR_0",100)
Parameter("NADPplus_0",100)
Parameter("G6PD_0",100)
Parameter("LOH_0",100)
Parameter("PPARg_0",100)


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
Initial(Glu(),Glu_0)
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

Parameter("k_Glu_Cystine",1)
Parameter("kf_Xc_Erastin",1)
Parameter("kr_Xc_Erastin",1)

Observable("LO_Obs",LO(ferrostatin=None))
Observable("GSH_Obs",GSH())
Observable("Cys_Obs",Cys())
Observable("Cystine_extra_Obs",Cystine_extra())
Observable("NADPH_Obs",NADPH())
Observable("LOOH_Obs",LOOH())
Observable("LOH_Obs",LOH())

# Glu (intracellular) + Cystine (extracellular) + System Xc--> Glu (extracellular) + Cystine (intracellular) + System Xc-
Rule("Glu_Cystine_Transport",Glu_intra() + Cystine_extra() + System_Xc(erastin=None) >>
     Glu_extra() + Cystine_intra() + System_Xc(erastin=None), k_Glu_Cystine)

# System Xc- + Erastin   <--> System Xc- : Erastin
Rule("Xc_Erastin",System_Xc(erastin=None) + Erastin(sys_xc=None) | System_Xc(erastin=1) % Erastin(sys_xc=1),kf_Xc_Erastin,kr_Xc_Erastin)

# Cystine (intracellular)  -> Cys
Parameter("k_Cystine_intra",1)
Rule("Cystine_intra_to_cys",Cystine_intra() >> Cys(),k_Cystine_intra)

# Cys  + Trx/TrxR -> ? + Trx/TrxR
Parameter("k_Cys_TrxR",1)
Rule("Cys_TrxR",Cys() + Trx_TrxR() >> Trx_TrxR(),k_Cys_TrxR)

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

# Glu + Cys + GCL -> Glu_Cys_GCL_Product + GCL
Parameter("k_Glu_Cys_GCL",1)
Rule("Glu_Cys_GCL",Glu() + Cys()+ GCL() >> Glu_Cys_GCL_Product() + GCL(),k_Glu_Cys_GCL)

# LOOH + Lipid_metab + Iron -> LO. + Iron
Parameter("k_LOOH_Lipid_metab_Iron",1)
Rule("LOOH_Lipid_metab_Iron",LOOH() + Lipid_metab() + Iron(b=None) >> LO(ferrostatin=None) + Iron(b=None),k_LOOH_Lipid_metab_Iron)

# Glu_Cys_GCL_Product + Gly + GSS -> GSH + GSS
Parameter("k_Gly_Glu_Cys_GCL_Product",1)
Rule("Gly_Glu_Cys_GCL_Product",Glu_Cys_GCL_Product() + Gly() + GSS() >> GSH() + GSS(),k_Gly_Glu_Cys_GCL_Product)

# Example: LOOH + GSH + GPX4 -> LOH + GSSG + GPX4
Parameter("k_LOOH_LOH",1)
Rule("LOOH_LOH",LOOH() + GSH() + GPX4(rsl3=None) >> LOH() + GSSG() + GPX4(rsl3=None),k_LOOH_LOH)

# RSL3 + GPX4 <--> RSL3:GPX4
Parameter("kf_RSL3_GPX4",1)
Parameter("kr_RSL3_GPX4",1)
Rule("RSL3_GPX4",RSL3(gpx4=None) + GPX4(rsl3=None) | RSL3(gpx4=1) % GPX4(rsl3=1),kf_RSL3_GPX4,kr_RSL3_GPX4)

# PPARG -> PPARG + Lipid Metabolism
Parameter("k_PPARG_Lipid_Metab",1)
Rule("PPARG_Lipid_metab",PPARg() >> PPARg() + Lipid_metab(),k_PPARG_Lipid_Metab)

# Example 2: GSSG + NADPH + GR -> GSH + NADP+ + GR
Parameter("k_GSSG_GSH",1)
Rule("GSSG_GSH",GSSG() + NADPH() + GR() >> GSH() + NADPplus() + GR(),k_GSSG_GSH)

# Example 3: NADP+ + G6PD -> NADPH + G6PD (only if one compound is required for reaction to occur)
Parameter("k_NADPplus_NADPH",1)
Rule("NADPplus_NADPH",NADPplus() + G6PD() >> NADPH() + G6PD(),k_NADPplus_NADPH)








Tspan=np.linspace(0,0.1,101)
print(Tspan)
sim=ScipyOdeSimulator(model,Tspan,verbose=True)
result=sim.run()

for obs in model.observables:
    plt.plot(Tspan,result.observables[obs.name],lw=2,label=obs.name)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend(loc="best")

plt.tight_layout()
plt.show()
