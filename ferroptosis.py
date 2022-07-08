from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer("Glu_intra")
Monomer("Cystine_extra")
Monomer("System_Xc",["b"])
Monomer("Glu_extra")
Monomer("Cystine_intra")
Monomer("Erastin",["b"])
Monomer("RSL3",["b"])
Monomer("GPX4",["b"])
Monomer("Cys")
Monomer("Trx_TrxR")
Monomer("Iron_storage",["b"])
Monomer("Iron",["b"])
Monomer("Iron_chelators",["b"])
Monomer("Ferrostatin",["b"])
Monomer("LO",["b"])
Monomer("Glu")
Monomer("GCL")
Monomer("Intermediate")
Monomer("LOOH")
Monomer("Lipid")
Monomer("Gly")
Monomer("GSS")
Monomer("GSH")
Monomer("GSSG")
Monomer("NADPH")
Monomer("GR")
Monomer("NADPplus")
Monomer("G6PD")
Monomer("LOH")

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
Parameter("Intermediate_0",100)
Parameter("LOOH_0",100)
Parameter("Lipid_0",100)
Parameter("Gly_0",100)
Parameter("GSS_0",100)
Parameter("GSH_0",100)
Parameter("GSSG_0",100)
Parameter("NADPH_0",100)
Parameter("GR_0",100)
Parameter("NADPplus_0",100)
Parameter("G6PD_0",100)
Parameter("LOH_0",100)


Initial(Glu_intra(),Glu_intra_0)
Initial(Cystine_extra(),Cystine_extra_0)
Initial(System_Xc(b=None),System_Xc_0)
Initial(Erastin(b=None),Erastin_0)
Initial(RSL3(b=None),RSL3_0)
Initial(GPX4(b=None),GPX4_0)
Initial(Cys(),Cys_0)
Initial(Trx_TrxR(),Trx_TrxR_0)
Initial(Iron_storage(b=None),Iron_storage_0)
Initial(Iron(b=None),Iron_0)
Initial(Iron_chelators(b=None),Iron_chelators_0)
Initial(Ferrostatin(b=None),Ferrostatin_0)
Initial(LO(b=None),LO_0)
Initial(Glu(),Glu_0)
Initial(GCL(),GCL_0)
Initial(Intermediate(),Intermediate_0)
Initial(LOOH(),LOOH_0)
Initial(Lipid(),Lipid_0)
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

Observable("LO_Obs",LO(b=None))
Observable("GSH_Obs",GSH())
Observable("Cys_Obs",Cys())
Observable("Cystine_extra_Obs",Cystine_extra())
Observable("NADPH_Obs",NADPH())
Observable("LOOH_Obs",LOOH())
Observable("LOH_Obs",LOH())

# Glu (intracellular) + Cystine (extracellular) + System Xc--> Glu (extracellular) + Cystine (intracellular) + System Xc-
Rule("Glu_Cystine_Transport",Glu_intra() + Cystine_extra() + System_Xc(b=None) >>
     Glu_extra() + Cystine_intra() + System_Xc(b=None), k_Glu_Cystine)

# System Xc- + Erastin   <--> System Xc- : Erastin
Rule("Xc_Erastin",System_Xc(b=None) + Erastin(b=None) | System_Xc(b=1) % Erastin(b=1),kf_Xc_Erastin,kr_Xc_Erastin)

# Cystine (intracellular)  -> Cys
Parameter("k_Cystine_intra",1)
Rule("Cystine_intra_to_cys",Cystine_intra() >> Cys(),k_Cystine_intra)

# Cys  + Trx/TrxR -> ? + Trx/TrxR
Parameter("k_Cys_TrxR",1)
Rule("Cys_TrxR",Cys() + Trx_TrxR() >> Trx_TrxR(),k_Cys_TrxR)

# Iron storage & trafficking heme metabolism + Iron <--> Iron storage & trafficking heme metabolism : Iron
Parameter("kf_Storage",1)
Parameter("kr_Storage",1)
Rule("Storage_hememetabolism",Iron_storage(b=None) + Iron(b=None) | Iron_storage(b=1) % Iron(b=1),kf_Storage,kr_Storage)

# Iron Chelators + Iron  <--> Iron Chelators : Iron
Parameter("kf_Chelators",1)
Parameter("kr_Chelators",1)
Rule("Chelators_Iron", Iron_chelators(b=None) + Iron(b=None) | Iron(b=1) % Iron_chelators(b=1),kf_Chelators,kr_Chelators)

# Ferrostatin-1 + LO. <--> Ferrostatin-1:LO.
Parameter("kf_Ferrostatin_LO",1)
Parameter("kr_Ferrostatin_LO",1)
Rule("Ferrostatin_LO",Ferrostatin(b=None) + LO(b=None) | Ferrostatin(b=1) % LO(b=1),kf_Ferrostatin_LO,kr_Ferrostatin_LO)

# Glu + Cys + GCL -> Intermediate + GCL
Parameter("k_Glu_Cys_GCL",1)
Rule("Glu_Cys_GCL",Glu() + Cys()+ GCL() >> Intermediate() + GCL(),k_Glu_Cys_GCL)

# LOOH + Lipid + Iron -> LO. + Iron
Parameter("k_LOOH_Lipid_Iron",1)
Rule("LOOH_Lipid_Iron",LOOH() + Lipid() + Iron(b=None) >> LO(b=None) + Iron(b=None),k_LOOH_Lipid_Iron)

# Intermediate + Gly + GSS -> GSH + GSS
Parameter("k_Gly_Intermediate",1)
Rule("Gly_Intermediate",Intermediate() + Gly() + GSS() >> GSH() + GSS(),k_Gly_Intermediate)

# Example: LOOH + GSH + GPX4 -> LOH + GSSG +GPX4
Parameter("k_LOOH_LOH",1)
Rule("LOOH_LOH",LOOH() + GSH() + GPX4(b=None) >> LOH() + GSSG() + GPX4(b=None),k_LOOH_LOH)

# RSL3 + GPX4 <--> RSL3:GPX4
Parameter("kf_RSL3_GPX4",1)
Parameter("kr_RSL3_GPX4",1)
Rule("RSL3_GPX4",RSL3(b=None) + GPX4(b=None) | RSL3(b=1) % GPX4(b=1),kf_RSL3_GPX4,kr_RSL3_GPX4)

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
