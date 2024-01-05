from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

#Metabolites
Monomer("Glu_intra") #done
Monomer("Cystine_extra")
Monomer("Glu_extra") #done
Monomer("Cystine_intra")
Monomer("Cys") #done
Monomer("Glu_Cys_GCL_Product")
Monomer("Gly")
Monomer("GSH")

#Enzymes
Monomer("System_Xc",["erastin"])
Monomer("Trx_TrxR")
Monomer("GCL")
Monomer("GSS")

#Observables
Observable("Cys_Obs",Cys())
Observable("Cystine_extra_Obs",Cystine_extra())
Observable("Cystine_intra_Obs",Cystine_intra())
Observable("Glu_intra_Obs",Glu_intra())
Observable("Glu_extra_Obs",Glu_extra())
Observable("Gly_Obs",Gly())
Observable("GSH_Obs",GSH())
Observable("Glu_Cys_GCL_Product_Obs",Glu_Cys_GCL_Product())

#Parameters
Parameter("System_Xc_0",100)
Parameter("Trx_TrxR_0",100)
Parameter("GCL_0",100)
Parameter("GSS_0",100)
Parameter("Gly_0",0)

#Initials
Initial(System_Xc(erastin=None),System_Xc_0)
Initial(Trx_TrxR(),Trx_TrxR_0)
Initial(GCL(),GCL_0)
Initial(GSS(),GSS_0)
Initial(Gly(),Gly_0)

#Rules
Parameter("k_Glu_intra",2e3) #1.95e5 originally
Rule("Glu_Intra_Synthesis",None >> Glu_intra(),k_Glu_intra)

Parameter("k_Glu_extra",0.1) #5 originally
Rule("Glu_Extra_Degradation",Glu_extra() >> None, k_Glu_extra)

Parameter("k_Cystine_extra",1e3) #1e5 originally
Rule("Cystine_Extra_Synthesis",None >> Cystine_extra(),k_Cystine_extra)

Parameter("kcat_Glu_Cystine",1e-7)#1e-3
#Parameter("km_Glu_Cystine",1e7)
#Expression("keff_Glu_Cystine",kcat_Glu_Cystine / (km_Glu_Cystine+Glu_intra_Obs+Cystine_extra_Obs))
Rule("Glu_Cystine_Transport",Glu_intra() + Cystine_extra() + System_Xc(erastin=None) >>
     Glu_extra() + Cystine_intra() + System_Xc(erastin=None), kcat_Glu_Cystine)

# Cystine (intracellular)  -> Cys #CHANGING MAKES THINGS GO WACK
Parameter("k_Cystine_intra",0.1) #originally 5
Rule("Cystine_intra_to_cys",Cystine_intra() >> Cys(),k_Cystine_intra)

# Cys  + Trx/TrxR -> ? + Trx/TrxR #Note: Giving value here is important to stabilize
# Cys_intra, Glu_intra and Cys_extra concentrations
Parameter("kcat_Cys_TrxR",0) #1e-4
#Parameter("km_Cys_TrxR",100)
#Expression("keff_Cys_TrxR",kcat_Cys_TrxR / (km_Cys_TrxR+Cys_Obs))
Rule("Cys_TrxR",Cys() + Trx_TrxR() >> Trx_TrxR(),kcat_Cys_TrxR)

# Glu + Cys + GCL -> Glu_Cys_GCL_Product + GCL DONE
Parameter("kcat_Glu_Cys_GCL",2e-7) #1e-4 originally
#Parameter("km_Glu_Cys_GCL",1e7) #100
#Expression("keff_Glu_Cys_GCL",kcat_Glu_Cys_GCL / (km_Glu_Cys_GCL+Cys_Obs+Glu_intra_Obs))
Rule("Glu_Cys_GCL",Glu_intra() + Cys() + GCL() >> Glu_Cys_GCL_Product() + GCL(),kcat_Glu_Cys_GCL)

# Glu_Cys_GCL_Product + Gly + GSS -> GSH + GSS DONE
Parameter("kcat_Gly_Glu_Cys_GCL_Product",0.001) #1 0.01
Parameter("km_Gly_Glu_Cys_GCL_Product",1e3)#100
Expression("keff_Gly_Glu_Cys_GCL_Product",kcat_Gly_Glu_Cys_GCL_Product / (km_Gly_Glu_Cys_GCL_Product+Gly_Obs+Glu_Cys_GCL_Product_Obs))
Rule("Gly_Glu_Cys_GCL_Product",Glu_Cys_GCL_Product() + Gly() + GSS() >> GSH() + GSS(),keff_Gly_Glu_Cys_GCL_Product)

Parameter("k_Gly",75)#100
Rule("Gly_Synthesis",None >> Gly(),k_Gly)

Parameter("k_GSH_deg",0.05)#100
Rule("GSH_Degradation",GSH() >> None, k_GSH_deg)

Tspan=np.linspace(0,200,1001)
print(Tspan)

for km in [100]:
    sim=ScipyOdeSimulator(model,Tspan,verbose=True)
    result=sim.run()#()(param_values={"km_LOOH_GSH": km})

    #obs2plot=["LOOH_Obs","LOH_Obs","LO_Obs"]
    obs2plot=[["Cystine_extra_Obs","Cystine_intra_Obs","Cys_Obs"],
              ["Glu_intra_Obs","Glu_extra_Obs"],
              ["Gly_Obs", "GSH_Obs"]]

    fig,axs=plt.subplots(nrows=len(obs2plot),ncols=1,sharex=True,figsize=(9.6,9.6))
    #print(len(axs))
    #quit()
    fig.suptitle("km_LOOH_GSH = %g" % km, fontsize=15)



    row=0
    col=0
    for obs_group in obs2plot:
        for obs_name in obs_group:

            axs[row].plot(Tspan, result.observables[obs_name], lw=2, label=obs_name)
            if row == 1:
                axs[row].set_xlabel("Time",fontsize=16)
            axs[row].set_ylabel("Concentration",fontsize=16)
            axs[row].legend(loc="best",fontsize=16)
            axs[row].xaxis.set_tick_params(labelsize=16)
            axs[row].yaxis.set_tick_params(labelsize=16)

        if (col+1) % 1 == 0:
            row += 1
            col = 0
        else:
            col += 1


    plt.tight_layout()
    #plt.savefig("km_Glu_Cys_%g.pdf" % km,format="pdf")

plt.show()