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

#Enzymes
Monomer("System_Xc",["erastin"])
Monomer("Trx_TrxR")
Monomer("GCL")

#Parameters
Parameter("System_Xc_0",100)
Parameter("Trx_TrxR_0",100)
Parameter("GCL_0",100)

#Initials
Initial(System_Xc(erastin=None),System_Xc_0)
Initial(Trx_TrxR(),Trx_TrxR_0)
Initial(GCL(),GCL_0)

#Rules
Parameter("k_Glu_intra",1e5)
Rule("Glu_Intra_Synthesis",None >> Glu_intra(),k_Glu_intra)

Parameter("k_Cystine_extra",1e5)
Rule("Cystine_Extra_Synthesis",None >> Cystine_extra(),k_Cystine_extra)

Parameter("kcat_Glu_Cystine",1e-3)
#Parameter("km_Glu_Cystine",1e7)
#Expression("keff_Glu_Cystine",kcat_Glu_Cystine / (km_Glu_Cystine+Glu_intra_Obs+Cystine_extra_Obs))
Rule("Glu_Cystine_Transport",Glu_intra() + Cystine_extra() + System_Xc(erastin=None) >>
     Glu_extra() + Cystine_intra() + System_Xc(erastin=None), kcat_Glu_Cystine)

# Cystine (intracellular)  -> Cys
Parameter("k_Cystine_intra",1)
Rule("Cystine_intra_to_cys",Cystine_intra() >> Cys(),k_Cystine_intra)

# Cys  + Trx/TrxR -> ? + Trx/TrxR
Parameter("kcat_Cys_TrxR",1)
#Parameter("km_Cys_TrxR",100)
#Expression("keff_Cys_TrxR",kcat_Cys_TrxR / (km_Cys_TrxR+Cys_Obs))
Rule("Cys_TrxR",Cys() + Trx_TrxR() >> Trx_TrxR(),kcat_Cys_TrxR)

# Glu + Cys + GCL -> Glu_Cys_GCL_Product + GCL DONE
Parameter("kcat_Glu_Cys_GCL",0)
#Parameter("km_Glu_Cys_GCL",1e7) #100
#Expression("keff_Glu_Cys_GCL",kcat_Glu_Cys_GCL / (km_Glu_Cys_GCL+Cys_Obs+Glu_intra_Obs))
Rule("Glu_Cys_GCL",Glu_intra() + Cys()+ GCL() >> GCL(),kcat_Glu_Cys_GCL)

#Observables
Observable("Cys_Obs",Cys())
Observable("Cystine_extra_Obs",Cystine_extra())
Observable("Cystine_intra_Obs",Cystine_intra())
Observable("Glu_intra_Obs",Glu_intra())
Observable("Glu_extra_Obs",Glu_extra())

Tspan=np.linspace(0,10,1001)
print(Tspan)

for km in [100]:
    sim=ScipyOdeSimulator(model,Tspan,verbose=True)
    result=sim.run()#()(param_values={"km_LOOH_GSH": km})

    #obs2plot=["LOOH_Obs","LOH_Obs","LO_Obs"]
    obs2plot=[["Cystine_extra_Obs","Cystine_intra_Obs","Cys_Obs"],
              ["Glu_intra_Obs"]]

    fig,axs=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(9.6,9.6))
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