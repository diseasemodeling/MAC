
import pandas as pd
import numpy as np
from Multiple_disease.Multiple_disease_module.agents import Agents
# import cProfile 
# import pstats
# from numba import njit,prange
import pdb
# import random
import os
# import pickle



simulation_step = 41
time_unit = 2
initial_year = 2006
dryrun_step = 70

 # increased risk =0 no risk, 1 all increased risk, 2 HPV, 3 progression
def runExperiment(is_screening, is_vaccine, increased_risk, random_seed, write, plot):
        np.random.seed(random_seed)
        kwargs = {"mode":"python",
        "increased_risk": increased_risk,
        "initial_year": initial_year,
        "dryrun_step":dryrun_step,
        "simulation_step": simulation_step,
        "model_screening": is_screening,
        "model_vaccine":is_vaccine,
        "time_unit": time_unit,
        "plot": plot
        }
        d2 = Agents(**kwargs)
        d2.dryrun()
        d2.disease_state_assignment()
        d2.run(time = simulation_step * d2.time_unit)
        d2.output.output_results(random_seed, time_unit,increased_risk, is_screening, is_vaccine, write)



runExperiment(is_screening = True,
              is_vaccine = True,
              plot = False,
              increased_risk = 0,
              random_seed = 1024,
              write = True)
os.system('osascript -e \'beep\'')

