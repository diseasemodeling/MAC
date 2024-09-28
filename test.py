import pynetlogo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from SALib.sample import sobol as sobolsample
from SALib.analyze import sobol
import ipyparallel as ipp
import os

netlogo = pynetlogo.NetLogoLink(
    netlogo_home ="/Applications/NetLogo 6.2.2",
    jvm_path="/Applications/NetLogo 6.2.2/JRE/Contents/Home/jre/lib/server/libjvm.dylib"
)

netlogo.load_model("/Users/eleanor/Documents/GitHub/PATH4-public-main/BaseModel.nlogo")
# netlogo.command("runExperiment")

