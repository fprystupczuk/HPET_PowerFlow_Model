#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import math
from collections import namedtuple

# Packages
from MainPackages.RunDSS import RunDSS
from MainPackages.PlotDSS import PlotDSS
from MainPackages.CreateDSS import CreateDSS
from MainPackages.Extract_OpenDSS_data import Extract_data_OpenDSS


#models = ['LFT']   # 'LFT','PET','HPET'
#models = ['PET']   # 'LFT','PET','HPET'
#models = ['HPET']   # 'LFT','PET','HPET'
#models = ['LFT','HPET']   # 'LFT','PET','HPET'
#models = ['PET','HPET']   # 'LFT','PET','HPET'
models = ['LFT','PET','HPET']   # 'LFT','PET','HPET'
#WhatToRun = 'Simulations' # DSS or Plot or DSSandPlot
#WhatToRun = 'Plots' # DSS or Plot or DSSandPlot
WhatToRun = 'SimulationsAndPlots' # DSS or Plot or DSSandPlot



MV_loads = 'n'  # Whether to connect MV loads or not
# Simulation time. Provided profiles have a 5 min resolution.
# First time step >= 5, last time step <= 1440
timeStruct = namedtuple("timeStruct", "start end step")
Time = timeStruct(start=5, end=1440, step=5)

# Create time set
timeArray = np.array(range(Time.start,Time.end+1,Time.step)) # Time set for simulations

# Study feeders - Including all feeders will simulate the complete LV Network - IMPORTANT: this must match the model in 'OPF_model' folder
Network = 'network_12'
Feeders = ['Feeder_1','Feeder_2','Feeder_3'] # ['Feeder1','Feeder2','Feeder3'] - can use all three or 1 or 2

# Upstream netwok line names for selected variables
SelectedLinesUS = np.array(['mv_line'])

# Downstream netwok line names for selected variables
#SelectedLinesDS = np.array(['feeder1_line1','feeder1_line237',
#                            'feeder2_line1','feeder2_line155',
#                            'feeder3_line1','feeder3_line149'])
SelectedLinesDS = ''

# Power ratings of the power electronics devices
PowerRatingPET = 0.5
FractionPowerHPET = 0.1

# PV characteristics
PV_penetration = 0.0 # [0.0 - 1.0] number_of_customers/customers_with_PV

# Voltage control mode
Vcontrol = np.array(['Vfixed',[1.0,1.0]]) # Mode: Vfixed | Vband | V_CVR , Value: [max, min] pu
#Vcontrol = np.array(['Vband',[0.95,1.05]]) # Mode: Vfixed | Vband | V_CVR , Value: [max, min] pu
#Vcontrol = np.array(['V_CVR', 1.0]) # Mode: Vfixed | V_CVR , Value: Vref pu

# Name of the electric variable to be collected per line or bus by Extract_DSSmonitors_data.py
ElectricVars = np.array(['Vmag_send','Vang_send','Vmag_rec','Vang_rec','Imag','P_rec','Q_rec','P_send','Q_send'])

# Nominal voltages
#VoltageBases = [10,0.4,0.04] # in kV line to line

# Folders
Main_path = os.path.dirname(os.path.realpath(__file__))
NetworkModels_path = Main_path + '/NetworkData/NetworkModels'
Profiles_path = Main_path + '/NetworkData/Profiles'
Main_Results_path = Main_path + '/Results'


for model in models:
    # Create OpenDSS object, import network and solve
    [DSSObj,DSSText,DSSCircuit,DSSElement,DSSSolution,DSSMonitors,CableData,TransformerRating] = CreateDSS(pd,math,model,Main_path,Network,Feeders,MV_loads,PowerRatingPET,FractionPowerHPET)
    # Extract Network data
    [Lines_set,Lines_upstream,Lines_downstream,Lines_upstream_data,Lines_downstream_data,Buses_set,Buses_upstream,Buses_downstream,Bus_Vnom,Nodes_set,Loads_LV,Load_phase,Load_bus,Load_Vnom] = Extract_data_OpenDSS(math,np,pd,DSSObj,model)
    # Run simulation if required
    if WhatToRun != 'Plots':
        RunDSS(os,np,pd,math,Main_path,model,Network,Feeders,Vcontrol,PV_penetration,ElectricVars,Profiles_path,Main_Results_path,timeArray,DSSObj,Lines_set,Loads_LV,Load_bus,Load_phase,Bus_Vnom,MV_loads,TransformerRating,PowerRatingPET,FractionPowerHPET)

# Plot saved results
if WhatToRun != 'Simulations':
    PlotDSS(os,np,Main_Results_path,Time,timeArray,models,Network,Feeders,Vcontrol,PV_penetration,Lines_set,SelectedLinesUS,SelectedLinesDS,MV_loads,Bus_Vnom,Lines_upstream_data,Lines_downstream_data,PowerRatingPET,FractionPowerHPET)
