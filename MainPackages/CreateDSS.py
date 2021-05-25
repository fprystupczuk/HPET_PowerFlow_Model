#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import platform

def CreateDSS(pd,math,model,Main_path,Network,Feeders,MV_loads,PowerRatingPET,FractionPowerHPET):

    # Create OpenDSS object
    if platform.system() == 'Linux':
        import dss # You need to have dss_python installed (https://github.com/dss-extensions/dss_python): on all major platforms, you can install it directly from pip: pip install dss_python. Or, if you’re using the Anaconda distribution, you can use: conda install -c pmeira dss_python
        DSSObj = dss.DSS
    
    if platform.system() == 'Windows':
        import win32com.client
        from win32com.client import makepy
        import sys
        sys.argv = ["makepy", "OpenDSSEngine.DSS"]
        makepy.main()
        DSSObj = win32com.client.Dispatch("OpenDSSEngine.DSS")
    
    DSSText = DSSObj.Text # to excecute OpenDSS text commands
    DSSCircuit = DSSObj.ActiveCircuit # Use it to access the elements in the circuit (e.g., capacitors, buses, etc.)
    DSSElement = DSSCircuit.ActiveCktElement
    DSSSolution = DSSCircuit.Solution
    DSSMonitors = DSSCircuit.Monitors
    
    # This is the datapath were the DSS file is. Results will also be saved here.
    DSSText.Command = 'set datapath=' + Main_path
    DSSText.Command = 'Clear' # clear any existing circuit in the engine
    DSSText.Command = 'Set DefaultBaseFrequency=50'
    
    # Vsource
    DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/OpenDSS_NewCircuit.txt'
    
    # Line codes
    DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/Cable_data/OpenDSS_line_types.txt'
    
    # Cable type data
    CableData = pd.read_excel(Main_path + '/NetworkData/NetworkModels/Cable_data/Cable_data.xlsx',index_col=0 )
    CableData.index = CableData.index.str.lower()
    
    # Lines and loads
    for i_feeder in Feeders:
        # Compile OpenDSS codes for lines, loads and monitors
        DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/' + Network + '/' + i_feeder + '/OpenDSS_Lines.txt'
        DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/' + Network + '/' + i_feeder + '/OpenDSS_Loads.txt'
        if MV_loads == 'y': DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/OpenDSS_MVloads.txt'

    # Transformer model
    if model == 'LFT':
        DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/Transformers/OpenDSS_Transformer.txt'
        DSSCircuit.Transformers.Wdg=1
        TransformerRating = DSSCircuit.Transformers.kva
        DSSText.Command ='set VoltageBases=[' + str(10.0) + ',' + str(0.4) + ']'
    elif model == 'PET':
        DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/Transformers/OpenDSS_PETmodel.txt'
        TransformerRating = PowerRatingPET*(800/3) # Per-phase power rating
        DSSText.Command ='set VoltageBases=[' + str(10.0) + ',' + str(0.4) + ']'
    elif model == 'HPET':
        DSSText.Command = 'Redirect ' + Main_path + '/NetworkData/NetworkModels/Transformers/OpenDSS_HPETmodel.txt'        
        DSSCircuit.Transformers.Name = 'Trafo1'
        DSSCircuit.Transformers.Wdg = 1
        TransformerRating = DSSCircuit.Transformers.kva
        Smax_auxiliary = (TransformerRating*FractionPowerHPET)/(FractionPowerHPET+1)
        
        # Redefine power ratings and nominal voltages for transformers, STvsources and STloads
        for phase in range(3):
            DSSText.Command = 'Edit Transformer.Trafo' + str(phase+1) + ' Wdg=3 kV=' + str(0.4*FractionPowerHPET) + ' kVA=' + str(Smax_auxiliary)
            DSSText.Command = 'Edit Vsource.STvsource' + str(phase+1) + ' BasekV=' + str(0.4*FractionPowerHPET/math.sqrt(3))
            DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kV=' + str(0.4*FractionPowerHPET)
            
        # Redefine voltage bases
        DSSText.Command = 'set voltageBases=[' + str(10.0) + ',' + str(0.4) + ',' + str(0.4*FractionPowerHPET) + ']'
    
    DSSText.Command = 'CalcVoltageBase'
    DSSText.Command = 'set controlmode=static'
    DSSText.Command = 'set mode=snapshot'
    
    DSSSolution.Solve()
    if DSSSolution.Converged:
        print('Solution has converged')
    else:
        raise ValueError('Solution has not converged')
    
    return [DSSObj,DSSText,DSSCircuit,DSSElement,DSSSolution,DSSMonitors,CableData,TransformerRating]
