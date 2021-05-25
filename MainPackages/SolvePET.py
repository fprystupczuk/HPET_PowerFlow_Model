#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy import interpolate

def SolvePET(Main_path,np,Vcontrol,DSSObj,Bus_Vnom,Loads_LV,Load_phase,TransformerRating):
    
    DSSText = DSSObj.Text # to excecute OpenDSS text commands
    DSSCircuit = DSSObj.ActiveCircuit # Use it to access the elements in the circuit (e.g., capacitors, buses, etc.)
    DSSSolution = DSSCircuit.Solution
    DSSMonitors = DSSCircuit.Monitors

    Vloads_LV = np.full([Loads_LV.size,3],np.nan) # LV load voltages split per phase
    Vloads_min = np.zeros(3)
    DeltaV = np.zeros(3) # voltage, angle for the STvsource elements
    V_LV_nom = Bus_Vnom.loc['secondary','Vnom_pn'] # Min statutory LV voltage limit    
    Mv = 0.002 # Security margin for minimum voltage regulation
    V_LV_min = (0.9 + Mv)*V_LV_nom # Min statutory LV voltage limit plus security margin Mv
    V_LV_max = (1.1 - Mv)*V_LV_nom # Max statutory LV voltage limit minus security margin Mv
    Ps = np.zeros(3)
    Qs = np.zeros(3)

    # Efficiency curves and interpolation function
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_MatlabLossModel_UnityPF.txt","rb"),delimiter=",")
    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_KACOblueplanet125TL3.txt","rb"),delimiter=",")
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_XantrexGT3.8Inverter.txt","rb"),delimiter=",")
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_Qin2010.txt","rb"),delimiter=",")
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_Adabi2018.txt","rb"),delimiter=",")
    eff_points[:,1] = (0.975/np.nanmax(eff_points[:,1])) * eff_points[:,1] # offset efficiency to a different maximum
    eff_interpFunction = interpolate.interp1d(eff_points[:,0],eff_points[:,1], kind='linear', fill_value="extrapolate")
    efficiency = np.zeros(3)
    
    
    if Vcontrol[0] == 'V_CVR': # otherwise the STvsources pu is always 1.0
        # In the first place: calculate converter's output voltage deppending on V_LV
        
        # Rseset STvsources
        for phase in range(3):
            DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
            DSSCircuit.Vsources.pu = 1.0

        ## Solve to get the current voltages
        DSSText.Command ='Reset Monitors'
        DSSSolution.Solve()
        if not(DSSSolution.Converged):
            raise ValueError('Solution did not Converge')
        DSSMonitors.SampleAll() 
        DSSMonitors.SaveAll()
        
        Pagregated = 0
        ZIPagregated = np.zeros(3)
    
        # Get all LV loads voltage, per phase
        for i_load in range(Loads_LV.size): # Loop thought all loads in the system
            DSSMonitors.Name = 'Monitor_' + Loads_LV[i_load] + '_vi'
            Vloads_LV[i_load,Load_phase.at[Loads_LV[i_load],'phase']-1] = DSSMonitors.Channel(1)[0]
            DSSCircuit.Loads.Name = Loads_LV[i_load]
            Pagregated += DSSCircuit.Loads.kW # Adds al the load power for the next calculation
            
        # Obtains voltage sensitivity VS by computing agregated load: equations (17) and (23) from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7725539
        for i_load in range(Loads_LV.size):
            DSSCircuit.Loads.Name = Loads_LV[i_load]
            ZIP = DSSCircuit.Loads.ZIPV[[0,1,2]] # Obtains only the first 3 components which correspond to the active power coefficients
            ZIPagregated += (DSSCircuit.Loads.kW/Pagregated)*ZIP
            VS = (2*ZIPagregated[0] + ZIPagregated[1] + ZIPagregated[2])/sum(ZIPagregated)

        # Deppending on the voltage sensitivity, regulates the output voltage for a min or max value depending on load voltages actual values
        if VS > 0: # Obtains the minimum voltage for each phase and calculates the available headroom
            Vloads_min = np.nanmin(Vloads_LV,0)
            DeltaV = (V_LV_min - Vloads_min)
        elif VS < 0: # Obtains the maximum voltage for each phase and calculates the available headroom
            Vloads_max = np.nanmax(Vloads_LV,0)
            DeltaV = (V_LV_max - Vloads_max)
        else: # If no sensitivity, then keeps the voltage at 1pu
            DeltaV = 0

    # Update STvsources pu and agnle parameters (this is needed here to obtain the power for this set-point)
    DSSCircuit.Vsources.Name = 'STvsource1'
    DSSCircuit.Vsources.pu = 1.0 + DeltaV[0]/V_LV_nom
    DSSCircuit.Vsources.Name = 'STvsource2'
    DSSCircuit.Vsources.pu = 1.0 + DeltaV[1]/V_LV_nom
    DSSCircuit.Vsources.Name = 'STvsource3'
    DSSCircuit.Vsources.pu = 1.0 + DeltaV[2]/V_LV_nom


    ## Solve to get the actual Pload and Qload for the operating point previously set
    DSSText.Command ='Reset Monitors'
    DSSSolution.Solve()
    if not(DSSSolution.Converged):
        raise ValueError('Solution did not Converge')
    DSSMonitors.SampleAll()
    DSSMonitors.SaveAll()

    # Read secondary P,Q per phase. Ps and Qs are measured at Vsource's output, thus they ar negative
    DSSMonitors.Name = 'STvsource1-PQ'
    Ps[0] = -DSSMonitors.Channel(1)[0]
    Qs[0] = -DSSMonitors.Channel(2)[0]
    DSSMonitors.Name = 'STvsource2-PQ'
    Ps[1] = -DSSMonitors.Channel(1)[0]
    Qs[1] = -DSSMonitors.Channel(2)[0]
    DSSMonitors.Name = 'STvsource3-PQ'
    Ps[2] = -DSSMonitors.Channel(1)[0]
    Qs[2] = -DSSMonitors.Channel(2)[0]
    
    Pp = np.zeros(3)
    Qp = np.zeros(3)
    
    # Normalization of S, per phase
    Snom = TransformerRating # per-phase PET nominal apparent power in kVA
    Ss = np.sqrt(Ps**2 + Qs**2)
    Spu = Ss / Snom
    PF = abs(Ps/Ss)
        
    # Calculate efficiency and active power in the primary side, per phase
    for phase in range(3):
        efficiency[phase] = eff_interpFunction(Spu[phase]) # interpolate actual efficiency
        if Ps[phase] >=0: # forward power flow
            Pp[phase] = Ss[phase] / efficiency[phase]
        else: # reverse power flow
            Pp[phase] = Ss[phase] * efficiency[phase]
    
    # This could be changed in the future, but here's explicitly set to zero for total Q compensation
    Qp = np.zeros(3)

    # Update STloads values (in this case the PET balance the 3 phases)
    DSSCircuit.Loads.Name = 'STload'
    DSSCircuit.Loads.kW = sum(Pp)
    DSSCircuit.Loads.kvar = sum(Qp)

    
    # SOLVE: final result for the current time step
    DSSText.Command ='Reset Monitors'
    DSSSolution.Solve()
    if not(DSSSolution.Converged):
        raise ValueError('Solution did not Converge')        
    DSSMonitors.SampleAll() 
    DSSMonitors.SaveAll()

    
    # Calculate total values and exit
    loss_total = sum(abs(Pp - Ps)) # Losses in Watts
    Ps = sum(Ps)
    Qs = sum(abs(Qs))
    efficiency_total = sum(efficiency)/3
    Spu_total = sum(Spu)/3
    PF_total = sum(PF)/3
    
    return [Pp,Qp,Ps,Qs,loss_total,efficiency_total,Spu_total,PF_total]
