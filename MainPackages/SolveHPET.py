#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy import interpolate
#import cmath

def SolveHPET(Main_path,t,np,math,Vcontrol,DSSObj,Loads_LV,Load_phase,TransformerRating,FractionPowerHPET):

    DSSText = DSSObj.Text # to excecute OpenDSS text commands
    DSSCircuit = DSSObj.ActiveCircuit # Use it to access the elements in the circuit (e.g., capacitors, buses, etc.)
    DSSSolution = DSSCircuit.Solution
    DSSMonitors = DSSCircuit.Monitors
    
    # Whole network variables
    Vloads_LV = np.full([Loads_LV.size,3],np.nan) # LV load voltages split per phase
    Vloads_min = np.zeros(3)    
    
    # Whole HPET variables
    V_LV = np.zeros((3,2)) # voltages,angles for each phase of the transformer's LV terminal
    V_LV_cmplx = np.zeros(3,dtype=np.complex_) # V_LV expressed as a complex in rectangular form
    I_LV_cmplx = np.zeros(3,dtype=np.complex_) # V_LV expressed as a complex in rectangular form
    Vtarget_cmplx = np.zeros(3,dtype=np.complex_) # Nominal values of V_LV expressed as a complex in rectangular form
    Vwin_nom = np.zeros(3)
    
    # Converter variables
    Vc2_pu = np.zeros((3,2)) # Delta voltages,angles in pu for each phase of the STvsources
    Vc2_cmplx = np.zeros(3,dtype=np.complex_) # V_LV expressed as a complex in rectangular form
    Pc2 = np.zeros(3)
    Qc2 = np.zeros(3)
    Sdc2 = np.zeros(3)
    Pdc = np.zeros(3)
    Pc1 = np.zeros(3)
    Qc1 = np.zeros(3)
    Sc1 = np.zeros(3)
    V_Seq = np.zeros(3,dtype=np.complex_)
    I_Seq = np.zeros(3,dtype=np.complex_)
    Qpos = np.zeros(3)
    Qneg = np.zeros(3)
    Qavail = np.zeros(3)

    #DeltaV_pu = np.zeros(3)

    # Get windings' nominal voltages
    DSSCircuit.Transformers.Name = 'Trafo1'
    DSSCircuit.Transformers.Wdg = 1
    Vwin_nom[0] = DSSCircuit.Transformers.kV * 1000
    DSSCircuit.Transformers.Wdg = 2
    Vwin_nom[1] = DSSCircuit.Transformers.kV * 1000
    DSSCircuit.Transformers.Wdg = 3
    Vwin_nom[2] = DSSCircuit.Transformers.kV * 1000
    
    # Create variable for STvsources nominal voltage (Base kV)
    Vc2_nom = Vwin_nom[1]*FractionPowerHPET
    Sc2_nom = (TransformerRating*FractionPowerHPET)/(1+FractionPowerHPET) # Per-phase PEM's nominal apparent power in kVA
    Sc1_nom  = Sc2_nom 
    
    # Set output voltage limits
    Mv = 0.002 # Security margin for minimum voltage regulation
    V_LV_min = (0.9 + Mv)*Vwin_nom[1] # Min statutory LV voltage limit plus security margin Mv
    V_LV_max = (1.1 - Mv)*Vwin_nom[1] # Max statutory LV voltage limit minus security margin Mv
    Vtarget = np.zeros(3)
    VregFlag = True
        
    Iter = 0
    Vtest = np.ones(3)
    Qtest = np.ones(3)
    Verr = np.ones(3)
    Qerr = np.ones(3)

    # Efficiency curves and interpolation function
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_MatlabLossModel_UnityPF.txt","rb"),delimiter=",")
    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_KACOblueplanet125TL3.txt","rb"),delimiter=",")
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_XantrexGT3.8Inverter.txt","rb"),delimiter=",")
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_Qin2010.txt","rb"),delimiter=",")
#    eff_points = np.loadtxt(open(Main_path + "/NetworkData/NetworkModels/Transformers/EffCurve_Adabi2018.txt","rb"),delimiter=",")
#    eff_points[:,1] = (0.992/np.nanmax(eff_points[:,1])) * eff_points[:,1] # offset efficiency to a different maximum
    eff_interpFunction = interpolate.interp1d(eff_points[:,0],eff_points[:,1], kind='linear', fill_value="extrapolate")
    efficiency1 = np.zeros(3)
    efficiency2 = np.zeros(3)
    
    
    # Pasivate STvsources
    for phase in range(3):
        DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
        DSSCircuit.Vsources.pu = 1e-3

        
    # 1st SOLVE: get actual demand
    DSSText.Command ='Reset Monitors'
    DSSSolution.Solve()
    if not(DSSSolution.Converged):
        raise ValueError('Solution did not Converge')
    DSSMonitors.SampleAll() 
    DSSMonitors.SaveAll()

    
    
    while (np.max(Qerr) > 0.05 or np.max(Verr) > 0.05) and Iter < 10:
        
        # Read secondary winding and STvsources actual voltages and angles per phase
        DSSMonitors.Name = 'Trafo_sec-VI'
        for phase in range(3):
            V_LV[phase,0] = DSSMonitors.Channel(phase * 2 + 1)[0]
            V_LV[phase,1] = DSSMonitors.Channel(phase * 2 + 2)[0]
        
                
        if Vcontrol[0] == 'Vfixed':
            Vtarget.fill(Vcontrol[1][0]) # Voltage reference [pu]
    
    
        if Vcontrol[0] == 'Vband':
            for phase in range(3):
                if V_LV[phase,0] < Vcontrol[1][0]*Vwin_nom[1]: # If V_LV is lower the allowed band, regulate voltage using the lower Vband
                    Vtarget.fill(Vcontrol[1][0])
                elif V_LV[phase,0] > Vcontrol[1][1]*Vwin_nom[1]: # If V_LV is greater the allowed band, regulate voltage using the upper Vband
                    Vtarget.fill(Vcontrol[1][1])
                else:
                    VregFlag = False # If V_LV is inside the allowed band, keep STvsources passivated
    
    
        if Vcontrol[0] == 'V_CVR': # the STvsources voltage is regulated to obtain the minimum possible demand depending on the voltage sensitivity
            Pagregated = 0
            ZIPagregated = np.zeros(3)
    
            # Get all LV loads voltage, per phase
            for i_load in range(Loads_LV.size):
                DSSMonitors.Name = 'Monitor_' + Loads_LV[i_load] + '_vi'
                Vloads_LV[i_load,Load_phase.at[Loads_LV[i_load],'phase']-1] = DSSMonitors.Channel(1)[0]
                DSSCircuit.Loads.Name = Loads_LV[i_load]
                Pagregated += DSSCircuit.Loads.kW # Adds al the load power for the next calculation
    
            # Obtains voltage sensitivity VS by computing agregated load: equations (17) and (23) from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7725539
            for i_load in range(Loads_LV.size):
                DSSCircuit.Loads.Name = Loads_LV[i_load]
                #print(Loads_LV[i_load])
                ZIP = DSSCircuit.Loads.ZIPV[[0,1,2]] # Obtains only the first 3 components which correspond to the active power coefficients
                ZIPagregated += (DSSCircuit.Loads.kW/Pagregated)*ZIP
                VS = (2*ZIPagregated[0] + ZIPagregated[1] + ZIPagregated[2])/sum(ZIPagregated)
    
            if VS > 0:
                Vloads_min = np.nanmin(Vloads_LV,0)
                Vtarget = (V_LV[:,0] - (Vloads_min - V_LV_min))/Vwin_nom[1]
            elif VS < 0:
                Vloads_max = np.nanmax(Vloads_LV,0)
                Vtarget = (V_LV[:,0] + (Vloads_max - V_LV_max))/Vwin_nom[1]
            else:
                Vtarget = V_LV[:,0]/Vwin_nom[1]
        

        # Voltage regulation calculation
        if VregFlag == True:
            for phase in range(3):
                V_LV_cmplx[phase] = V_LV[phase,0] * np.exp(1j*math.radians(V_LV[phase,1]))
                Vtarget_cmplx[phase] = Vtarget[phase] * Vwin_nom[1] * np.exp(1j*math.radians(30 - 120.0*phase))
                DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
                Vc2_cmplx[phase] = DSSCircuit.Vsources.pu * Vc2_nom * np.exp(1j*math.radians(DSSCircuit.Vsources.AngleDeg))
            
            Vc2_cmplx = Vtarget_cmplx - V_LV_cmplx + Vc2_cmplx
            Vc2_pu[:,0] = np.absolute(Vc2_cmplx)/Vc2_nom
            Vc2_pu[:,1] = np.degrees(np.angle(Vc2_cmplx))
    
            # Update STvsources pu and agnle parameters
            for phase in range(3):
                DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
                if abs(Vc2_pu[phase,0]) > 1.0: # Voltage limitation of the BtB converter
                    DSSCircuit.Vsources.pu = 1.0
                    print('Voltage regulation exceeded ' + str(Vc2_pu[phase,0]))
                else:
                    DSSCircuit.Vsources.pu = abs(Vc2_pu[phase,0])
                # The angle defines the sign of the compensation
                DSSCircuit.Vsources.AngleDeg = Vc2_pu[phase,1]    
    
        
            # 3rd SOLVE: final result for the current time step
            DSSText.Command ='Reset Monitors'
            DSSSolution.Solve()
            if not(DSSSolution.Converged):
                raise ValueError('Solution did not Converge')
            DSSMonitors.SampleAll() 
            DSSMonitors.SaveAll()
        
    
        # Efficiency calculation
        
        # Read Pc2 and Qc2 at the converter's output per phase. Those powers are measured at Vsource's output, thus they are negative
        DSSMonitors.Name = 'STvsources-PQ'
        for phase in range(3):
            Pc2[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            Qc2[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
    
        # Normalization of S, per phase
        Sc2 = np.sqrt(Pc2**2 + Qc2**2)
        Sc2_pu = Sc2/Sc2_nom
        PFc2 = abs(Pc2/Sc2)
        
        # Calculate efficiency for Module 2 and apparent power in the DC-link, per phase
        for phase in range(3):
            efficiency2[phase] = eff_interpFunction(Sc2_pu[phase])
        
        if Pc2[0] >= 0:
            Sdc2 = Sc2 / efficiency2 # forward power flow
            Ploss2 = Sdc2 - Sc2
            Pdc = Pc2 + Ploss2
        else:
            Sdc2 = Sc2 * efficiency2 # reverse power flow
            Ploss2 = Sdc2 - Sc2
            Pdc = Pc2 - Ploss2
            
        Sdc1 = np.sqrt(Pdc**2 + Qc1**2)
        Sdc1_pu = Sdc1/Sc1_nom
        PFc1 = abs(Pc1/Sdc1)
        
        # Calculate efficiency for Module 1 and input apparent power, per phase
        for phase in range(3):
            efficiency1[phase] = eff_interpFunction(Sdc1_pu[phase])
        
        if Pdc[0] >= 0:
            Sc1 = Sdc1 / efficiency1 # forward power flow
            Ploss1 = Sc1 - Sdc1
            Pc1 = Pdc + Ploss1
        else:
            Sc1 = Sdc1 * efficiency1 # forward power flow
            Ploss1 = Sc1 - Sdc1
            Pc1 = Pdc - Ploss1

        # Update STloads values
        for phase in range(3):
            DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kW=' + str(Pc1[phase]) + ' kvar=' + str(Qc1[phase])
  

            
        
        # 2nd SOLVE
        DSSText.Command ='Reset Monitors'
        DSSSolution.Solve()
        if not(DSSSolution.Converged):
            raise ValueError('Solution did not Converge')
        DSSMonitors.SampleAll()
        DSSMonitors.SaveAll()    
        


        # Reactive power compensation for the unbalanced voltage case

        # Q compensation measuring secondary side
        DSSMonitors.Name = 'STvsources-PQ'
        for phase in range(3):
            Pc2[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            Qc2[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
    
        # Unbalanced reactive power compensation calculations with unbalanced voltages
        # See Blasco, P. et al. Compensation of reactive power and unbalanced power in three-phase three-wire systems connected to an infinite power network (Year 2020)
        DSSMonitors.Name = 'Trafo_sec-VI'
        for phase in range(3):
            V_LV_cmplx[phase] = DSSMonitors.Channel(phase * 2 + 1)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 2)[0]))
            I_LV_cmplx[phase] = DSSMonitors.Channel(phase * 2 + 7)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 8)[0]))
        
        V_LVpp = np.array(((V_LV_cmplx[0]-V_LV_cmplx[1]),(V_LV_cmplx[1]-V_LV_cmplx[2]),(V_LV_cmplx[2]-V_LV_cmplx[0])),dtype=np.complex_)
        
        # Calculate symetrical sequence components from phase voltages and currents (the components for phase 2 and 3 must be obtained by using the properties of the positive-, negative-, and zero-sequence sets which they represent (Bergen1999, pag 448))
        a = np.exp(1j*math.radians(120.0))
        b = math.sqrt(3) * np.exp(1j*math.radians(30.0))
        Ainv = (1/3) * np.array(([1,1,1],[1,a,a**2],[1,a**2,a]),dtype=np.complex_)
        
        V_Seq = np.dot(Ainv,V_LVpp)
        I_Seq = np.dot(Ainv,I_LV_cmplx)        
        V_s1 = np.array((V_Seq[1], a**2*V_Seq[1], a*V_Seq[1]),dtype=np.complex_)/b
        V_s2 = np.array((V_Seq[2], a*V_Seq[2], a**2*V_Seq[2]),dtype=np.complex_)/np.conjugate(b)
        I_s1 = np.array((I_Seq[1], a**2*I_Seq[1], a*I_Seq[1]),dtype=np.complex_)
        I_s2 = np.array((I_Seq[2], a*I_Seq[2], a**2*I_Seq[2]),dtype=np.complex_)
        
        Qpos = np.abs(V_s1) * np.abs(I_s1) * np.sin(np.angle(V_s1) - np.angle(I_s1))
        Qneg = np.abs(V_s2) * np.abs(I_s2) * np.sin(np.angle(V_s2) - np.angle(I_s2))
        
        delta = np.abs(V_s2[0])/np.abs(V_s1[0])
        Qpos_comp = (Qneg - Qpos)/(1 - delta**2)
        Icomp_neg = I_s2 + V_s2/(1j*(np.absolute(V_s1)**2/Qpos_comp))
        Qnegpos_comp = -2*np.abs(V_s1) * np.abs(Icomp_neg) * np.sin(np.angle(V_s1) - np.angle(Icomp_neg))
        
        Qt = -Qpos_comp/1000 - Qnegpos_comp/1000 - Qc2

            
        for phase in range(3):
            if abs(Pc2[phase]) > Sc2_nom:
                Qavail[phase] = 0
            else:
                Qavail[phase] = math.sqrt(Sc2_nom**2 - Pc2[phase]**2)   
                if abs(Qt[phase]) > Qavail[phase]:
                    Qc1[phase] = -(abs(Qt[phase])/Qt[phase]) * Qavail[phase] # Take the actual sign of Q_LV[phase] and invert it for compensation
                    print('Q compensation exceeded on phase ' + str(phase+1))
                else:
                    Qc1[phase] = -Qt[phase] # Total compensation: Qp inyects directly to the transformer's core
        
        
        # Efficiency calculation
        
        # Read Pc2 and Qc2 at the converter's output per phase. Those powers are measured at Vsource's output, thus they are negative
        DSSMonitors.Name = 'STvsources-PQ'
        for phase in range(3):
            Pc2[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            Qc2[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
    
        # Normalization of S, per phase
        Sc2 = np.sqrt(Pc2**2 + Qc2**2)
        Sc2_pu = Sc2/Sc2_nom
        PFc2 = abs(Pc2/Sc2)
        
        # Calculate efficiency for Module 2 and apparent power in the DC-link, per phase
        for phase in range(3):
            efficiency2[phase] = eff_interpFunction(Sc2_pu[phase])
        
        if Pc2[0] >= 0:
            Sdc2 = Sc2 / efficiency2 # forward power flow
            Ploss2 = Sdc2 - Sc2
            Pdc = Pc2 + Ploss2
        else:
            Sdc2 = Sc2 * efficiency2 # reverse power flow
            Ploss2 = Sdc2 - Sc2
            Pdc = Pc2 - Ploss2

        Sdc1 = np.sqrt(Pdc**2 + Qc1**2)
        Sdc1_pu = Sdc1/Sc1_nom
        PFc1 = abs(Pc1/Sdc1)
        
        # Calculate efficiency for Module 1 and input apparent power, per phase
        for phase in range(3):
            efficiency1[phase] = eff_interpFunction(Sdc1_pu[phase])
        
        if Pdc[0] >= 0:
            Sc1 = Sdc1 / efficiency1 # forward power flow
            Ploss1 = Sc1 - Sdc1
            Pc1 = Pdc + Ploss1
        else:
            Sc1 = Sdc1 * efficiency1 # forward power flow
            Ploss1 = Sc1 - Sdc1
            Pc1 = Pdc - Ploss1

        # Update STloads values
        for phase in range(3):
            DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kW=' + str(Pc1[phase]) + ' kvar=' + str(Qc1[phase])
            
            
        # 3rd SOLVE
        DSSText.Command ='Reset Monitors'
        DSSSolution.Solve()
        if not(DSSSolution.Converged):
            raise ValueError('Solution did not Converge')
        DSSMonitors.SampleAll()
        DSSMonitors.SaveAll()  

        
        for phase in range(3):
            DSSMonitors.Name = 'Trafo_sec-VI'
            Verr[phase] = abs((DSSMonitors.Channel(phase * 2 + 1)[0] - Vtest[phase])/Vtest[phase])
            Vtest[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            DSSMonitors.Name = 'Trafo_prim-PQ'
            Qerr[phase] = abs((DSSMonitors.Channel(phase * 2 + 2)[0] - Qtest[phase])/Qtest[phase])
            Qtest[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
            
        Iter += 1

    
    # Calculate total values and exit
    loss_total = sum(abs(Pc1 - Pc2)) # Losses in kWatts
    efficiency1_avg = sum(efficiency1)/3
    efficiency2_avg = sum(efficiency2)/3
    Sdc1_pu_total = sum(Sdc1_pu)/3
    Sc2_pu_total = sum(Sc2_pu)/3
    PFc1_total = sum(PFc1)/3
    PFc2_total = sum(PFc2)/3
    
    return [Pc1/Sc1_nom, Qc1/Sc1_nom, Pc2/Sc2_nom, Qc2/Sc2_nom, loss_total, efficiency1_avg, efficiency2_avg, Sdc1_pu_total, Sc2_pu_total, PFc1_total, PFc2_total]
