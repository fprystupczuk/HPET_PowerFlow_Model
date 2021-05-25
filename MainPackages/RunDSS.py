#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from MainPackages.SolvePET import SolvePET
from MainPackages.SolveHPET import SolveHPET
from MainPackages.Extract_DSSmonitors_data import Extract_DSSmonitors_data

def RunDSS(os,np,pd,math,Main_path,model,Network,Feeders,Vcontrol,PV_penetration,ElectricVars,Profiles_path,Main_Results_path,timeArray,DSSObj,Lines_set,Loads_LV,Load_bus,Load_phase,Bus_Vnom,MV_loads,TransformerRating,PowerRatingPET,FractionPowerHPET):

    DSSText = DSSObj.Text # to excecute OpenDSS text commands
    DSSCircuit = DSSObj.ActiveCircuit # Use it to access the elements in the circuit (e.g., capacitors, buses, etc.)
    DSSElement = DSSCircuit.ActiveCktElement
    DSSSolution = DSSCircuit.Solution
    DSSMonitors = DSSCircuit.Monitors

    """ DEMAND AND PV ALLOCATION PROFILES
    As written, the script randomly allocates demand and PV active power profiles 
    to the network using the profiles provided in NetworkData/Profiles """
    
    # MV Demand
    MV_Dem_data = {} # defined as dictionary
    for i in range(3):
        i_MV = 'MV_' + str(i+1)
        csv_data = pd.read_csv(Profiles_path + '/Load_profiles/Profile_MV_PZIP_QZIP.csv',index_col=0)
        MV_Dem_data[i_MV,'P_Profile'] = csv_data['P_' + str(i+1)] # Active power profile [W]
        MV_Dem_data[i_MV,'Q_Profile'] = csv_data['Q_' + str(i+1)] # Reactive power profile [var]
        MV_Dem_data[i_MV,'ZIP_P'] = csv_data[['P_Z','P_I','P_P']]
        MV_Dem_data[i_MV,'ZIP_Q'] = csv_data[['Q_Z','Q_I','Q_P']]
        MV_Dem_data[i_MV,'ZIP_P'].columns = ['Z','I','P']
        MV_Dem_data[i_MV,'ZIP_Q'].columns = ['Z','I','P']
    
    
    # LV Demand
    np.random.seed(0) # for repeatability
    House_Dem_data = {} # defined as dictionary
    h_profiles = np.random.randint(1,100,size=Loads_LV.size) # define which random profiles will be allocated to each customer
    for i_house in range(Loads_LV.size): 
        csv_data = pd.read_csv(Profiles_path + '/Load_profiles/Profile_' + str(h_profiles[i_house]) + '_PZIP_QZIP.csv',index_col=0)    
        House_Dem_data[Loads_LV[i_house],'P_Profile'] = csv_data['P'] # Active power profile [W]
        House_Dem_data[Loads_LV[i_house],'Q_Profile'] = csv_data['Q'] # Reactive power profile [var] - here, assumed a PF of
        House_Dem_data[Loads_LV[i_house],'ZIP_P'] = csv_data[['P_Z','P_I','P_P']]
        House_Dem_data[Loads_LV[i_house],'ZIP_Q'] = csv_data[['Q_Z','Q_I','Q_P']]    
        House_Dem_data[Loads_LV[i_house],'ZIP_P'].columns = ['Z','I','P']
        House_Dem_data[Loads_LV[i_house],'ZIP_Q'].columns = ['Z','I','P']
    
    
    # PV generation
    np.random.seed(1) # for repeatability
    random.seed(1)
    PV_Gen_data = {} # defined as dictionary
    PV_set = np.array(random.sample(list(Loads_LV), int(np.rint(len(Loads_LV) * PV_penetration)))) # Random PV alocation
    PV_profiles = np.random.randint(1,100,size=PV_set.size) # Random selection of PV profiles from the pool
    PV_profiles_data = pd.read_excel(Profiles_path + '/PV_profiles/Summer_PV_Profiles.xlsx',header=0,index_col=0) # Read selected profiles (data in kW)
    PV_rating_data = pd.read_csv(Profiles_path + '/PV_profiles/PV_rating.csv')
    for i_PV in range(PV_set.size): 
        PV_Gen_data[PV_set[i_PV],'Profile'] = PV_profiles_data[PV_profiles[i_PV]]*1000 # Active power output for this PV [W]
        PV_Gen_data[PV_set[i_PV],'Rating'] = PV_rating_data.loc[PV_profiles[i_PV],'Rating [kW]'] # PV rating in kW
    
    # Create PVs in OpenDSS
    for PV in PV_set:
        DSSText.Command = 'New Load.PV_' + PV + ' Phases=1 Bus1=' + Load_bus.loc[PV,'Bus'] + '.' + str(Load_phase.loc[PV,'phase']) + ' kV=' + str(0.4/math.sqrt(3)) + ' kW=0 Kvar=0 Model=1 Vminpu=0.7 Vmaxpu=1.3'

    
    """ 
    RUN POWER FLOW AND OBTAIN RESULTS
    """
    
    # OUTPUT DATA CONTAINERS    
    Data = np.zeros((timeArray.size,ElectricVars.size,Lines_set.size,3))
    ## Order of ElectricVars in Data array ##
    # 0 <- Vmag_send
    # 1 <- Vang_send
    # 2 <- Vmag_rec
    # 3 <- Vang_rec
    # 4 <- Imag
    # 5 <- P_rec
    # 6 <- Q_rec
    # 7 <- P_send
    # 8 <- Q_send

    # Name of variables to be collected for not per line but only for each time-step    
    P_LV = np.zeros(timeArray.size) # P delivered by the transformer's secondary side (sum of phases)
    Q_LV = np.zeros(timeArray.size) # Q delivered by the transformer's secondary side (sum of phases)    
    TotalLosses = np.zeros((timeArray.size,2)) # Line losses + Transformer losses
    
    LFTlosses = np.zeros(timeArray.size) # Line losses + Transformer losses
    
    PETlosses = np.zeros(timeArray.size) # Losses of the Power Electronic Module
    PETefficiency = np.zeros(timeArray.size) # Efficiency of the Power Electronic Module
    PET_Spu = np.zeros(timeArray.size) # Apparent power of the Power Electronic Module in pu
    PET_PF_LV = np.zeros(timeArray.size) # Power Factor at the transformer's secondary side

    PEM_Pc1_pu = np.zeros((timeArray.size,3)) # Active power of Module 1 of the HPET in pu
    PEM_Qc1_pu = np.zeros((timeArray.size,3)) # Reactive power of Module 1 of the HPET in pu
    PEM_Pc2_pu = np.zeros((timeArray.size,3)) # Active power of Module 2 of the HPET in pu
    PEM_Qc2_pu = np.zeros((timeArray.size,3)) # Reactive power of Module 2 of the HPET in pu    
    HPETlosses = np.zeros((timeArray.size,2)) # HPET losses: PEM losses,Transformer losses
    PEMefficiency1 = np.zeros(timeArray.size) # Efficiency of Module 1 of the HPET
    PEMefficiency2 = np.zeros(timeArray.size) # Efficiency of Module 2 of the HPET
    PEM_Sc1_pu = np.zeros(timeArray.size) # Apparent power of Module 1 of the HPET in pu
    PEM_Sc2_pu = np.zeros(timeArray.size) # Apparent power of Module 2 of the HPET in pu
    PEM_PFc1 = np.zeros(timeArray.size) # Power Factor at Module 1 of the HPET
    PEM_PFc2 = np.zeros(timeArray.size) # Power Factor at Module 2 of the HPET

    # Path to store results
    fileName_start = Main_Results_path + '/DataResults/' + Network + '/'
    try: os.makedirs(fileName_start)
    except FileExistsError: pass



    # Run time series of Power Flow analysis            
    for t in range(timeArray.size):
        
        # OpenDSS MV loads
        if MV_loads == 'y':
            for i in range(3):
                i_MV = 'MV_' + str(i+1)
                DSSCircuit.Loads.Name = i_MV
                DSSCircuit.Loads.kW = MV_Dem_data[i_MV,'P_Profile'].loc[timeArray[t]]/1000.0
                DSSCircuit.Loads.kvar = MV_Dem_data[i_MV,'Q_Profile'].loc[timeArray[t]]/1000.0       
                Z_p = float(MV_Dem_data[i_MV,'ZIP_P'].loc[timeArray[t],'Z'])
                I_p = float(MV_Dem_data[i_MV,'ZIP_P'].loc[timeArray[t],'I'])
                Z_q = float(MV_Dem_data[i_MV,'ZIP_Q'].loc[timeArray[t],'Z'])
                I_q = float(MV_Dem_data[i_MV,'ZIP_Q'].loc[timeArray[t],'I'])
                DSSCircuit.Loads.ZIPV = (Z_p,I_p,1-Z_p-I_p,Z_q,I_q,1-Z_q-I_q,0.8) # Last coefficient: voltage in pu from wich the load model changes to constant impedance to facilitate convergency of OpenDSS
        
        # OpenDSS LV feeder house demand
        for i_house in Loads_LV:
            DSSCircuit.Loads.Name = i_house
            DSSCircuit.Loads.kW = House_Dem_data[i_house,'P_Profile'].loc[timeArray[t]]/1000.0
            DSSCircuit.Loads.kvar = House_Dem_data[i_house,'Q_Profile'].loc[timeArray[t]]/1000.0       
            Z_p = float(House_Dem_data[i_house,'ZIP_P'].loc[timeArray[t],'Z'])
            I_p = float(House_Dem_data[i_house,'ZIP_P'].loc[timeArray[t],'I'])
            Z_q = float(House_Dem_data[i_house,'ZIP_Q'].loc[timeArray[t],'Z'])
            I_q = float(House_Dem_data[i_house,'ZIP_Q'].loc[timeArray[t],'I'])
            DSSCircuit.Loads.ZIPV = (Z_p,I_p,1-Z_p-I_p,Z_q,I_q,1-Z_q-I_q,0.8) # Last coefficient: voltage in pu from wich the load model changes to constant impedance to facilitate convergency of OpenDSS
    
        # OpenDSS PVs
        for PV in PV_set:
            DSSCircuit.Loads.Name = 'PV_' + PV
            DSSCircuit.Loads.kW = -1 * PV_Gen_data[PV,'Profile'].loc[timeArray[t]]/1000.0
            DSSCircuit.Loads.kvar = -1 * 0.0 # Assumes PF=1 for all the PVs        

        # Solve time step for different transformer models
        if model == 'LFT':
                
            # SOLVE
            DSSText.Command ='Reset Monitors'
            DSSSolution.Solve()
            if not(DSSSolution.Converged):
                raise ValueError('Solution did not Converge')        
            DSSMonitors.SampleAll() 
            DSSMonitors.SaveAll()
            
            DSSMonitors.Name = 'LFT-PQ'
            for phase in range(3):
                P_LV[t] = P_LV[t] - DSSMonitors.Channel(phase * 2 + 1)[0]
                Q_LV[t] = Q_LV[t] + abs(DSSMonitors.Channel(phase * 2 + 2)[0])

            # Compute transformer losses and total losses            
            DSSCircuit.Transformers.Name = 'TR1'
            LFTlosses[t] = DSSElement.Losses[0]/1000 # DSSElement.Losses returns Watts and VAr
            TotalLosses[t,:] = np.asarray(DSSCircuit.Losses)/1000 # DSSCircuit.Losses returns Watts and VAr

            

        if model == 'PET':
            # Regulate output voltage, compute losses and update active power in the primary side
            [Pp,Qp,P_LV[t],Q_LV[t],PETlosses[t],PETefficiency[t],PET_Spu[t],PET_PF_LV[t]] = SolvePET(Main_path,np,Vcontrol,DSSObj,Bus_Vnom,Loads_LV,Load_phase,TransformerRating)

            # Compute total losses and add power electronics losses
            TotalLosses[t,:] = np.asarray(DSSCircuit.Losses)/1000 # DSSCircuit.Losses returns Watts and VAr
            TotalLosses[t,0] = TotalLosses[t,0] + PETlosses[t]

            

        if model == 'HPET':
            # compute losses and update active power in the primary side
            [PEM_Pc1_pu[t,:],PEM_Qc1_pu[t,:],PEM_Pc2_pu[t,:],PEM_Qc2_pu[t,:],HPETlosses[t,0],PEMefficiency1[t],PEMefficiency2[t],PEM_Sc1_pu[t],PEM_Sc2_pu[t],PEM_PFc1[t],PEM_PFc2[t]] = SolveHPET(Main_path,t,np,math,Vcontrol,DSSObj,Loads_LV,Load_phase,TransformerRating,FractionPowerHPET)
            
            # Obtain losses at the LFT of the HPET
            for phase in range(3):
                DSSCircuit.Transformers.Name = 'Trafo' + str(phase+1)
                HPETlosses[t,1] = HPETlosses[t,1] + DSSElement.Losses[0]/1000 # HPET losses: PEM losses,Transformer losses
                
            # Compute total losses and add power electronics losses
            TotalLosses[t,:] = np.asarray(DSSCircuit.Losses)/1000 # DSSCircuit.Losses returns Watts and VAr
            TotalLosses[t,0] = TotalLosses[t,0] + HPETlosses[t,0] # HPET losses: PEM losses,Transformer losses

            
            # Obtains the total power delivered by the secondary winding of the transformer
            DSSMonitors.Name = 'Trafo_sec-PQ'
            for phase in range(3):
                P_LV[t] = P_LV[t] + DSSMonitors.Channel(phase * 2 + 1)[0]
                Q_LV[t] = Q_LV[t] + abs(DSSMonitors.Channel(phase * 2 + 2)[0])

        

        # Store main results into a numpy array
        Data[t,:,:,:] = Extract_DSSmonitors_data(np,ElectricVars,DSSMonitors,Lines_set)
    
    
    # Define final paths to store results for each transformer model
    if model == 'LFT':
        fileName_end = model + '_PV' + str(PV_penetration) + '_' + '-'.join(Feeders)
        np.save(fileName_start + 'LFTlosses_' + fileName_end,LFTlosses)

    if model == 'PET':
        fileName_end = model + str(PowerRatingPET) + '_PV' + str(PV_penetration) + '_' + str(Vcontrol[0]) + str(Vcontrol[1][0]) + '_' + str(Vcontrol[1][1]) + '_' + '-'.join(Feeders)
        np.save(fileName_start + 'PETlosses_' + fileName_end,PETlosses)
        np.save(fileName_start + 'PETefficiency_' + fileName_end,PETefficiency)
        np.save(fileName_start + 'PET_Spu_' + fileName_end,PET_Spu)
        np.save(fileName_start + 'PET_PF_LV_' + fileName_end,PET_PF_LV)

    if model == 'HPET':
        fileName_end = model + str(FractionPowerHPET) + '_PV' + str(PV_penetration) + '_' + str(Vcontrol[0]) + str(Vcontrol[1][0]) + '_' + str(Vcontrol[1][1]) + '_' + '-'.join(Feeders)
        np.save(fileName_start + 'HPETlosses_' + fileName_end,HPETlosses) # HPET losses: PEM losses,Transformer losses
        np.save(fileName_start + 'PEM_Pc1_pu_' + fileName_end,PEM_Pc1_pu)
        np.save(fileName_start + 'PEM_Qc1_pu_' + fileName_end,PEM_Qc1_pu)
        np.save(fileName_start + 'PEM_Pc2_pu_' + fileName_end,PEM_Pc2_pu)
        np.save(fileName_start + 'PEM_Qc2_pu_' + fileName_end,PEM_Qc2_pu)
        np.save(fileName_start + 'PEMefficiency1_' + fileName_end,PEMefficiency1)
        np.save(fileName_start + 'PEMefficiency2_' + fileName_end,PEMefficiency2)
        np.save(fileName_start + 'PEM_Sc1_pu_' + fileName_end,PEM_Sc1_pu)
        np.save(fileName_start + 'PEM_Sc2_pu_' + fileName_end,PEM_Sc2_pu)
        np.save(fileName_start + 'PEM_PFc1_' + fileName_end,PEM_PFc1)
        np.save(fileName_start + 'PEM_PFc2_' + fileName_end,PEM_PFc2)

    
    # Save final results into the destination npy files
    if MV_loads == 'y': fileName_end = fileName_end + '_MVloads'
    
    np.save(fileName_start + 'Data_' + fileName_end,Data)
    np.save(fileName_start + 'P_LV_' + fileName_end,P_LV)
    np.save(fileName_start + 'Q_LV_' + fileName_end,Q_LV)
    np.save(fileName_start + 'TotalLosses_' + fileName_end,TotalLosses)
    