#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

def Extract_data_OpenDSS(math,np,pd,DSSObj,model):
    
    DSSText = DSSObj.Text # to excecute OpenDSS text commands
    DSSCircuit = DSSObj.ActiveCircuit # Use it to access the elements in the circuit (e.g., capacitors, buses, etc.)
    DSSElement = DSSCircuit.ActiveCktElement
    DSSBus = DSSCircuit.ActiveBus
    
    # special circuit elements sets (upstream network elements)
    Lines_upstream = np.array(['mv_line'])
    Buses_upstream = np.array(['primary','secondary'])
    
    
    # downstream network elements sets
    Lines_set_all = np.array(DSSCircuit.Lines.AllNames)
    Lines_set = Lines_set_all

    if model == 'HT':
        Lines_HT = np.array(['htline1','htline2','htline3','htline4'])
        for i_line in range(Lines_HT.size):
            Lines_set = np.delete(Lines_set, np.where(Lines_set == Lines_HT[i_line])[0][0])
    
    Lines_downstream = Lines_set
    for i_line in range(Lines_upstream.size):
        Lines_downstream = np.delete(Lines_downstream, np.where(Lines_downstream == Lines_upstream[i_line])[0][0])
    
    for i_line in Lines_set:
        DSSText.Command = 'New Monitor.' + i_line + '_VI_sending Line.' + i_line + ' Terminal=1 Mode=0 VIpolar=yes'
        DSSText.Command = 'New Monitor.' + i_line + '_VI_receiving Line.' + i_line + ' Terminal=2 Mode=0 VIpolar=yes'
        DSSText.Command = 'New Monitor.' + i_line + '_PQ_sending Line.' + i_line + ' Terminal=1 Mode=1 ppolar=no'
        DSSText.Command = 'New Monitor.' + i_line + '_PQ_receiving Line.' + i_line + ' Terminal=2 Mode=1 ppolar=no'
        
    
    Buses_set = np.array(DSSCircuit.AllBusNames)
    Buses_downstream = np.delete(Buses_set, np.where(Buses_set == 'slack')[0][0])
    for i_bus in range(Buses_upstream.size):
        Buses_downstream = np.delete(Buses_downstream, np.where(Buses_downstream == Buses_upstream[i_bus])[0][0])
    
    Nodes_set = np.array(DSSCircuit.AllNodeNames)


    # Get buses nominal voltages
    Bus_Vnom = pd.DataFrame(0,index=Buses_set,columns=['Vnom_pp','Vnom_pn'])
    for i_bus in Buses_set:
        DSSCircuit.SetActiveBus(i_bus)
        Bus_Vnom.loc[i_bus,'Vnom_pp'] = DSSBus.kVBase*1000*math.sqrt(3)
        Bus_Vnom.loc[i_bus,'Vnom_pn'] = DSSBus.kVBase*1000

    
    # Get line downstream data and create line monitors
    Lines_downstream_data = pd.DataFrame(index=Lines_downstream, columns=['Sending bus','Receiving bus','Cable code'])
    for i_line in Lines_downstream:
        DSSCircuit.Lines.Name = i_line
        Lines_downstream_data.loc[i_line,'Sending bus']=DSSCircuit.Lines.Bus1
        Lines_downstream_data.loc[i_line,'Receiving bus']=DSSCircuit.Lines.Bus2
        Lines_downstream_data.loc[i_line,'Cable code']=DSSCircuit.Lines.LineCode

    # Get line upstream data and create line monitors
    Lines_upstream_data = pd.DataFrame(index=Lines_upstream, columns=['Sending bus','Receiving bus','Cable code'])
    for i_line in Lines_upstream:
        DSSCircuit.Lines.Name = i_line
        Lines_upstream_data.loc[i_line,'Sending bus']=DSSCircuit.Lines.Bus1
        Lines_upstream_data.loc[i_line,'Receiving bus']=DSSCircuit.Lines.Bus2
        #Lines_upstream_data.loc[i_line,'Cable code']=DSSCircuit.Lines.LineCode
        
    # Get load data
    Loads_set_all = np.array(DSSCircuit.Loads.AllNames)
    Loads_set = Loads_set_all
    for i_load in range(Loads_set_all.size):
        if re.fullmatch('feeder.*', Loads_set_all[i_load]):
            pass
        else:
            Loads_set = np.delete(Loads_set, np.where(Loads_set == Loads_set_all[i_load])[0][0])
    
    Load_phase = pd.DataFrame(0,index=Loads_set, columns=["phase"]) 
    Load_bus = pd.DataFrame(index=Loads_set, columns=["Bus"]) # you cannot mix string with floats
    Load_Vnom = pd.DataFrame(index=Loads_set, columns=["Vnom [V]"]) # you cannot mix string with floats
    
    for i_load in Loads_set: # Loop thought all loads in the system  and create load monitors
        DSSText.Command = 'New Monitor.Monitor_' + i_load + '_VI Load.' + i_load + ' Terminal=1 Mode=0 VIpolar=yes'
        DSSText.Command = 'New Monitor.Monitor_' + i_load + '_PQ Load.' + i_load + ' Terminal=1 Mode=1 ppolar=no' 
        DSSCircuit.SetActiveElement("Load." + i_load) # Defines the active circuit element. Allows to access some properties
        DSSCircuit.Loads.Name = i_load # Active load element in DSS
        if DSSElement.NumPhases > 1:
            raise ValueError('not single phase load')
        
        # Get connection phases
        BusNames = DSSElement.BusNames[0] # Get the bus connection of this load
        Load_bus.loc[i_load]["Bus"] = BusNames.partition(".")[0]
        Phases_name = BusNames.partition(".")[2]
        Load_phase.loc[i_load]["phase"] = int(Phases_name)
        Load_Vnom.loc[i_load]["Vnom [V]"] = DSSCircuit.Loads.kV * 1000


    return [Lines_set,Lines_upstream,Lines_downstream,Lines_upstream_data,Lines_downstream_data,Buses_set,Buses_upstream,Buses_downstream,Bus_Vnom,Nodes_set,Loads_set,Load_phase,Load_bus,Load_Vnom]
