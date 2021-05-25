#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def PlotDSS(os,np,Main_Results_path,Time,timeArray,models,Network,Feeders,Vcontrol,PV_penetration,Lines_set,SelectedLinesUS,SelectedLinesDS,MV_loads,Bus_Vnom,Lines_upstream_data,Lines_downstream_data,PowerRatingPET,FractionPowerHPET):
    
    figure_size = (12,6)
    font_size = 20
    plt.rcParams.update({'font.size': font_size})
    
    # Load saved data from files for each electricVarName and model
    fileName_start = Main_Results_path + '/DataResults/' + Network + '/'
    if MV_loads == 'y': fileName_end = '_PV' + str(PV_penetration) + '_' + '-'.join(Feeders) + '_MVloads.npy'
    
    if 'LFT' in models:
        fileName_end = '_PV' + str(PV_penetration) + '_' + '-'.join(Feeders)
        LFTlosses = np.load(fileName_start + 'LFTlosses_LFT' + fileName_end + '.npy')
        Data_LFT = np.load(fileName_start + 'Data_LFT' + fileName_end + '.npy')
        P_LV_LFT = np.load(fileName_start + 'P_LV_LFT' + fileName_end + '.npy')
        Q_LV_LFT = np.load(fileName_start + 'Q_LV_LFT' + fileName_end + '.npy')
        TotalLosses_LFT = np.load(fileName_start + 'TotalLosses_LFT' + fileName_end + '.npy')
    if 'PET' in models:
        fileName_end = '_PV' + str(PV_penetration) + '_' + Vcontrol[0] + str(Vcontrol[1][0]) + '_' + str(Vcontrol[1][1]) + '_' + '-'.join(Feeders)
        PETlosses = np.load(fileName_start + 'PETlosses_PET' + str(PowerRatingPET) + fileName_end + '.npy')
        PETefficiency = np.load(fileName_start + 'PETefficiency_PET' + str(PowerRatingPET) + fileName_end + '.npy')
        PET_Spu = np.load(fileName_start + 'PET_Spu_PET' + str(PowerRatingPET) + fileName_end + '.npy')
        Data_PET = np.load(fileName_start + 'Data_PET' + str(PowerRatingPET) + fileName_end + '.npy')
        P_LV_PET = np.load(fileName_start + 'P_LV_PET' + str(PowerRatingPET) + fileName_end + '.npy')
        Q_LV_PET = np.load(fileName_start + 'Q_LV_PET' + str(PowerRatingPET) + fileName_end + '.npy')
        TotalLosses_PET = np.load(fileName_start + 'TotalLosses_PET' + str(PowerRatingPET) + fileName_end + '.npy')
    if 'HPET' in models:
        fileName_end = '_PV' + str(PV_penetration) + '_' + Vcontrol[0] + str(Vcontrol[1][0]) + '_' + str(Vcontrol[1][1]) + '_' + '-'.join(Feeders)
        HPETlosses = np.load(fileName_start + 'HPETlosses_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEMefficiency1 = np.load(fileName_start + 'PEMefficiency1_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEMefficiency2 = np.load(fileName_start + 'PEMefficiency2_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEM_Pc1_pu = np.load(fileName_start + 'PEM_Pc1_pu_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEM_Qc1_pu = np.load(fileName_start + 'PEM_Qc1_pu_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEM_Pc2_pu = np.load(fileName_start + 'PEM_Pc2_pu_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEM_Qc2_pu = np.load(fileName_start + 'PEM_Qc2_pu_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEM_Sc1_pu = np.load(fileName_start + 'PEM_Sc1_pu_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        PEM_Sc2_pu = np.load(fileName_start + 'PEM_Sc2_pu_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        Data_HPET = np.load(fileName_start + 'Data_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        P_LV_HPET = np.load(fileName_start + 'P_LV_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        Q_LV_HPET = np.load(fileName_start + 'Q_LV_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
        TotalLosses_HPET = np.load(fileName_start + 'TotalLosses_HPET' + str(FractionPowerHPET) + fileName_end + '.npy')
    
    # Name of the folder to store all the results
    versus = ''
    for i in range(len(models)):
        if models[i] == 'LFT':
            versus += '_LFT'
        elif models[i] == 'PET':
            versus += '_PET' + str(PowerRatingPET)
        if models[i] == 'HPET':
            versus += '_HPET' + str(FractionPowerHPET)

    ResultsFolder = Main_Results_path + '/Plots/' + Network + '/PV' + str(PV_penetration) + '_' + Vcontrol[0] + str(Vcontrol[1][0]) + '_' + str(Vcontrol[1][1]) + versus + '_' + '-'.join(Feeders)
    try: os.makedirs(ResultsFolder)
    except FileExistsError: pass

    ## ElectricVars order in the Data array ##
    # 0 <- Vmag_send
    # 1 <- Vang_send
    # 2 <- Vmag_rec
    # 3 <- Vang_rec
    # 4 <- Imag
    # 5 <- P_rec
    # 6 <- Q_rec
    # 7 <- P_send
    # 8 <- Q_send
    
    # Energy calculation in kWh
    E_SubStation_LFT = 0
    E_SubStation_PET = 0
    E_SubStation_HPET = 0
    E_LFT_sec = 0
    E_PET_sec = 0
    E_HPET_sec = 0
    E_SystemLosses_LFT = 0
    E_SystemLosses_PET = 0
    E_SystemLosses_HPET = 0
    E_TrafoLosses_LFT = 0
    E_TrafoLosses_PET = 0
    E_TrafoLosses_HPET = 0
    for i in range(len(models)):
        for t in range(timeArray.size):
            if models[i] == 'LFT':
                E_SubStation_LFT += (Time.step/60.0) * (Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],0] + Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],1] + Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],2])
                E_LFT_sec += (Time.step/60.0) * P_LV_LFT[t]
                E_SystemLosses_LFT += (Time.step/60.0) * TotalLosses_LFT[t,0]
                E_TrafoLosses_LFT += (Time.step/60.0) * LFTlosses[t]
            if models[i] == 'PET':
                E_SubStation_PET += (Time.step/60.0) * (Data_PET[t,7,np.where(Lines_set=='mv_line')[0][0],0] + Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],1] + Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],2])
                E_PET_sec += (Time.step/60.0) * P_LV_PET[t]
                E_SystemLosses_PET += (Time.step/60.0) * TotalLosses_PET[t,0]
                E_TrafoLosses_PET += (Time.step/60.0) * PETlosses[t]
            if models[i] == 'HPET':
                E_SubStation_HPET += (Time.step/60.0) * (Data_HPET[t,7,np.where(Lines_set=='mv_line')[0][0],0] + Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],1] + Data_LFT[t,7,np.where(Lines_set=='mv_line')[0][0],2])
                E_HPET_sec += (Time.step/60.0) * P_LV_HPET[t]
                E_SystemLosses_HPET += (Time.step/60.0) * TotalLosses_HPET[t,0]
                E_TrafoLosses_HPET += (Time.step/60.0) * (HPETlosses[t,0] + HPETlosses[t,1])
    
    print('\t\t\tLFT\t\tPET\t\tHPET')
    print('E at substation [kWh]\t' + str("{:.1f}".format(E_SubStation_LFT)) + '\t\t' + str("{:.1f}".format(E_SubStation_PET)) + '(' + str("{:.2f}".format(E_SubStation_PET/E_SubStation_LFT)) + 'x)\t' + str("{:.1f}".format(E_SubStation_HPET)) + '(' + str("{:.2f}".format(E_SubStation_HPET/E_SubStation_LFT)) + 'x)')
    print('E at secondary [kWh]\t' + str("{:.1f}".format(E_LFT_sec)) + '\t\t' + str("{:.1f}".format(E_PET_sec)) + '(' + str("{:.2f}".format(E_PET_sec/E_LFT_sec)) + 'x)\t' + str("{:.1f}".format(E_HPET_sec)) + '(' + str("{:.2f}".format(E_HPET_sec/E_LFT_sec)) + 'x)')
    print('Trafo Losses [kWh]\t' + str("{:.1f}".format(E_TrafoLosses_LFT)) + '\t\t' + str("{:.1f}".format(E_TrafoLosses_PET)) + '(' + str("{:.2f}".format(E_TrafoLosses_PET/E_TrafoLosses_LFT)) + 'x)\t' + str("{:.1f}".format(E_TrafoLosses_HPET)) + '(' + str("{:.2f}".format(E_TrafoLosses_HPET/E_TrafoLosses_LFT)) + 'x)')
    print('System Losses [kWh]\t' + str("{:.1f}".format(E_SystemLosses_LFT)) + '\t\t' + str("{:.1f}".format(E_SystemLosses_PET)) + '(' + str("{:.1f}".format(E_SystemLosses_PET/E_SystemLosses_LFT)) + 'x)\t' + str("{:.1f}".format(E_SystemLosses_HPET)) + '(' + str("{:.1f}".format(E_SystemLosses_HPET/E_SystemLosses_LFT)) + 'x)\n')
    
    # LaTex Table
    print('\\begin{table}[h]')
    print('    \centering')
    print('    \\begin{tabular}[center]{lccc}')
    print('        \hline & \\textbf{LFT} & \\textbf{PET} & \\textbf{HPET}\\\\\hline')
    print('        Energy from Substation [kWh] & ' + str("{:.1f}".format(E_SubStation_LFT)) + ' & ' + str("{:.1f}".format(E_SubStation_PET)) + ' & ' + str("{:.1f}".format(E_SubStation_HPET)) + '\\\\')
    print('        Energy from Transformer [kWh] & ' + str("{:.1f}".format(E_LFT_sec)) + ' & ' + str("{:.1f}".format(E_PET_sec)) + ' & ' + str("{:.1f}".format(E_HPET_sec)) + '\\\\')
    print('        Transformer Losses [kWh] & ' + str("{:.1f}".format(E_TrafoLosses_LFT)) + ' & ' + str("{:.1f}".format(E_TrafoLosses_PET)) + ' & ' + str("{:.1f}".format(E_TrafoLosses_HPET)) + '\\\\')
    print('        Total System Losses [kWh] & ' + str("{:.1f}".format(E_SystemLosses_LFT)) + ' & ' + str("{:.1f}".format(E_SystemLosses_PET)) + ' & ' + str("{:.1f}".format(E_SystemLosses_HPET)) + '\\\\\hline')
    print('    \end{tabular}')
    print('    \caption{Resulting computations of energy and losses in the power flow simulation.}')
    print('    \label{tab:Losses}')
    print('\end{table}')
    
    

    # PET Efficiency, Load [pu], PF
    if 'PET' in models:
        fig, ax = plt.subplots(1, figsize=figure_size)
        fig.subplots_adjust(right=0.85)
        ax2 = ax.twinx() # instantiate a second axes that shares the same x-axis        
        line1, = ax.plot(timeArray/60.0,PETefficiency,color='purple',linestyle='solid',linewidth=1.2, label='Efficiency')
        line2, = ax2.plot(timeArray/60.0,PET_Spu,color='crimson',linestyle='solid',linewidth=1.2, label='Load [pu]')    
        ax.grid(True)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set_ylim(0,1)
        ax2.set_ylim(0,1)        
        ax.set(ylabel = 'Efficiency', xlabel= 'Time [h]')
        ax2.set(ylabel = 'Load [pu]')        
        ax.yaxis.label.set_color(line1.get_color())
        ax2.yaxis.label.set_color(line2.get_color())        
        tkw = dict(size=4, width=1.5)
        ax.tick_params(axis='y', colors=line1.get_color(), **tkw)
        ax2.tick_params(axis='y', colors=line2.get_color(), **tkw)
        ax.tick_params(axis='x', **tkw)
        fig.savefig(ResultsFolder + '/Eff_Load' + fileName_end + '_PET' + str(PowerRatingPET) + '.pdf')
    
    
    # HPET Efficiency, Load [pu], PF
    if 'HPET' in models:
        fig, ax = plt.subplots(1, figsize=figure_size)
        fig.subplots_adjust(right=0.85)
        ax2 = ax.twinx() # instantiate a second axes that shares the same x-axis        
        line1, = ax.plot(timeArray/60.0,PEMefficiency1,color='purple',linestyle='solid',linewidth=1.2, label='Efficiency 1')
        line1, = ax.plot(timeArray/60.0,PEMefficiency2,color='purple',linestyle='dashed',linewidth=1.2, label='Efficiency 2')
        line2, = ax2.plot(timeArray/60.0,PEM_Sc1_pu,color='crimson',linestyle='solid',linewidth=1.2, label='Load 1 [pu]')
        line2, = ax2.plot(timeArray/60.0,PEM_Sc2_pu,color='crimson',linestyle='dashed',linewidth=1.2, label='Load 2 [pu]')        
        ax.grid(True)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set_ylim(0,1)
        ax2.set_ylim(0,1)        
        ax.set(ylabel = 'Efficiency', xlabel= 'Time [h]')
        ax2.set(ylabel = 'Load [pu]')
        ax.yaxis.label.set_color(line1.get_color())
        ax2.yaxis.label.set_color(line2.get_color())        
        tkw = dict(size=4, width=1.5)
        ax.tick_params(axis='y', colors=line1.get_color(), **tkw)
        ax2.tick_params(axis='y', colors=line2.get_color(), **tkw)
        ax.tick_params(axis='x', **tkw)
        fig.savefig(ResultsFolder + '/Eff_Load' + fileName_end + '_HPET' + str(FractionPowerHPET) + '.pdf')
        
        # Module 1 active and reactive power per unit
        fig, ax = plt.subplots(1, figsize=figure_size)
        for phase in range(3):
            line1, = ax.plot(timeArray/60.0,PEM_Pc1_pu[:,phase],color='red',linestyle='solid',linewidth=1.2)
            line2, = ax.plot(timeArray/60.0,PEM_Qc1_pu[:,phase],color='tab:blue',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'Power [pu]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend((r'$P_{C1}$',r'$Q_{C1}$'))
        fig.savefig(ResultsFolder + '/HPET_Pc1_Qc1_pu' + fileName_end + versus + '.pdf')
        
        # Module 2 active and reactive power per unit
        fig, ax = plt.subplots(1, figsize=figure_size)
        for phase in range(3):
            line1, = ax.plot(timeArray/60.0,PEM_Pc2_pu[:,phase],color='red',linestyle='solid',linewidth=1.2)
            line2, = ax.plot(timeArray/60.0,PEM_Qc2_pu[:,phase],color='tab:blue',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'Power [pu]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend((r'$P_{C2}$' , r'$Q_{C2}$'))
        fig.savefig(ResultsFolder + '/HPET_Pc2_Qc2_pu' + fileName_end + versus + '.pdf')



    # Total Losses [kW]
    fig, ax = plt.subplots(1, figsize=figure_size)
    if 'LFT' in models:
        line1, = ax.plot(timeArray/60.0,TotalLosses_LFT[:,0],color='forestgreen',linestyle='solid',linewidth=1.2)
    if 'PET' in models:
        line2, = ax.plot(timeArray/60.0,TotalLosses_PET[:,0],color='tab:blue',linestyle='solid',linewidth=1.5)
    if 'HPET' in models:
        line3, = ax.plot(timeArray/60.0,TotalLosses_HPET[:,0],color='red',linestyle='solid',linewidth=1.2)
    ax.set_xlim(0,24)
    plt.xticks(range(0,24+2,2),range(0,24+2,2))
    ax.set(ylabel = 'System Losses [kW]', xlabel= 'Time [h]')
    plt.grid(True)
    ax.legend(models)
    fig.savefig(ResultsFolder + '/Losses_kW' + fileName_end + versus + '.pdf')
    

    # Losses (reactive power)
    fig, ax = plt.subplots(1, figsize=figure_size)
    if 'LFT' in models:
        line1, = ax.plot(timeArray/60.0,TotalLosses_LFT[:,1],color='forestgreen',linestyle='solid',linewidth=1.2)
    if 'PET' in models:
        line2, = ax.plot(timeArray/60.0,TotalLosses_PET[:,1],color='tab:blue',linestyle='solid',linewidth=1.5)
    if 'HPET' in models:
        line2, = ax.plot(timeArray/60.0,TotalLosses_HPET[:,1],color='red',linestyle='solid',linewidth=1.2)
    ax.set_xlim(0,24)
    plt.xticks(range(0,24+2,2),range(0,24+2,2))
    ax.set(ylabel = 'System Reactive Losses [kVAr]', xlabel= 'Time [h]')
    plt.grid(True)
    ax.legend(models)
    fig.savefig(ResultsFolder + '/Losses_kVAr' + fileName_end + versus + '.pdf')
    
    
    # LFT, PET and HPET Losses [kW]
    fig, ax = plt.subplots(1, figsize=figure_size)
    if 'LFT' in models:
        line1, = ax.plot(timeArray/60.0,LFTlosses,color='forestgreen',linestyle='solid',linewidth=1.2)
    if 'PET' in models:
        line2, = ax.plot(timeArray/60.0,PETlosses,color='tab:blue',linestyle='solid',linewidth=1.5)
    if 'HPET' in models:
        line3, = ax.plot(timeArray/60.0,HPETlosses[:,0]+HPETlosses[:,1],color='red',linestyle='solid',linewidth=1.2)
    ax.set_xlim(0,24)
    plt.xticks(range(0,24+2,2),range(0,24+2,2))
    ax.set(ylabel = 'Transformer Losses [kW]', xlabel= 'Time [h]')
    plt.grid(True)
    ax.legend(models)
    fig.savefig(ResultsFolder + '/Device_Losses_kW' + fileName_end + versus + '.pdf')



    # Active Power at the secondary side (sum of the three phases)
    fig, ax = plt.subplots(1, figsize=figure_size)
    if 'LFT' in models:
        line1, = ax.plot(timeArray/60.0,P_LV_LFT,color='forestgreen',linestyle='solid',linewidth=1.2)
    if 'PET' in models:
        line2, = ax.plot(timeArray/60.0,P_LV_PET,color='tab:blue',linestyle='solid',linewidth=1.5)
    if 'HPET' in models:
        line3, = ax.plot(timeArray/60.0,P_LV_HPET,color='red',linestyle='solid',linewidth=1.2)
    ax.set_xlim(0,24)
    plt.xticks(range(0,24+2,2),range(0,24+2,2))
    ax.set(ylabel = 'P [kW]', xlabel= 'Time [h]')
    plt.grid(True)
    ax.legend(models)
    fig.savefig(ResultsFolder + '/P_TransformersOutput' + fileName_end + versus + '.pdf')
    

    # Reactive Power at the secondary side (sum of the three phases)
    fig, ax = plt.subplots(1, figsize=figure_size)
    if 'LFT' in models:
        line1, = ax.plot(timeArray/60.0,Q_LV_LFT,color='forestgreen',linestyle='solid',linewidth=1.2)
    if 'PET' in models:
        line2, = ax.plot(timeArray/60.0,Q_LV_PET,color='tab:blue',linestyle='solid',linewidth=1.5)
    if 'HPET' in models:
        line3, = ax.plot(timeArray/60.0,Q_LV_HPET,color='red',linestyle='solid',linewidth=1.2)
    ax.set_xlim(0,24)
    plt.xticks(range(0,24+2,2),range(0,24+2,2))
    ax.set(ylabel = 'Q [kVAr]', xlabel= 'Time [h]')
    plt.grid(True)
    ax.legend(models)
    fig.savefig(ResultsFolder + '/Q_TransformersOutput' + fileName_end + versus + '.pdf')

    
    # UPSTREAM NETWORK
    for i_line in SelectedLinesUS:
        
        # Voltage
        bus = Lines_upstream_data.loc[i_line,'Receiving bus']
        fig, ax = plt.subplots(1, figsize=figure_size)
        Vnom = Bus_Vnom.loc[bus,'Vnom_pn']
        for phase in range(3):
            if 'LFT' in models:
                line1, = ax.plot(timeArray/60.0,Data_LFT[:,2,np.where(Lines_set==i_line)[0][0],phase]/Vnom,color='forestgreen',linestyle='solid',linewidth=1.2)
            if 'PET' in models:
                line2, = ax.plot(timeArray/60.0,Data_PET[:,2,np.where(Lines_set==i_line)[0][0],phase]/Vnom,color='tab:blue',linestyle='solid',linewidth=1.5)
            if 'HPET' in models:
                line3, = ax.plot(timeArray/60.0,Data_HPET[:,2,np.where(Lines_set==i_line)[0][0],phase]/Vnom,color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'V [pu]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Vmag_Bus_' + bus + fileName_end + versus + '.pdf')
        
        
        # Current
        fig, ax = plt.subplots(1, figsize=figure_size)
        for phase in range(3):
            if 'LFT' in models:
                line1, = ax.plot(timeArray/60.0,Data_LFT[:,4,np.where(Lines_set==i_line)[0][0],phase],color='forestgreen',linestyle='solid',linewidth=1.2)
            if 'PET' in models:
                line2, = ax.plot(timeArray/60.0,Data_PET[:,4,np.where(Lines_set==i_line)[0][0],phase],color='tab:blue',linestyle='solid',linewidth=1.5)
            if 'HPET' in models:
                line3, = ax.plot(timeArray/60.0,Data_HPET[:,4,np.where(Lines_set==i_line)[0][0],phase],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'I [kA]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Imag_Line_' + i_line + fileName_end + versus + '.pdf')
    

        # Active Power in the primary side (sum of the three phases)
        fig, ax = plt.subplots(1, figsize=figure_size)
        if 'LFT' in models:
            line1, = ax.plot(timeArray/60.0,Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],0] + Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],1] + Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],2],color='forestgreen',linestyle='solid',linewidth=1.2)
        if 'PET' in models:
            line2, = ax.plot(timeArray/60.0,Data_PET[:,5,np.where(Lines_set==i_line)[0][0],0] + Data_PET[:,5,np.where(Lines_set==i_line)[0][0],1] + Data_PET[:,5,np.where(Lines_set==i_line)[0][0],2],color='tab:blue',linestyle='solid',linewidth=1.5)
        if 'HPET' in models:
            line3, = ax.plot(timeArray/60.0,Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],0] + Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],1] + Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],2],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'P [kW]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/P_Line_' + i_line + fileName_end + versus + '.pdf')
        

        # Active Power in the primary side (per phase)
        fig, ax = plt.subplots(1, figsize=figure_size)
        for phase in range(3):
            if 'LFT' in models:
                line1, = ax.plot(timeArray/60.0,Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],phase],color='forestgreen',linestyle='solid',linewidth=1.2)
            if 'PET' in models:
                line2, = ax.plot(timeArray/60.0,Data_PET[:,5,np.where(Lines_set==i_line)[0][0],phase],color='tab:blue',linestyle='solid',linewidth=1.5)
            if 'HPET' in models:
                line3, = ax.plot(timeArray/60.0,Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],phase],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'P [kW]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/P_Line_' + i_line + '_PerPhase' + fileName_end + versus + '.pdf')
    

        # Reactive Power (sum of the three phases)
        fig, ax = plt.subplots(1, figsize=figure_size)
        if 'LFT' in models:
            line1, = ax.plot(timeArray/60.0,Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],0] + Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],1] + Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],2],color='forestgreen',linestyle='solid',linewidth=1.2)
        if 'PET' in models:
            line2, = ax.plot(timeArray/60.0,Data_PET[:,6,np.where(Lines_set==i_line)[0][0],0] + Data_PET[:,6,np.where(Lines_set==i_line)[0][0],1] + Data_PET[:,6,np.where(Lines_set==i_line)[0][0],2],color='tab:blue',linestyle='solid',linewidth=2.0)
        if 'HPET' in models:
            line3, = ax.plot(timeArray/60.0,Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],0] + Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],1] + Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],2],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'Q [kVAr]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Q_Line_' + i_line + fileName_end + versus + '.pdf')
    
        # Reactive Power (per phase)
        fig, ax = plt.subplots(1, figsize=figure_size)
        for phase in range(3):
            if 'LFT' in models:
                line1, = ax.plot(timeArray/60.0,Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],phase],color='forestgreen',linestyle='solid',linewidth=1.2)
            if 'PET' in models:
                line2, = ax.plot(timeArray/60.0,Data_PET[:,6,np.where(Lines_set==i_line)[0][0],phase],color='tab:blue',linestyle='solid',linewidth=1.5)
            if 'HPET' in models:
                line3, = ax.plot(timeArray/60.0,Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],phase],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'Q [kVAr]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Q_Line_' + i_line + '_PerPhase' + fileName_end + versus + '.pdf')



    # DOWNSTREAM NETWORK
    
    # Transformer's secondary voltage
    bus = 'secondary'
    fig, ax = plt.subplots(1, figsize=figure_size)
    Vnom = Bus_Vnom.loc[bus,'Vnom_pn']
    for phase in range(3):
        if 'LFT' in models:
            line1, = ax.plot(timeArray/60.0,Data_LFT[:,0,1,phase]/Vnom,color='forestgreen',linestyle='solid',linewidth=1.2)
        if 'PET' in models:
            line2, = ax.plot(timeArray/60.0,Data_PET[:,0,1,phase]/Vnom,color='tab:blue',linestyle='solid',linewidth=1.5)
        if 'HPET' in models:
            line3, = ax.plot(timeArray/60.0,Data_HPET[:,0,1,phase]/Vnom,color='red',linestyle='solid',linewidth=1.2)
    ax.set_xlim(0,24)
    plt.xticks(range(0,24+2,2),range(0,24+2,2))
    ax.set(ylabel = 'V [pu]', xlabel= 'Time [h]')
    plt.grid(True)
    ax.legend(models)
    fig.savefig(ResultsFolder + '/Vmag_Bus_' + bus + fileName_end + versus + '.pdf')
    
    
    for i_line in SelectedLinesDS:
            
        # Voltage
        bus = Lines_downstream_data.loc[i_line,'Receiving bus']
        fig, ax = plt.subplots(1, figsize=figure_size)
        Vnom = Bus_Vnom.loc[bus,'Vnom_pn']
        for phase in range(3):
            if 'LFT' in models:
                line1, = ax.plot(timeArray/60.0,Data_LFT[:,2,np.where(Lines_set==i_line)[0][0],phase]/Vnom,color='forestgreen',linestyle='solid',linewidth=1.2)
            if 'PET' in models:
                line2, = ax.plot(timeArray/60.0,Data_PET[:,2,np.where(Lines_set==i_line)[0][0],phase]/Vnom,color='tab:blue',linestyle='solid',linewidth=1.5)
            if 'HPET' in models:
                line3, = ax.plot(timeArray/60.0,Data_HPET[:,2,np.where(Lines_set==i_line)[0][0],phase]/Vnom,color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'V [pu]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Vmag_Bus_' + bus + fileName_end + versus + '.pdf')
        

        # Current
        fig, ax = plt.subplots(1, figsize=figure_size)
        for phase in range(3):
            if 'LFT' in models:
                line1, = ax.plot(timeArray/60.0,Data_LFT[:,4,np.where(Lines_set==i_line)[0][0],phase],color='forestgreen',linestyle='solid',linewidth=1.2)
            if 'PET' in models:
                line2, = ax.plot(timeArray/60.0,Data_PET[:,4,np.where(Lines_set==i_line)[0][0],phase],color='tab:blue',linestyle='solid',linewidth=1.5)
            if 'HPET' in models:
                line3, = ax.plot(timeArray/60.0,Data_HPET[:,4,np.where(Lines_set==i_line)[0][0],phase],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'I [kA]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Imag_Line_' + i_line + fileName_end + versus + '.pdf')
    

        # Active Power (sum of the three phases)
        fig, ax = plt.subplots(1, figsize=figure_size)
        if 'LFT' in models:
            line1, = ax.plot(timeArray/60.0,Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],0] + Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],1] + Data_LFT[:,5,np.where(Lines_set==i_line)[0][0],2],color='forestgreen',linestyle='solid',linewidth=1.2)
        if 'PET' in models:
            line2, = ax.plot(timeArray/60.0,Data_PET[:,5,np.where(Lines_set==i_line)[0][0],0] + Data_PET[:,5,np.where(Lines_set==i_line)[0][0],1] + Data_PET[:,5,np.where(Lines_set==i_line)[0][0],2],color='tab:blue',linestyle='solid',linewidth=1.5)
        if 'HPET' in models:
            line3, = ax.plot(timeArray/60.0,Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],0] + Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],1] + Data_HPET[:,5,np.where(Lines_set==i_line)[0][0],2],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'P [kW]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/P_Line_' + i_line + fileName_end + versus + '.pdf')
    

        # Reactive Power (sum of the three phases)
        fig, ax = plt.subplots(1, figsize=figure_size)
        if 'LFT' in models:
            line1, = ax.plot(timeArray/60.0,Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],0] + Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],1] + Data_LFT[:,6,np.where(Lines_set==i_line)[0][0],2],color='forestgreen',linestyle='solid',linewidth=1.2)
        if 'PET' in models:
            line2, = ax.plot(timeArray/60.0,Data_PET[:,6,np.where(Lines_set==i_line)[0][0],0] + Data_PET[:,6,np.where(Lines_set==i_line)[0][0],1] + Data_PET[:,6,np.where(Lines_set==i_line)[0][0],2],color='tab:blue',linestyle='solid',linewidth=1.5)
        if 'HPET' in models:
            line3, = ax.plot(timeArray/60.0,Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],0] + Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],1] + Data_HPET[:,6,np.where(Lines_set==i_line)[0][0],2],color='red',linestyle='solid',linewidth=1.2)
        ax.set_xlim(0,24)
        plt.xticks(range(0,24+2,2),range(0,24+2,2))
        ax.set(ylabel = 'Q [kVAr]', xlabel= 'Time [h]')
        plt.grid(True)
        ax.legend(models)
        fig.savefig(ResultsFolder + '/Q_Line_' + i_line + fileName_end + versus + '.pdf')

    plt.close('all')
