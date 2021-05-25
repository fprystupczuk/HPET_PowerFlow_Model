import os
import platform
import math
import numpy as np
from scipy import interpolate
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter


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

# OpenDSS objects
DSSText = DSSObj.Text # to excecute OpenDSS text commands
DSSCircuit = DSSObj.ActiveCircuit # Use it to access the elements in the circuit (e.g., capacitors, buses, etc.)
DSSElement = DSSCircuit.ActiveCktElement
DSSSolution = DSSCircuit.Solution
DSSMonitors = DSSCircuit.Monitors
DSSBus = DSSCircuit.ActiveBus



# Whole HT variables
V_LV = np.zeros((3,2)) # voltages,angles for each phase of the transformer's LV terminal
V_LV_cmplx = np.zeros(3,dtype=np.complex_) # V_LV expressed as a complex in rectangular form
V_LV_nom_cmplx = np.zeros(3,dtype=np.complex_) # Nominal values of V_LV expressed as a complex in rectangular form
Vwin_nom = np.zeros(3)
P_LV = np.zeros(3)
Q_LV = np.zeros(3)

# Converter variables
Vc2_pu = np.zeros((3,2)) # Delta voltages,angles in pu for each phase of the STvsources
Vc2_cmplx = np.zeros(3,dtype=np.complex_) # V_LV expressed as a complex in rectangular form
Pc2 = np.zeros(3)
Qc2 = np.zeros(3)
Sc2 = np.zeros(3)
Qavail = np.zeros(3)
Pc1 = np.zeros(3)
Qc1 = np.zeros(3)

# Efficiency curves
load_points = np.array([0.025714285714286, 0.041428571428572, 0.057142857142857, 0.072857142857143, 0.088571428571429, 0.1, 0.112857142857143, 0.137142857142857, 0.161428571428572, 0.175714285714286, 0.192857142857143, 0.217142857142857, 0.25, 0.325714285714286, 0.4, 0.475714285714286, 0.55, 0.624285714285714, 0.7, 0.774285714285714, 0.85, 0.924285714285714, 1])
PF_points = np.array([0.7, 0.86, 1.0])
eff_points = np.zeros((PF_points.size,load_points.size))
eff_points[0,:] = [0.257971014492754, 0.359420289855072, 0.434782608695652, 0.492753623188406, 0.539130434782609, 0.568115942028985, 0.594202898550725, 0.634782608695652, 0.668115942028985, 0.684057971014493, 0.701449275362319, 0.723188405797101, 0.746376811594203, 0.784057971014493, 0.808695652173913, 0.827536231884058, 0.840579710144928, 0.850724637681159, 0.859420289855072, 0.865217391304348, 0.871014492753623, 0.87536231884058, 0.879710144927536]
eff_points[1,:] = [0.298550724637681, 0.407246376811594, 0.484057971014493, 0.543478260869565, 0.58695652173913, 0.615942028985507, 0.640579710144928, 0.679710144927536, 0.710144927536232, 0.72463768115942, 0.740579710144928, 0.759420289855072, 0.779710144927536, 0.81304347826087, 0.834782608695652, 0.850724637681159, 0.86231884057971, 0.871014492753623, 0.876811594202898, 0.882608695652174, 0.88695652173913, 0.889855072463768, 0.892753623188406]
eff_points[2,:] = [0.330434782608696, 0.443478260869565, 0.521739130434783, 0.578260869565217, 0.623188405797102, 0.649275362318841, 0.672463768115942, 0.710144927536232, 0.73768115942029, 0.752173913043478, 0.766666666666667, 0.784057971014493, 0.802898550724638, 0.831884057971014, 0.852173913043478, 0.865217391304348, 0.87536231884058, 0.882608695652174, 0.888405797101449, 0.892753623188406, 0.897101449275362, 0.9, 0.901449275362319]
eff_points = (0.97/eff_points[2,-1]) * eff_points # offset efficiency to a different maximum
efficiency = np.zeros(3)


# Load OpenDSS circuit and solve it
DSSText.Command = 'set datapath=' + os.path.dirname(os.path.realpath(__file__))
DSSText.Command = 'Redirect ./test_HTlines_QLinearSweep_Delta.txt'

DSSText.Command = 'set voltageBases=[' + str(10.0) + ',' + str(0.4) + ',' + str(0.04) + ']'
DSSText.Command = 'CalcVoltageBase'

DSSText.Command = 'set controlmode=static'
DSSText.Command = 'set mode=snapshot'

DSSText.Command = 'Reset Monitors'

DSSSolution.Solve()
if not(DSSSolution.Converged):
    raise ValueError('Solution did not Converge')
DSSMonitors.SampleAll()
DSSMonitors.SaveAll()  


# Get windings' power rating and nominal voltages
DSSCircuit.Transformers.Name = 'Trafo1'
DSSCircuit.Transformers.Wdg = 1
S_max = DSSCircuit.Transformers.kVA
Vwin_nom[0] = DSSCircuit.Transformers.kV*1000
DSSCircuit.Transformers.Wdg = 2
Vwin_nom[1] = DSSCircuit.Transformers.kV*1000
DSSCircuit.Transformers.Wdg = 3
Vwin_nom[2] = DSSCircuit.Transformers.kV*1000



# CALCULATIONS (all quantities in p.u.)
P2_points = 200
PF2_min = 0.8
Q2 = np.linspace(0,S_max*math.sin(math.acos(PF2_min)),P2_points) # Linear sweep of Q2
P2 = np.sqrt(S_max**2 - Q2**2)
PF2 = P2/S_max

FractionPowerHT = np.array([0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 1e-6]) # FractionPowerHT = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
Qc_max = (S_max*FractionPowerHT)/(1+FractionPowerHT) # Is this causing that the lines in the graph Q1 vs. Q2 are not equally separated?
# S1 = S2 + S3
# S3 = f * S2
# S1 = S2 + f*S2 = S2*(1+f)
# S2 = S1/(1+f)
# S3 = f*S1/(1+f)

Q1 = np.zeros([P2_points,FractionPowerHT.size])
PF1 = np.zeros([P2_points,FractionPowerHT.size])

# Ideal reactive power compensation
for j in range(FractionPowerHT.size):
    for i in range(P2_points):
        if(Q2[i] > Qc_max[j]):
            Q1[i,j] = Q2[i] - Qc_max[j] # Reactive power seen by the secondary winding and transmitted (ideally) to the primary winding
        else:
            Q1[i,j] = 0
    PF1[:,j] = np.cos(np.arctan(Q1[:,j]/P2)) # P2 is here because we're supposing that P1 = P2

# Variables to store results
P_prim = np.zeros((PF2.size,FractionPowerHT.size))
Q_prim = np.zeros((PF2.size,FractionPowerHT.size))
P_sec = np.zeros((PF2.size,FractionPowerHT.size))
Q_sec = np.zeros((PF2.size,FractionPowerHT.size))
P_aux = np.zeros((PF2.size,FractionPowerHT.size))
Q_aux = np.zeros((PF2.size,FractionPowerHT.size))
PF1_model = np.zeros((PF2.size,FractionPowerHT.size))


for j in range(FractionPowerHT.size):

    # Redefine secondary and auxiliary power ratings
    Smax_secondary = S_max/(FractionPowerHT[j]+1)
    Smax_auxiliary = (S_max*FractionPowerHT[j])/(FractionPowerHT[j]+1)    
    for phase in range(3):        
        DSSCircuit.Transformers.Wdg = 1
        Vwin_nom_primary = DSSCircuit.Transformers.kV
        DSSCircuit.Transformers.Wdg = 2
        Vwin_nom_secondary = DSSCircuit.Transformers.kV
        Vwin_nom_auxiliary = Vwin_nom_secondary * FractionPowerHT[j] * math.sqrt(3)
        
        DSSText.Command = 'Edit Transformer.Trafo' + str(phase+1) + ' Wdg=3 kV=' + str(Vwin_nom_auxiliary) + ' kVA=' + str(Smax_auxiliary)
        DSSText.Command = 'Edit Vsource.STvsource' + str(phase+1) + ' BasekV=' + str(Vwin_nom_auxiliary/math.sqrt(3))
        DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kV=' + str(Vwin_nom_auxiliary)

    # Redefine voltage bases
    DSSText.Command ='set VoltageBases=[' + str(Vwin_nom_primary) + ',' + str(Vwin_nom_secondary * math.sqrt(3)) + ',' + str(Vwin_nom_auxiliary) + ']'
    DSSText.Command = 'CalcVoltageBase'
    DSSText.Command = 'set controlmode=static'
    DSSText.Command = 'set mode=snapshot'


    # Create variable for STvsources nominal voltage (Base kV)
    Vc2_nom = Vwin_nom[1]*FractionPowerHT[j]
    Sc2_nom = (S_max*FractionPowerHT[j])/(1+FractionPowerHT[j]) # Per-phase PEM's nominal apparent power in kVA

    
    for i in range(P2_points):
        
        # Assign load in function of the load sweep through all the values of Q2
        for phase in range(3):
            DSSCircuit.Loads.Name = 'Feeder1_LOAD' + str(phase+1)
            DSSCircuit.Loads.kW = P2[i]
            DSSCircuit.Loads.kvar = Q2[i]
            
#        # Rseset STvsources
#        for phase in range(3):
#            DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
#            DSSCircuit.Vsources.pu = 0.0
            
        # SOLVE: get actual voltages and calculate voltage regulation
        DSSText.Command ='Reset Monitors'
        DSSSolution.Solve()
        if not(DSSSolution.Converged):
            raise ValueError('Solution did not Converge')
        DSSMonitors.SampleAll()
        DSSMonitors.SaveAll()        
        

        ## Read secondary winding actual voltages and angles per phase
        DSSMonitors.Name = 'Trafo_sec-VI'
        for phase in range(3):
            V_LV[phase,0] = DSSMonitors.Channel(phase * 2 + 1)[0]
            V_LV[phase,1] = DSSMonitors.Channel(phase * 2 + 2)[0]
            V_LV_cmplx[phase] = V_LV[phase,0] * np.exp(1j*math.radians(V_LV[phase,1]))
            V_LV_nom_cmplx[phase] = Vwin_nom[1] * np.exp(1j*math.radians(30 - 120.0*phase))
            DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
            Vc2_cmplx[phase] = DSSCircuit.Vsources.pu * Vc2_nom * np.exp(1j*math.radians(DSSCircuit.Vsources.AngleDeg))
        
        #Vc2_cmplx = V_LV_nom_cmplx - V_LV_cmplx # In this situation V_LV = Vt because the STvsources are passivated
        Vc2_cmplx = V_LV_nom_cmplx - V_LV_cmplx + Vc2_cmplx
        
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
            
        
        # SOLVE: get the actual Pc1 and Qc1 corresponding to the operating point previously set
        DSSText.Command ='Reset Monitors'
        DSSSolution.Solve()
        if not(DSSSolution.Converged):
            raise ValueError('Solution did not Converge')
        DSSMonitors.SampleAll()
        DSSMonitors.SaveAll()
        
#        DSSMonitors.Name = 'Trafo_sec-VI'
#        for phase in range(3):
#            print('V_LV' + str(phase+1) + ' = ' + str(DSSMonitors.Channel(phase * 2 + 1)[0]/Vwin_nom[1]) + '/_' + str(DSSMonitors.Channel(phase * 2 + 2)[0]))

        
        # Reactive power compensation
        
        ## Read P_LV and Q_LV at the secondary winding. Those powers are measured at Transformer's output, thus they are negative
        DSSMonitors.Name = 'Trafo_sec-PQ'
        for phase in range(3):
            P_LV[phase] = DSSMonitors.Channel(phase * 2 + 1)[0] # P_LV is not used
            Q_LV[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
        
        ## Read Pc2 and Qc2 at the converter's output per phase
        ## Those powers are measured at Vsource's output, thus they are negative
        DSSMonitors.Name = 'STvsources-PQ'
        for phase in range(3):
            Pc2[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            Qc2[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
        
        # Normalization of S, per phase
        Sc2 = np.sqrt(Pc2**2 + Qc2**2)
        Sc2_pu = Sc2/Sc2_nom
        PFvsource = abs(Pc2/Sc2)
        
        # Reactive power compensation capability
        Qt = Q_LV - Qc2
        for phase in range(3):
            if abs(Pc2[phase]) > Sc2_nom:
                Qavail[phase] = 0
            else:
                Qavail[phase] = math.sqrt(Sc2_nom**2 - Pc2[phase]**2)   
                if abs(Qt[phase]) > Qavail[phase]:
                    Qc1[phase] = -(abs(Qt[phase])/Qt[phase]) * Qavail[phase] # Take the actual sign of Q_LV[phase] and invert it for compensation
                    #print('Q compensation exceeded')
                else:
                    Qc1[phase] = -Qt[phase] # Total compensation: Qp inyects directly to the transformer's core


        ## Efficiency calculation
        
        # Efficiency curve for linear interpolation: efficiency = f(S,PF)
        eff_interpFunction = interpolate.interp2d(load_points, PF_points, eff_points, kind='linear')
        
        # Calculate efficiency and active power in the primary side, per phase
        for phase in range(3):
            efficiency[phase] = eff_interpFunction(Sc2_pu[phase], PFvsource[phase]) # interpolate actual efficiency
            if Pc2[phase] >= 0: # forward power flow
                Pc1[phase] = Pc2[phase] / efficiency[phase]
            else: # reverse power flow
                Pc1[phase] = Pc2[phase] * efficiency[phase]
    
        # Update STloads values
        for phase in range(3):
#            DSSCircuit.Loads.Name = 'STload' + str(phase+1)
#            DSSCircuit.Loads.kW = Pc1[phase]
#            DSSCircuit.Loads.kvar = Qc1[phase] # In this case Qc1 is not cero, to compensate primary side Q
            DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kW=' + str(Pc1[phase]) + ' kvar=' + str(Qc1[phase])
        
        
        # SOLVE: final result for the current time step
        DSSText.Command ='Reset Monitors'
        DSSSolution.Solve()
        if not(DSSSolution.Converged):
            raise ValueError('Solution did not Converge')        
        DSSMonitors.SampleAll() 
        DSSMonitors.SaveAll()
        
    
        # Read only the first phase since the system is balanced
        DSSMonitors.Name = 'Trafo_prim-PQ'
        P_prim[i,j] = DSSMonitors.Channel(1)[0]
        Q_prim[i,j] = DSSMonitors.Channel(2)[0]
        
        DSSMonitors.Name = 'Trafo_sec-PQ'
        P_sec[i,j] = DSSMonitors.Channel(1)[0]
        Q_sec[i,j] = DSSMonitors.Channel(2)[0]
        
        DSSMonitors.Name = 'Trafo_aux-PQ'
        P_aux[i,j] = DSSMonitors.Channel(1)[0]
        Q_aux[i,j] = DSSMonitors.Channel(2)[0]
    
        PF1_model[i,j] = np.cos(np.arctan(Q_prim[i,j]/P_prim[i,j]))



# FIGURES
figure_size = (10, 9)
font_size = 20
plt.rcParams.update({'font.size': font_size, 'mathtext.fontset':'cm', 'mathtext.rm': 'serif'})

sns.set_palette("tab10",12) # See colormaps on https://matplotlib.org/tutorials/colors/colormaps.html

# Q1 vs. P1
fig, ax = plt.subplots(1, figsize=figure_size)
for j in range(FractionPowerHT.size-1):
    ax.plot(P_prim[:,j]/S_max,Q_prim[:,j]/S_max, label=r'$\alpha$ = ' + str(FractionPowerHT[j]))
ax.plot(P_prim[:,-1]/S_max,Q_prim[:,-1]/S_max, '--', color='silver', label='No PEC')
plt.grid(True)
ax.set(xlabel = r'$P_{MV}\mathrm{\ [pu]}$', ylabel= r'$Q_{MV}\mathrm{\ [pu]}$')
ax.legend()
#plt.title(r'$Q_1$ vs. $Q_2$ @ Noloadloss=0.08125, Loadloss=0.875')
fig.savefig('./test_HTlines_QLinearSweep_Graph/Q1vsP1_Noloadloss0.08125_Loadloss0.875.pdf')


sns.set_palette("tab20",12) # See colormaps on https://matplotlib.org/tutorials/colors/colormaps.html

# Q1 vs. Q2
fig, ax = plt.subplots(1, figsize=figure_size)
for j in range(FractionPowerHT.size-1):
    ax.plot(Q_sec[:,j]/S_max,Q_prim[:,j]/S_max, label=r'$\alpha$ = ' + str(FractionPowerHT[j]))
    ax.plot(Q2/S_max,Q1[:,j]/S_max, '--')
ax.plot(Q_sec[:,-1]/S_max,Q_prim[:,-1]/S_max, '--', color='silver', label='No PEC')
plt.grid(True)
ax.set(xlabel = r'$Q_{LV}\mathrm{\ [pu]}$', ylabel= r'$Q_{MV}\mathrm{\ [pu]}$')
ax.legend()
#plt.title(r'$Q_1$ vs. $Q_2$ @ Noloadloss=0.08125, Loadloss=0.875')
fig.savefig('./test_HTlines_QLinearSweep_Graph/Q1vsQ2_Noloadloss0.08125_Loadloss0.875.pdf')
#plt.title(r'$Q_1$ vs. $Q_2$ @ Noloadloss=0, Loadloss=0')
#fig.savefig('./test_HTlines_QLinearSweep_Graph/Q1vsQ2_Noloadloss0.0_Loadloss0.0.pdf')


# PF1 vs. Q2
fig, ax = plt.subplots(1, figsize=figure_size)
for j in range(FractionPowerHT.size-1):
    ax.plot(Q2/S_max,PF1_model[:,j], label=r'$\alpha$ = ' + str(FractionPowerHT[j]))
    ax.plot(Q2/S_max,PF1[:,j], '--')
ax.plot(Q2/S_max,PF1_model[:,-1], '--', color='silver', label='No PEC')
plt.grid(True)
ax.set(xlabel = r'$Q_{LV}\mathrm{\ [pu]}$', ylabel= r'$PF_{MV}$')
#ax.set(xlabel = r'Test $\alpha$', ylabel= 'PF1')
ax.legend()
#plt.title(r'$PF_1$ vs. $Q_2$ @ Noloadloss=0.08125, Loadloss=0.875')
fig.savefig('./test_HTlines_QLinearSweep_Graph/PF1vsQ2_Noloadloss0.08125_Loadloss0.875.pdf')
#plt.title(r'$PF_1$ vs. $Q_2$ @ Noloadloss=0, Loadloss=0')
#fig.savefig('./test_HTlines_QLinearSweep_Graph/PF1vsQ2_Noloadloss0.0_Loadloss0.0.pdf')


# PF1 vs. PF2
fig, ax = plt.subplots(1, figsize=figure_size)
for j in range(FractionPowerHT.size-1):
    ax.plot(PF2,PF1_model[:,j], label=r'$\alpha$ = ' + str(FractionPowerHT[j]))
    ax.plot(PF2,PF1[:,j], '--')
ax.plot(PF2,PF2, '--', color='silver', label='No PEC')
ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
plt.grid(True)
ax.set(xlabel = r'$PF_{LV}$', ylabel= r'$PF_{MV}$')
ax.legend()
#plt.title(r'$PF_1$ vs. $PF_2$ @ Noloadloss=0.08125, Loadloss=0.875')
fig.savefig('./test_HTlines_QLinearSweep_Graph/PF1vsPF2_Noloadloss0.08125_Loadloss0.875.pdf')
#plt.title(r'$PF_1$ vs. $PF_2$ @ Noloadloss=0, Loadloss=0')
#fig.savefig('./test_HTlines_QLinearSweep_Graph/PF1vsPF2_Noloadloss0.0_Loadloss0.0.pdf')

plt.close('all')
