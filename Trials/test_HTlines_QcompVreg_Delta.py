import platform
import math
import numpy as np
from scipy import interpolate

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


# Whole HT variables
V_LV = np.zeros((3,2)) # voltages,angles for each phase of the transformer's LV terminal
V_LV_cmplx = np.zeros(3,dtype=np.complex_) # V_LV expressed as a complex in rectangular form
V_LV_nom_cmplx = np.zeros(3,dtype=np.complex_) # Nominal values of V_LV expressed as a complex in rectangular form
Vwin_nom = np.zeros(3)
P_LV = np.zeros(3)
Q_LV = np.zeros(3)

# Converter variables
Vc2_pu = np.zeros((3,2)) # Delta voltages,angles in pu for each phase of the STvsources
Vc2_cmplx = np.zeros(3) # Vvsources expressed as a complex in rectangular form
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
DSSText.Command = 'Redirect ./test_HTlines_QcompVreg_Delta.txt'

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

FractionPowerHT = 0.1

# Get windings' power rating and nominal voltages
DSSCircuit.Transformers.Name = 'Trafo1'
DSSCircuit.Transformers.Wdg = 1
Vwin_nom_primary = DSSCircuit.Transformers.kV
Vwin_nom[0] = DSSCircuit.Transformers.kV*1000
S_max = DSSCircuit.Transformers.kVA
Smax_auxiliary = (S_max*FractionPowerHT)/(FractionPowerHT+1)
DSSCircuit.Transformers.Wdg = 2
Vwin_nom[1] = DSSCircuit.Transformers.kV*1000
Vwin_nom_secondary = DSSCircuit.Transformers.kV
Vwin_nom_auxiliary = Vwin_nom_secondary * FractionPowerHT * math.sqrt(3)
DSSCircuit.Transformers.Wdg = 3
Vwin_nom[2] = DSSCircuit.Transformers.kV*1000

# Redefine voltages and power ratings
for phase in range(3):
    DSSText.Command = 'Edit Transformer.Trafo' + str(phase+1) + ' Wdg=3 kV=' + str(Vwin_nom_auxiliary) + ' kVA=' + str(Smax_auxiliary)
    DSSText.Command = 'Edit Vsource.STvsource' + str(phase+1) + ' BasekV=' + str(Vwin_nom_auxiliary/math.sqrt(3))
    DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kV=' + str(Vwin_nom_auxiliary)

# Redefine voltage bases
DSSText.Command ='set VoltageBases=[' + str(Vwin_nom_primary) + ',' + str(Vwin_nom_secondary * math.sqrt(3)) + ',' + str(Vwin_nom_auxiliary) + ']'
DSSText.Command = 'CalcVoltageBase'
DSSText.Command = 'set controlmode=static'
DSSText.Command = 'set mode=snapshot'


# Create variable for STvsources nominal voltage (Base kV)
Vvsources_nom = Vwin_nom[1]*FractionPowerHT
Sc2_nom = (S_max*FractionPowerHT)/(1+FractionPowerHT) # Per-phase PEM's nominal apparent power in kVA


# Rseset STvsources
for phase in range(3):
    DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
    DSSCircuit.Vsources.pu = 0.0
    
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

Vc2_cmplx = V_LV_nom_cmplx - V_LV_cmplx # In this situation V_LV = Vt because the STvsources are passivated

Vc2_pu[:,0] = np.absolute(Vc2_cmplx)/Vvsources_nom
Vc2_pu[:,1] = np.degrees(np.angle(Vc2_cmplx))
print('V_LV: ' + str(V_LV[:,0]/Vwin_nom[1]) + '/_' + str(V_LV[:,1]))
print('V_vs: ' + str(Vc2_pu[:,0]) + '/_' + str(Vc2_pu[:,1]))

# Update STvsources pu and agnle parameters
for phase in range(3):
    DSSCircuit.Vsources.Name = 'STvsource' + str(phase+1)
    if abs(Vc2_pu[phase,0]) > 1.0: # Voltage limitation of the BtB converter
        DSSCircuit.Vsources.pu = 1.0
        print('Voltage regulation exceeded ' + str(Vc2_pu[phase,:]))
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

DSSMonitors.Name = 'Trafo_sec-VI'
for phase in range(3):
    V_LV[phase,0] = DSSMonitors.Channel(phase * 2 + 1)[0]
    V_LV[phase,1] = DSSMonitors.Channel(phase * 2 + 2)[0]
print('V_LV: ' + str(V_LV[:,0]/Vwin_nom[1]) + '/_' + str(V_LV[:,1]))


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

## Reactive power compensation capability
Qt = Q_LV - Qc2
for phase in range(3):
    if abs(Pc2[phase]) > Sc2_nom:
        Qavail[phase] = 0
    else:
        Qavail[phase] = math.sqrt(Sc2_nom**2 - Pc2[phase]**2)   
        if abs(Qt[phase]) > Qavail[phase]:
            Qc1[phase] = -(abs(Qt[phase])/Qt[phase]) * Qavail[phase] # Take the actual sign of Q_LV[phase] and invert it for compensation
            print('Q compensation exceeded')
        else:
            Qc1[phase] = -Qt[phase] # Total compensation: Qp inyects directly to the transformer's core

# Efficiency calculation
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
    DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kW=' + str(Pc1[phase]) + ' kvar=' + str(Qc1[phase])


# SOLVE: final result for the current time step
DSSText.Command ='Reset Monitors'
DSSSolution.Solve()
if not(DSSSolution.Converged):
    raise ValueError('Solution did not Converge')        
DSSMonitors.SampleAll() 
DSSMonitors.SaveAll()



print ('Powers:\t\tP [kW]\t\tQ [kVAr]\tS [kVA]')
for phase in range(3):
    print('\tPhase ' + str(phase+1))
    DSSMonitors.Name = 'Trafo_prim-PQ'
    print('Prim Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))
    DSSMonitors.Name = 'Trafo_sec-PQ'
    print('Sec Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))
    DSSMonitors.Name = 'Trafo_aux-PQ'
    print('Aux Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))) + '\n')

for phase in range(3):
    DSSMonitors.Name = 'Feeder1_LOAD' + str(phase+1) + '-PQ'
    print('\tLoad' + str(phase+1) + '\t' + str("{:.2f}".format(DSSMonitors.Channel(1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(1)[0]**2 + DSSMonitors.Channel(2)[0]**2))))

print()
DSSMonitors.Name = 'STvsources-PQ'
for phase in range(3):
    print('\tSTvs\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))
print()
DSSMonitors.Name = 'STloads-PQ'
for phase in range(3):
    print('\tSTload\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))

print('\nVoltages:')
for phase in range(3):
    print('\tPhase ' + str(phase+1))
    DSSMonitors.Name = 'Trafo_prim-VI'
    print('Prim Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    DSSMonitors.Name = 'Trafo_sec-VI'
    print('Sec Wdg (V)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('Sec Wdg (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg')
    DSSMonitors.Name = 'Trafo_aux-VI'
    print('Aux Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg' + '\n')

    DSSMonitors.Name = 'STvsources-VI'
    print('STvs' + str(phase+1) + ' (V)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('STvs' + str(phase+1) + ' (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg' + '\n')

    DSSMonitors.Name = 'STloads-VI'
    print('STload' + str(phase+1) + ' (V)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('STload' + str(phase+1) + ' (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg' + '\n')

for phase in range(3):
    DSSMonitors.Name = 'Feeder1_LOAD' + str(phase+1) + '-VI'
    print('\tLoad' + str(phase+1) + '\t' + str("{:.2f}".format(DSSMonitors.Channel(1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(2)[0])) + ' deg')
