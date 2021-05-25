import platform
import math
import numpy as np


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

DSSText.Command = 'Redirect ./test_HTlines_Delta.txt'

## Line codes
#DSSText.Command = 'Redirect ' + '/home/federico/MEGAsync/UCD/PhD/OpenDSS/OpenDSS_LFTvsSSTvsHT/NetworkData/NetworkModel/CableData/OpenDSS_line_types.txt'
## Cable type data
#CableData = pd.read_excel('/home/federico/MEGAsync/UCD/PhD/OpenDSS/OpenDSS_LFTvsSSTvsHT/NetworkData/NetworkModel/CableData/CableData.xlsx',index_col=0 )
#CableData.index = CableData.index.str.lower()
#Feeders = ['Feeder1','Feeder2','Feeder3'] # ['Feeder1','Feeder2','Feeder3'] - can use all three or 1 or 2
## Lines and loads
#for i_feeder in Feeders:
#    # Compile OpenDSS codes for lines, loads and monitors
#    DSSText.Command = 'Redirect ' + '/home/federico/MEGAsync/UCD/PhD/OpenDSS/OpenDSS_LFTvsSSTvsHT/NetworkData/NetworkModel/' + i_feeder + '/OpenDSS_Lines.txt'
#    DSSText.Command = 'Redirect ' + '/home/federico/MEGAsync/UCD/PhD/OpenDSS/OpenDSS_LFTvsSSTvsHT/NetworkData/NetworkModel/' + i_feeder + '/OpenDSS_LVloads_flat.txt'
#    DSSText.Command = 'Redirect ' + '/home/federico/MEGAsync/UCD/PhD/OpenDSS/OpenDSS_LFTvsSSTvsHT/NetworkData/NetworkModel/' + i_feeder + '/OpenDSS_monitors.txt'
#    DSSText.Command = 'Redirect ' + '/home/federico/MEGAsync/UCD/PhD/OpenDSS/OpenDSS_LFTvsSSTvsHT/NetworkData/NetworkModel/' + i_feeder + '/OpenDSS_LoadMonitors.txt'




# Extract transformer power rating
DSSCircuit.Transformers.Name = 'Trafo1'
DSSCircuit.Transformers.Wdg = 1
TransformerRating = DSSCircuit.Transformers.kVA
FractionPowerHT = 0.10
Smax_auxiliary = (TransformerRating*FractionPowerHT)/(FractionPowerHT+1)
print('PEM rating: ' + str("{:.2f}".format(Smax_auxiliary)))


# Redefine power ratings and nominal voltages for transformers, STvsources and STloads
for phase in range(3):
    DSSText.Command = 'Edit Transformer.Trafo' + str(phase+1) + ' Wdg=3 kV=' + str(0.4*FractionPowerHT) + ' kVA=' + str(Smax_auxiliary)
    DSSText.Command = 'Edit Vsource.STvsource' + str(phase+1) + ' BasekV=' + str(0.4*FractionPowerHT/math.sqrt(3))
    DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kV=' + str(0.4*FractionPowerHT)
    
# Redefine voltage bases
DSSText.Command = 'set voltageBases=[' + str(10.0) + ',' + str(0.4) + ',' + str(0.4*FractionPowerHT) + ']'
DSSText.Command = 'CalcVoltageBase'

DSSText.Command = 'set controlmode=static'
DSSText.Command = 'set mode=snapshot'

DSSText.Command = 'Reset Monitors'
DSSSolution.Solve()
if not(DSSSolution.Converged):
    raise ValueError('Solution did not Converge')
DSSMonitors.SampleAll()
DSSMonitors.SaveAll()  

V_LV = np.zeros(3)
I_LV = np.zeros(3)
P_LV = np.zeros(3)
Q_LV = np.zeros(3)
P_MV = np.zeros(3)
Q_MV = np.zeros(3)
V_Seq = np.zeros(3,dtype=np.complex_)
I_Seq = np.zeros(3,dtype=np.complex_)
Qneg = np.zeros(3)
Qpos = np.zeros(3)
Pc2 = np.zeros(3)
Qc2 = np.zeros(3)
Qavail = np.zeros(3)
Pc1 = np.zeros(3)
Qc1 = np.zeros(3)

DSSMonitors.Name = 'Trafo_prim-PQ'
print('Primary Q before compensation: ' + str("{:.2f}".format(DSSMonitors.Channel(0 * 2 + 2)[0])) + ' | '+ str("{:.2f}".format(DSSMonitors.Channel(1 * 2 + 2)[0])) + ' | '+ str("{:.2f}".format(DSSMonitors.Channel(2 * 2 + 2)[0])))

# Read P_LV and Q_LV at the secondary winding. Those powers are measured at Transformer's output, thus they are negative
DSSMonitors.Name = 'Trafo_sec-PQ'
for phase in range(3):
    P_LV[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
    Q_LV[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]

DSSMonitors.Name = 'Trafo_prim-PQ'
for phase in range(3):
    P_MV[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
    Q_MV[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]

DSSMonitors.Name = 'STvsources-PQ'
for phase in range(3):
    Pc2[phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
    Qc2[phase] = DSSMonitors.Channel(phase * 2 + 2)[0]



# Reactive power compensation for unbalanced voltage case
    
## Take measurements on primary side
#V_MVpn = np.zeros(3,dtype=np.complex_)
#I_MV_phase = np.zeros(3,dtype=np.complex_)
#for phase in range(3):
#    DSSMonitors.Name = 'Trafo' + str(phase+1) + '_prim-VI'
#    I_MV_phase[phase] = DSSMonitors.Channel(5)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(6)[0]))
#    V_MVpn[phase] = DSSMonitors.Channel(1)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(2)[0]))
#V_MVpp = np.array(((V_MVpn[0]-V_MVpn[1]),(V_MVpn[1]-V_MVpn[2]),(V_MVpn[2]-V_MVpn[0])),dtype=np.complex_)
#I_MV_line = np.zeros(3,dtype=np.complex_)
#DSSMonitors.Name = 'Trafo_prim-VI'
#for phase in range(3):
#    I_MV_line[phase] = DSSMonitors.Channel(phase * 2 + 7)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 8)[0]))
#

# Take measurements on secondary side
V_LVpn = np.zeros(3,dtype=np.complex_)
I_LV_line = np.zeros(3,dtype=np.complex_)
DSSMonitors.Name = 'Trafo_sec-VI'
for phase in range(3):
    I_LV_line[phase] = DSSMonitors.Channel(phase * 2 + 7)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 8)[0]))
    V_LVpn[phase] = DSSMonitors.Channel(phase * 2 + 1)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 2)[0]))
V_LVpp = np.array(((V_LVpn[0]-V_LVpn[1]),(V_LVpn[1]-V_LVpn[2]),(V_LVpn[2]-V_LVpn[0])),dtype=np.complex_)

## Numerical example from Blasco, P. et al. (Year 2020)
#V_LVpp = np.zeros(3,dtype=np.complex_)
#I_LV_line = np.zeros(3,dtype=np.complex_)
#V_LVpp = (405.15 * np.exp(1j*math.radians(38.5)), 353.49 * np.exp(-1j*math.radians(95.68)), 299.1 * np.exp(1j*math.radians(160.56)))
#I_LV_line = (14.88 * np.exp(-1j*math.radians(138.77)), 25.51 * np.exp(-1j*math.radians(84.49)), 36.27 * np.exp(1j*math.radians(76.05)))

# Calculate symetrical sequence components from phase voltages and currents (the components for phase 2 and 3 must be obtained by using the properties of the positive-, negative-, and zero-sequence sets which they represent (Bergen1999, pag 448))
a = np.exp(1j*math.radians(120.0))
b = math.sqrt(3) * np.exp(1j*math.radians(30.0))
Ainv = (1/3) * np.array(([1,1,1],[1,a,a**2],[1,a**2,a]),dtype=np.complex_)

V_Seq = np.dot(Ainv,V_LVpp)
I_Seq = np.dot(Ainv,I_LV_line)
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
    


## Compensation for balanced voltage case and using Symetrical Components read from OpenDSS in the secondary side
#DSSMonitors.Name = 'Trafo_sec-VI_s'
#for phase in range(3):
#    V_Seq[phase] = DSSMonitors.Channel(phase * 2 + 1)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 2)[0]))
#    I_Seq[phase] = DSSMonitors.Channel(phase * 2 + 7)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 8)[0]))        
#a = np.exp(1j*math.radians(120.0))
#V_s1 = np.array((V_Seq[1], a**2*V_Seq[1], a*V_Seq[1]),dtype=np.complex_)
##V_s2 = np.array((V_Seq[2], a*V_Seq[2], a**2*V_Seq[2]),dtype=np.complex_)
#I_s1 = np.array((I_Seq[1], a**2*I_Seq[1], a*I_Seq[1]),dtype=np.complex_)
#I_s2 = np.array((I_Seq[2], a*I_Seq[2], a**2*I_Seq[2]),dtype=np.complex_)        
#Qpos = np.abs(I_s1) * np.abs(V_s1) * np.sin(np.angle(V_s1) - np.angle(I_s1))
#Qneg = 2 * np.abs(I_s2) * np.abs(V_s1) * np.sin(np.angle(V_s1) - np.angle(I_s2))        
#Qt = Qpos/1000 + Qneg/1000 - Qc2


## Compensation for balanced voltage case and taking measurements from primary side
#V_MVpn = np.zeros(3,dtype=np.complex_)
#I_MVpp = np.zeros(3,dtype=np.complex_)
#for phase in range(3):
#    DSSMonitors.Name = 'Trafo' + str(phase+1) + '_prim-VI'
#    V_MVpn[phase] = DSSMonitors.Channel(1)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(2)[0]))
#    I_MVpp[phase] = DSSMonitors.Channel(5)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(6)[0]))
#
#V_MVpp = np.array(((V_MVpn[0]-V_MVpn[1]),(V_MVpn[1]-V_MVpn[2]),(V_MVpn[2]-V_MVpn[0])),dtype=np.complex_)
#
#I_MV_line = np.zeros(3,dtype=np.complex_)
#DSSMonitors.Name = 'Trafo_prim-VI'
#for phase in range(3):
#    I_MV_line[phase] = DSSMonitors.Channel(phase * 2 + 7)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 8)[0]))
#
## Calculate symetrical sequence components from phase voltages and currents (the components for phase 2 and 3 must be obtained by using the properties of the positive-, negative-, and zero-sequence sets which they represent (Bergen1999, pag 448))
#a = np.exp(1j*math.radians(120.0))
#Ainv = (1/3) * np.array(([1,1,1],[1,a,a**2],[1,a**2,a]),dtype=np.complex_)
#V_MV_Seq = np.dot(Ainv,V_MVpp)
#I_MV_Seq = np.dot(Ainv,I_MVpp)
#
#V_MV_s1 = np.array((V_MV_Seq[1], a**2*V_MV_Seq[1], a*V_MV_Seq[1]),dtype=np.complex_)
#V_MV_s2 = np.array((V_MV_Seq[2], a*V_MV_Seq[2], a**2*V_MV_Seq[2]),dtype=np.complex_)
#I_MV_s1 = np.array((I_MV_Seq[1], a**2*I_MV_Seq[1], a*I_MV_Seq[1]),dtype=np.complex_)
#I_MV_s2 = np.array((I_MV_Seq[2], a*I_MV_Seq[2], a**2*I_MV_Seq[2]),dtype=np.complex_)
#
#Qpos = np.abs(I_MV_s1) * np.abs(V_MV_s1) * np.sin(np.angle(V_MV_s1) - np.angle(I_MV_s1))
#Qneg = np.abs(V_MV_s2) * np.abs(I_MV_s2) * np.sin(np.angle(V_MV_s2) - np.angle(I_MV_s2))
#
#delta = np.abs(V_MV_s2[0])/np.abs(V_MV_s1[0])
#Qpos_comp = (Qpos - Qneg)/(1 - delta**2)
#Icomp_neg = I_MV_s2 + V_MV_s2/(1j*(np.absolute(V_MV_s1)**2/Qpos_comp))
#Qnegpos_comp = 2*np.abs(V_MV_s1) * np.abs(Icomp_neg) * np.sin(np.angle(V_MV_s1) - np.angle(Icomp_neg))        
#Qt = Qpos_comp/1000 + Qnegpos_comp/1000


# Normalization of S, per phase
Sc2_nom = Smax_auxiliary
Sc2 = np.sqrt(Pc2**2 + Qc2**2)
Sc2_pu = Sc2/Sc2_nom
PFc2 = abs(Pc2/Sc2)

# Reactive power compensation algorithm
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

# Active power in the primary side
Pc1 = Pc2/1.0
#Qc1 = [0,0,0]

print('Compensation reactive power: ' + str("{:.2f}".format(Qc1[0])) + ' | '+ str("{:.2f}".format(Qc1[1])) + ' | '+ str("{:.2f}".format(Qc1[2])) + '\n')

# Update STloads values
for phase in range(3):
    DSSText.Command = 'Edit Load.STload' + str(phase+1) + ' kW=' + str(Pc1[phase]) + ' kvar=' + str(Qc1[phase])


DSSText.Command ='Reset Monitors'
DSSSolution.Solve()
if not(DSSSolution.Converged):
    raise ValueError('Solution did not Converge')        
DSSMonitors.SampleAll() 
DSSMonitors.SaveAll()
DSSText.Command = 'Export Yprim Yprims.csv'

#for phase in range(3):
#    DSSMonitors.Name = 'Trafo' + str(phase+1) + '_prim-VI'
#    V_MV_delta[phase] = DSSMonitors.Channel(1)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(2)[0]))
#    I_MV_delta[phase] = DSSMonitors.Channel(5)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(6)[0]))
#
#DSSMonitors.Name = 'Trafo_prim-VI'
#for phase in range(3):
#    I_MV_line[phase] = DSSMonitors.Channel(phase * 2 + 7)[0] * np.exp(1j*math.radians(DSSMonitors.Channel(phase * 2 + 8)[0]))


print ('Powers:\t\tP [kW]\t\tQ [kVAr]\tS [kVA]')
for phase in range(3):
    print('\tPhase ' + str(phase+1))
    DSSMonitors.Name = 'Trafo_prim-PQ'
    print('Prim Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))
    DSSMonitors.Name = 'Trafo_sec-PQ'
    print('Sec Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))
    DSSMonitors.Name = 'Trafo_aux-PQ'
    print('Aux Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))) + '\n')

print()
DSSMonitors.Name = 'STvsources-PQ'
for phase in range(3):
    print('STvsources\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))
print()
DSSMonitors.Name = 'STloads-PQ'
for phase in range(3):
    print('STload\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(phase * 2 + 1)[0]**2 + DSSMonitors.Channel(phase * 2 + 2)[0]**2))))


print()
for phase in range(3):
    DSSMonitors.Name = 'LumpedLOAD' + str(phase+1) + '-PQ'
    print('Load' + str(phase+1) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(1)[0])) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(2)[0])) + '\t\t' + str("{:.2f}".format(math.sqrt(DSSMonitors.Channel(1)[0]**2 + DSSMonitors.Channel(2)[0]**2))))

for phase in range(3):
    print('STload\t\t' + str("{:.2f}".format(Pc1[phase])) + '\t\t' + str("{:.2f}".format(Qc1[phase])) + '\t\t' + str("{:.2f}".format(math.sqrt(Pc1[phase]**2 + Qc1[phase]**2))))

print('\nVoltages:')
for phase in range(3):
    print('\tPhase ' + str(phase+1))
    DSSMonitors.Name = 'Trafo_prim-VI'
    print('Prim Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('Prim Wdg (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg')
    DSSMonitors.Name = 'Trafo_sec-VI'
    print('Sec Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('Sec Wdg (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg')
    DSSMonitors.Name = 'Trafo_aux-VI'
    print('Aux Wdg \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg' + '\n')

    DSSMonitors.Name = 'STvsources-VI'
    print('STvs' + str(phase+1) + '\t\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('STvs' + str(phase+1) + ' (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg' + '\n')

    DSSMonitors.Name = 'STloads-VI'
    print('STload' + str(phase+1) + ' \t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 2)[0])) + ' deg')
    print('STload' + str(phase+1) + ' (I)\t' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 7)[0])) + ' A, ' + str("{:.2f}".format(DSSMonitors.Channel(phase * 2 + 8)[0])) + ' deg' + '\n')

for phase in range(3):
    DSSMonitors.Name = 'LumpedLOAD' + str(phase+1) + '-VI'
    print('\tLoad' + str(phase+1) + '\t' + str("{:.2f}".format(DSSMonitors.Channel(1)[0])) + ' V, ' + str("{:.2f}".format(DSSMonitors.Channel(2)[0])) + ' deg')
