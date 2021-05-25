#import time

def Extract_DSSmonitors_data(np,ElectricVars,DSSMonitors,Lines_set):
    
    ## Monitor channels ##
    # VI monitors: V1 ang1 V2 ang2 V3 ang3 I1 angI1 I2 angI2 I3 angI3 [V & deg]
    # PQ monitors: P1 Q1 P2 Q2 P3 Q3 [kW & kVAr]
    
    ## ElectricVars order ##
    # 0 <- Vmag_send
    # 1 <- Vang_send
    # 2 <- Vmag_rec
    # 3 <- Vang_rec
    # 4 <- Imag
    # 5 <- P_rec
    # 6 <- Q_rec
    # 7 <- P_send
    # 8 <- Q_send

    data = np.zeros((ElectricVars.size,Lines_set.size,3))
    
    # ALL LINES - sending and receiving buses
    for i_line in range(Lines_set.size):

        # Voltages and currents at receiving nodes
        DSSMonitors.Name = Lines_set[i_line] + '_VI_sending'
        for phase in range(3):
            data[0,i_line,phase] = DSSMonitors.Channel(phase * 2 + 1)[0] # Vmag_send
            data[1,i_line,phase] = DSSMonitors.Channel(phase * 2 + 2)[0] # Vang_send
            data[4,i_line,phase] = DSSMonitors.Channel(phase * 2 + 7)[0] # Imag

        # Voltages at receiving nodes
        DSSMonitors.Name = Lines_set[i_line] + '_VI_receiving'
        for phase in range(3):
            data[2,i_line,phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            data[3,i_line,phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
        
        # Active and reactive powers at receiving nodes
        DSSMonitors.Name = Lines_set[i_line] + '_PQ_receiving'
        for phase in range(3):
            data[5,i_line,phase] = -DSSMonitors.Channel(phase * 2 + 1)[0]
            data[6,i_line,phase] = -DSSMonitors.Channel(phase * 2 + 2)[0]

        # Active and reactive powers at sending nodes
        DSSMonitors.Name = Lines_set[i_line] + '_PQ_sending'
        for phase in range(3):
            data[7,i_line,phase] = DSSMonitors.Channel(phase * 2 + 1)[0]
            data[8,i_line,phase] = DSSMonitors.Channel(phase * 2 + 2)[0]
            
    return data
